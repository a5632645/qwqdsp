#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>
#include <span>

//------------------------------------------------------------------------------
// 双声道 FDN：正交反馈矩阵 + 每条延迟线一个一阶高搁架反馈滤波器。
//
// 参数化与 gdalsanto/diff-fdn-colorless（MIT License，见 license/fdn_colorless）一致，
// 但参数固定：延迟长度、反馈矩阵、输入/输出增益由该工程的 CLI 优化器（solver.py，
// gain_per_sample = 0.9999）离线求解，作为预设 struct 注入模板参数，运行期不可改。
// 默认 RT60 1439ms 即对应 γ = 0.9999（48kHz）。
//
// 该工程的参数划分与本实现的对应关系：
//   feedbackMatrix -> TPreset::kFeedbackMatrix
//   inputGain      -> TPreset::kInputGain     （该工程为单声道 N×1）
//   outputGain     -> TPreset::kOutputGain    （该工程为单声道 1×N）
//   Γ（每条线一个衰减滤波器）-> SetDecay / SetHighDecay + SetShelfFrequency
//
// 与该工程的差异：
//   - 该工程用 --delays 传入延迟线长度；此处固定为素数组合，作为预设的一部分。
//   - 该工程的 Γ 是 DecayFitNet + GEQ 设计的（1 增益 + 8 峰值 + 2 个二阶搁架）；
//     此处只用**一个一阶高搁架**：Γ 即两段常数衰减，搁架转折频率以下 T60 = SetDecay、
//     以上 T60 = SetHighDecay，每条线的总增益按延迟长度取幂以保证各线衰减一致。
//   - 该工程为单通道（B 为 N×1、C 为 1×N）；此处双声道：左声道按 (B 注入, C 读出)、
//     右声道按 (C 注入, B 读出)，两个向量都取自优化器，因此左右尾音不同，
//     且 L→L 通路的传递函数就是优化器的 C (D − VΓ)^{-1} B。
//     两个读出向量的模长不同（左右音量会差十几 dB），故各配一个归一化系数
//     kGainLeft / kGainRight，按湿信号能量把左右两路对齐、并把三档也对齐。
//------------------------------------------------------------------------------

namespace qwqdsp_fx {

namespace detail::fdn {

/// @brief 返回不小于 num 的素数（把按采样率缩放的延迟长度落回素数，保持互质）
inline size_t NextPrime(size_t num) noexcept {
    auto is_prime = [](size_t v) noexcept {
        if (v < 2) {
            return false;
        }
        if (v % 2 == 0) {
            return v == 2;
        }
        for (size_t d = 3; d * d <= v; d += 2) {
            if (v % d == 0) {
                return false;
            }
        }
        return true;
    };
    while (!is_prime(num)) {
        ++num;
    }
    return num;
}

// ------------------------------------------------------------
// Delay — 精确长度的循环延迟线
// ------------------------------------------------------------
/// @brief 精确样本长度的循环延迟线，缓冲区不补齐到 2 的幂
class Delay {
public:
    /// @brief 分配长度为 num_samples 的缓冲区
    /// @param num_samples 延迟长度（样本），不足 1 按 1 处理
    void Init(size_t num_samples) {
        size_ = std::max<size_t>(num_samples, 1);
        buffer_.reset(new float[size_]);
        std::fill_n(buffer_.get(), size_, 0.0f);
        write_idx_ = 0;
        read_idx_ = 0; // 读写指针同步前进且读先于写，读出的是 size_ 个样本前的值
    }

    /// @brief 读出延迟了 num_samples 的样本
    float Read() const noexcept {
        return buffer_[read_idx_];
    }

    /// @brief 写入当前样本
    void Push(float x) noexcept {
        buffer_[write_idx_] = x;
        if (++write_idx_ == size_) {
            write_idx_ = 0;
        }
        if (++read_idx_ == size_) {
            read_idx_ = 0;
        }
    }

    void Reset() noexcept {
        std::fill_n(buffer_.get(), size_, 0.0f);
        write_idx_ = 0;
        read_idx_ = 0;
    }

    /// @brief 延迟长度（样本）
    size_t GetDelay() const noexcept {
        return size_;
    }
private:
    std::unique_ptr<float[]> buffer_;
    size_t size_ = 1;
    size_t write_idx_ = 0;
    size_t read_idx_ = 0;
};

// ------------------------------------------------------------
// HighShelf — 一阶高搁架滤波器
// ------------------------------------------------------------
/// @brief 一阶高搁架滤波器，直流增益 1、Nyquist 增益 gain，转折频率 freq_hz
class HighShelf {
public:
    /// @param freq_hz 转折频率（Hz）
    /// @param sample_rate 采样率（Hz）
    /// @param gain Nyquist 处的线性增益，须为正数
    void Set(float freq_hz, float sample_rate, float gain) noexcept {
        float const t = std::tan(std::numbers::pi_v<float> * freq_hz / sample_rate);
        float const inv = 1.0f / (1.0f + t);
        b0_ = (gain + t) * inv;
        b1_ = (t - gain) * inv;
        a1_ = (t - 1.0f) * inv;
    }

    float Process(float x) noexcept {
        float const y = b0_ * x + z1_;
        z1_ = b1_ * x - a1_ * y;
        return y;
    }

    void Reset() noexcept {
        z1_ = 0.0f;
    }
private:
    float b0_ = 1.0f;
    float b1_ = 0.0f;
    float a1_ = 0.0f;
    float z1_ = 0.0f;
};

} // namespace detail::fdn

//------------------------------------------------------------------------------
// FdnReverb
//
// 每样本：
//   y_i = 第 i 条延迟线的输出（FDN 状态）
//   Γ_i(y_i) = line_gain_i · 高搁架_i(y_i)
//   延迟线输入_i = Σ_j V[i][j] · Γ_j(y_j) + Σ_c B[i][c] · x_c
//   输出_c = Σ_i C[c][i] · y_i
//
// 与参考实现的传递函数 H(z) = C (D(z) − V·Γ)^{-1} B 一致（D = diag(z^{m_i})）。
//------------------------------------------------------------------------------

/// @brief FDN 预设须提供延迟长度（48kHz 样本数）、优化后的反馈矩阵、单声道输入/输出增益，
///        以及左右两个输出向量各自的归一化系数
template <typename T>
concept CFdnPreset = requires {
    T::kDelays;
    T::kFeedbackMatrix;
    T::kInputGain;
    T::kOutputGain;
    T::kGainLeft;
    T::kGainRight;
};

/// @brief Fdn4PresetSmall：1129/1423/1783/2243 样本（24~47ms），延迟线合计 137ms / 26KB，
///        模态数 Σm=6578（模态密度 0.14/Hz，间隔 7.30Hz，模态重叠 M=0.21）
/// @note 反馈矩阵与输入/输出增益由 diff-fdn-colorless 的 CLI 优化器求得（γ = 0.9999）；
///       kGainLeft/kGainRight 分别归一化左右声道的输出向量（左读 C、右读 B），使左右音量一致并与其它预设对齐
struct Fdn4PresetSmall {
    static constexpr std::array<float, 4> kDelays = {1129.0f, 1423.0f, 1783.0f, 2243.0f};
    static constexpr std::array<float, 16> kFeedbackMatrix = {
        +0.491219f, -0.528313f, +0.530753f, -0.444849f,
        +0.521529f, +0.449300f, -0.486430f, -0.538073f,
        -0.484274f, +0.492863f, +0.506527f, -0.515747f,
        +0.502186f, +0.525453f, +0.474465f, +0.496580f,
    };
    static constexpr std::array<float, 4> kInputGain = {-0.130385f, -0.535105f, +0.166409f, -0.636976f};
    static constexpr std::array<float, 4> kOutputGain = {-0.734588f, -0.621857f, -0.804314f, -0.561190f};
    static constexpr float kGainLeft = 1.0194f;
    static constexpr float kGainRight = 1.6483f;
};

/// @brief Fdn4PresetMedium：1499/1889/2381/2999 样本（31~62ms），延迟线合计 183ms / 34KB，
///        模态数 Σm=8768（模态密度 0.18/Hz，间隔 5.47Hz，模态重叠 M=0.28）
/// @note 反馈矩阵与输入/输出增益由 diff-fdn-colorless 的 CLI 优化器求得（γ = 0.9999）；
///       kGainLeft/kGainRight 分别归一化左右声道的输出向量（左读 C、右读 B），使左右音量一致并与其它预设对齐
struct Fdn4PresetMedium {
    static constexpr std::array<float, 4> kDelays = {1499.0f, 1889.0f, 2381.0f, 2999.0f};
    static constexpr std::array<float, 16> kFeedbackMatrix = {
        +0.654276f, -0.660835f, -0.265084f, +0.254853f,
        +0.238213f, +0.623663f, -0.683896f, +0.294254f,
        +0.251043f, +0.263978f, +0.625164f, +0.690263f,
        -0.672425f, -0.323508f, -0.266808f, +0.609919f,
    };
    static constexpr std::array<float, 4> kInputGain = {-0.387659f, -1.359821f, -0.092109f, -0.415752f};
    static constexpr std::array<float, 4> kOutputGain = {+0.334917f, -0.329772f, -0.463563f, -0.617267f};
    static constexpr float kGainLeft = 1.5964f;
    static constexpr float kGainRight = 0.7455f;
};

/// @brief Fdn4PresetLarge：2011/2539/3203/4021 样本（42~84ms），延迟线合计 245ms / 46KB，
///        模态数 Σm=11774（模态密度 0.25/Hz，间隔 4.08Hz，模态重叠 M=0.37）
/// @note 反馈矩阵与输入/输出增益由 diff-fdn-colorless 的 CLI 优化器求得（γ = 0.9999）；
///       kGainLeft/kGainRight 分别归一化左右声道的输出向量（左读 C、右读 B），使左右音量一致并与其它预设对齐
struct Fdn4PresetLarge {
    static constexpr std::array<float, 4> kDelays = {2011.0f, 2539.0f, 3203.0f, 4021.0f};
    static constexpr std::array<float, 16> kFeedbackMatrix = {
        +0.655825f, +0.336935f, +0.291682f, +0.609336f,
        +0.349429f, +0.657733f, -0.293187f, -0.599440f,
        -0.578283f, +0.585028f, +0.567980f, +0.027023f,
        -0.336723f, +0.334073f, -0.711591f, +0.518316f,
    };
    static constexpr std::array<float, 4> kInputGain = {+0.028872f, -0.128786f, -1.672240f, -0.603617f};
    static constexpr std::array<float, 4> kOutputGain = {+0.629581f, -0.634752f, -0.288931f, +0.405541f};
    static constexpr float kGainLeft = 1.4481f;
    static constexpr float kGainRight = 0.5633f;
};

//==============================================================================
/** 反馈延迟网络混响器，延迟线数量与参数由预设给出。
 */
template <CFdnPreset TPreset>
class FdnReverb {
public:
    static constexpr size_t kNumLines = TPreset::kDelays.size();
    static constexpr size_t kNumChannels = 2;

    // 素数延迟长度以 48kHz 下的样本数给出，其它采样率按比例缩放后取最近的素数
    static constexpr float kPresetSampleRate = 48000.0f;

    static constexpr float kDefaultDecayMs = 1439.0f;
    static constexpr float kDefaultHighDecayMs = 900.0f;
    static constexpr float kDefaultShelf = 6000.0f;
    static constexpr float kDefaultMix = 0.35f;

    /// @brief 初始化采样率并分配延迟线
    /// @param sample_rate 采样率（Hz）
    void Init(float sample_rate) {
        sample_rate_ = sample_rate;
        float const scale = sample_rate_ / kPresetSampleRate;
        for (size_t i = 0; i < kNumLines; ++i) {
            size_t const samples = detail::fdn::NextPrime(static_cast<size_t>(std::lround(TPreset::kDelays[i] * scale)));
            lines_[i].Init(samples);
            delay_samples_[i] = static_cast<float>(samples);
        }
        UpdateFeedbackFilters();
        Reset();
    }

    /// @brief 设置搁架转折频率以下的 RT60
    /// @param rt60_ms RT60 时间（ms），范围 [1, 100000]
    void SetDecay(float rt60_ms) noexcept {
        decay_ms_ = std::clamp(rt60_ms, 1.0f, 100000.0f);
        UpdateFeedbackFilters();
    }

    /// @brief 设置搁架转折频率以上的 RT60
    /// @param rt60_ms RT60 时间（ms），范围 [1, 100000]
    void SetHighDecay(float rt60_ms) noexcept {
        high_decay_ms_ = std::clamp(rt60_ms, 1.0f, 100000.0f);
        UpdateFeedbackFilters();
    }

    /// @brief 设置高搁架转折频率
    /// @param freq_hz 转折频率（Hz），范围 [16, 20000]
    void SetShelfFrequency(float freq_hz) noexcept {
        shelf_hz_ = std::clamp(freq_hz, 16.0f, 20000.0f);
        UpdateFeedbackFilters();
    }

    /// @brief 设置干湿比
    /// @param m 干湿比 [0, 1]
    void SetMix(float m) noexcept {
        mix_ = std::clamp(m, 0.0f, 1.0f);
    }

    /// @brief 处理一对立体声样本
    /// @param dry_left 干信号左声道
    /// @param dry_right 干信号右声道
    /// @param[out] left_out 处理后左声道
    /// @param[out] right_out 处理后右声道
    void Process(float dry_left, float dry_right, float* left_out, float* right_out) noexcept {
        float const x[kNumChannels] = {dry_left, dry_right};

        // 1) FDN 状态 = 延迟线输出
        float y[kNumLines];
        for (size_t i = 0; i < kNumLines; ++i) {
            y[i] = lines_[i].Read();
        }

        // 2) 每条线上的衰减滤波器（一阶高搁架）
        float filtered[kNumLines];
        for (size_t i = 0; i < kNumLines; ++i) {
            filtered[i] = line_gain_[i] * shelf_[i].Process(y[i]);
        }

        // 3) 反馈矩阵混合 + 输入注入
        for (size_t i = 0; i < kNumLines; ++i) {
            float acc = 0.0f;
            for (size_t j = 0; j < kNumLines; ++j) {
                acc += TPreset::kFeedbackMatrix[i * kNumLines + j] * filtered[j];
            }
            for (size_t c = 0; c < kNumChannels; ++c) {
                acc += kInputGain[i * kNumChannels + c] * x[c];
            }
            lines_[i].Push(acc);
        }

        // 4) 输出抽头
        float wet[kNumChannels] = {0.0f, 0.0f};
        for (size_t c = 0; c < kNumChannels; ++c) {
            for (size_t i = 0; i < kNumLines; ++i) {
                wet[c] += kOutputGain[c * kNumLines + i] * y[i];
            }
        }

        *left_out = dry_left * (1 - mix_) + wet[0] * mix_;
        *right_out = dry_right * (1 - mix_) + wet[1] * mix_;
    }

    /// @brief 批量处理立体声样本（in-place）
    /// @param l 左声道输入/输出缓冲
    /// @param r 右声道输入/输出缓冲
    /// @param num 样本帧数
    void Process(float* l, float* r, size_t num) noexcept {
        for (size_t i = 0; i < num; ++i) {
            Process(l[i], r[i], &l[i], &r[i]);
        }
    }

    /// @brief 复位所有延迟线与滤波器状态
    void Reset() noexcept {
        for (auto& line : lines_) {
            line.Reset();
        }
        for (auto& shelf : shelf_) {
            shelf.Reset();
        }
    }

    /// @brief 各延迟线的实际长度（样本数，已按采样率缩放并吸附到素数）
    std::span<const float, kNumLines> GetDelaySamples() const noexcept {
        return delay_samples_;
    }
private:
    static_assert(TPreset::kFeedbackMatrix.size() == kNumLines * kNumLines);
    static_assert(TPreset::kInputGain.size() == kNumLines);
    static_assert(TPreset::kOutputGain.size() == kNumLines);

    // 搁架的高频增益是相对低频段的比值，限制幅度以免极端 T60 下系数发散
    static constexpr float kMinShelfGain = 0.05f;
    static constexpr float kMaxShelfGain = 20.0f;

    /// @brief 双声道注入增益：左声道用优化器的 B、右声道用优化器的 C
    static constexpr std::array<float, kNumLines * kNumChannels> MakeInputGain() noexcept {
        std::array<float, kNumLines * kNumChannels> gain{};
        for (size_t i = 0; i < kNumLines; ++i) {
            gain[i * kNumChannels] = TPreset::kInputGain[i];
            gain[i * kNumChannels + 1] = TPreset::kOutputGain[i];
        }
        return gain;
    }

    /// @brief 双声道输出增益：左声道读优化器的 C、右声道读优化器的 B，各自乘该档的归一化系数
    static constexpr std::array<float, kNumChannels * kNumLines> MakeOutputGain() noexcept {
        std::array<float, kNumChannels * kNumLines> gain{};
        for (size_t i = 0; i < kNumLines; ++i) {
            gain[i] = TPreset::kOutputGain[i] * TPreset::kGainLeft;
            gain[kNumLines + i] = TPreset::kInputGain[i] * TPreset::kGainRight;
        }
        return gain;
    }

    static constexpr std::array<float, kNumLines * kNumChannels> kInputGain = MakeInputGain();
    static constexpr std::array<float, kNumChannels * kNumLines> kOutputGain = MakeOutputGain();

    /// @brief 每采样点的衰减系数，使该频段 T60 后衰减 60dB
    float PerSampleGain(float rt60_ms) const noexcept {
        return std::pow(10.0f, -3.0f / (sample_rate_ * rt60_ms * 0.001f));
    }

    /// @brief 按当前 T60 与搁架频率重算每线的总增益与高搁架
    void UpdateFeedbackFilters() noexcept {
        float const g_ref = PerSampleGain(decay_ms_);
        float const g_high = PerSampleGain(high_decay_ms_);
        for (size_t i = 0; i < kNumLines; ++i) {
            float const m = delay_samples_[i];
            line_gain_[i] = std::pow(g_ref, m);
            float const high_gain = std::clamp(std::pow(g_high / g_ref, m), kMinShelfGain, kMaxShelfGain);
            shelf_[i].Set(shelf_hz_, sample_rate_, high_gain);
        }
    }

    float sample_rate_ = 1.0f;
    float mix_ = 0.0f;
    float decay_ms_ = kDefaultDecayMs;
    float high_decay_ms_ = kDefaultHighDecayMs;
    float shelf_hz_ = kDefaultShelf;

    std::array<detail::fdn::Delay, kNumLines> lines_;
    std::array<detail::fdn::HighShelf, kNumLines> shelf_;
    std::array<float, kNumLines> delay_samples_ = {};
    std::array<float, kNumLines> line_gain_ = {};
};

} // namespace qwqdsp_fx
