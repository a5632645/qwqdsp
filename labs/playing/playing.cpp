#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iostream>
#include <numbers>
#include <string>

#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/blackman_harris_3term.hpp>

namespace playing {

static constexpr float kSampleRate = 48000.0f;
static constexpr int kWindowWidth = 900;
static constexpr int kWindowHeight = 480;
static constexpr const char* kWindowTitle = "Playing - PGHI";

// ------------------------------------------------------------
// TransientMode
// ------------------------------------------------------------

/**
 * @brief 瞬态检测算法选择
 *
 * 移调核心始终为 PGHI 相位传播，检测器只在检测到瞬态时重置合成相位。
 */
enum class TransientMode {
    None,      ///< 无瞬态检测，纯 PGHI
    Flux,      ///< 线性频谱通量 + 中位数基线自适应阈值
    SuperFlux, ///< 对数频带三角滤波器组 + 频率最大滤波 + 固定阈值
    Vocoder,   ///< log1p 谱通量 / 能量归一化 + EMA 均值方差双阈值
    DsPark,    ///< DSPark 谱通量：对数滤波器组 + 运行中位数 + 固定增量
};

static constexpr std::array<const char*, 5> kTransientModeNames{"None", "Flux", "SuperFlux", "Vocoder", "DSPark"};

// ------------------------------------------------------------
// RealtimePitchShifter
// ------------------------------------------------------------

namespace detail {

// 周期三项 Blackman-Harris 窗，六倍重叠。
constexpr std::size_t kOverlap = 6;
// 六倍重叠下周期三项 Blackman-Harris 窗平方和为常数。
constexpr float kOlaGain = 0.5431781f;
using WindowType = qwqdsp_window::BlackmanHarrisThreeTerm;

constexpr std::size_t nextPowerOfTwo(std::size_t value) noexcept {
    std::size_t result = 1;
    while (result < value) {
        result <<= 1;
    }
    return result;
}

static inline float wrapToPi(float phase) noexcept {
    constexpr float kTwoPi = 2.0f * std::numbers::pi_v<float>;
    return phase - kTwoPi * std::round(phase / kTwoPi);
}

// sinc 重采样每侧抽头数（共 2 * kSincRadius 个）
constexpr std::size_t kSincRadius = 8;
// sinc 核查找表长度（覆盖 [0, kSincRadius]）
constexpr std::size_t kSincTableSize = 2048;

/**
 * @brief 归一化 sinc 乘对称 Blackman 窗的内插核
 * @param[in] x 距插值中心的样本距离，单位采样
 */
static inline float sincKernel(float x) noexcept {
    if (std::abs(x) < 1.0e-6f) {
        return 1.0f;
    }
    const float px = std::numbers::pi_v<float> * x;
    const float window = 0.42f + 0.5f * std::cos(std::numbers::pi_v<float> * x / static_cast<float>(kSincRadius))
                       + 0.08f * std::cos(2.0f * std::numbers::pi_v<float> * x / static_cast<float>(kSincRadius));
    return std::sin(px) / px * window;
}

} // namespace detail

/**
 * @brief 使用固定容量缓冲的实时 PGHI 移调器
 *
 * 处理器采用 8192 点 STFT（窗口 3072，六倍重叠）。可变分析 hop 完成时间
 * 拉伸，固定合成 hop 保证 COLA，随后对连续 OLA 信号用窗函数 sinc 插值
 * 重采样完成移调，并静音移调后超出奈奎斯特频率的 bin 防止混叠。
 * 瞬态检测器（Flux / SuperFlux / TransientVocoder）仅在检测到瞬态时重置合成相位，
 * 其余帧统一使用 PGHI 相位传播。初始化完成后音频处理路径不进行动态内存分配。
 */
class RealtimePitchShifter {
public:
    static constexpr std::size_t kHopSize = 512;
    static constexpr std::size_t kWindowSize = kHopSize * detail::kOverlap;
    static constexpr std::size_t kFrameSize = detail::nextPowerOfTwo(kWindowSize * 2);
    static constexpr std::size_t kNumBins = kFrameSize / 2 + 1;
    static constexpr std::size_t kMaxTransientMarks = 32;
    // 波形图与瞬态检测点的时间范围，单位秒
    static constexpr float kDisplayTimeSeconds = 3.0f;
    // 波形显示缓冲长度（kDisplayTimeSeconds 秒采样数）
    static constexpr std::size_t kDisplaySize = static_cast<std::size_t>(kSampleRate * kDisplayTimeSeconds);

    /**
     * @brief 初始化 FFT、窗函数和处理状态
     * @param[in] sample_rate 输入输出采样率，单位 Hz
     */
    void init(float sample_rate) {
        sample_rate_ = sample_rate;
        fft_.Init(kFrameSize);
        detail::WindowType::Window(window_, true);
        buildSuperFluxFilterBank();
        buildSincTable();
        reset();
    }

    /**
     * @brief 设置音高移动量
     * @param[in] semitones 半音数，正值升高音高，负值降低音高
     */
    void setPitchShift(float semitones) {
        target_ratio_ = std::exp2(std::clamp(semitones, -12.0f, 12.0f) / 12.0f);
        if (first_frame_) {
            current_ratio_ = target_ratio_;
        }
    }

    /**
     * @brief 切换瞬态检测算法
     * @param[in] mode 要使用的瞬态检测算法
     */
    void setTransientMode(TransientMode mode) {
        if (transient_mode_ == mode) {
            return;
        }
        transient_mode_ = mode;
        resetPhaseState();
    }

    /**
     * @brief 清空流式缓冲和全部算法状态
     */
    void reset() {
        input_ring_.fill(0.0f);
        synthesis_ring_.fill(0.0f);
        input_write_ = 0;
        collected_samples_ = 0;
        samples_since_frame_ = 0;
        samples_to_frame_ = 0.0f;
        analysis_started_ = false;
        synthesis_frame_start_ = 0;
        synthesis_available_ = 0;
        synthesis_cleared_until_ = 0;
        resample_position_ = 0.0;
        resampler_started_ = false;
        current_ratio_ = target_ratio_;
        input_position_ = 0;
        display_position_.store(0, std::memory_order_relaxed);
        display_ring_.fill(0.0f);
        display_write_ = 0;
        for (auto& mark : transient_marks_) {
            mark.store(0, std::memory_order_relaxed);
        }
        transient_mark_write_.store(0, std::memory_order_relaxed);
        resetPhaseState();
    }

    /**
     * @brief 处理一段单声道音频
     * @param[in] input 单声道浮点输入缓冲
     * @param[out] output 单声道浮点输出缓冲
     * @param[in] frame_count 缓冲中的采样数
     */
    void process(const float* input, float* output, std::size_t frame_count) {
        if (output == nullptr) {
            return;
        }
        if (input == nullptr) {
            std::fill_n(output, frame_count, 0.0f);
            return;
        }

        for (std::size_t i = 0; i < frame_count; ++i) {
            input_ring_[input_write_] = input[i];
            input_write_ = (input_write_ + 1) % kWindowSize;
            collected_samples_ = std::min(collected_samples_ + 1, kWindowSize);
            ++samples_since_frame_;
            ++input_position_;
            display_position_.store(input_position_, std::memory_order_relaxed);
            display_ring_[display_write_] = input[i];
            display_write_ = (display_write_ + 1) % kDisplaySize;

            if (!analysis_started_ && collected_samples_ == kWindowSize) {
                processFrame(kHopSize);
                analysis_started_ = true;
                samples_since_frame_ = 0;
                samples_to_frame_ = static_cast<float>(kHopSize) / current_ratio_;
            }
            else if (analysis_started_) {
                samples_to_frame_ -= 1.0f;
                if (samples_to_frame_ <= 0.0f) {
                    processFrame(std::max<std::size_t>(1, samples_since_frame_));
                    samples_since_frame_ = 0;
                    samples_to_frame_ += static_cast<float>(kHopSize) / current_ratio_;
                }
            }
            output[i] = renderOutputSample();
        }
    }

    /**
     * @brief 返回最新输入采样在整段信号中的绝对位置
     */
    std::uint64_t displayPosition() const noexcept {
        return display_position_.load(std::memory_order_relaxed);
    }

    /**
     * @brief 读取指定绝对位置的输入采样（供 GUI 绘制 1 秒波形）
     * @param[in] absolute_position 输入采样绝对位置
     */
    float displaySampleAt(std::uint64_t absolute_position) const noexcept {
        return display_ring_[absolute_position % kDisplaySize];
    }

    /**
     * @brief 读取第 index 个瞬态检测点的采样绝对位置
     */
    std::uint64_t transientMark(std::size_t index) const noexcept {
        return transient_marks_[index % kMaxTransientMarks].load(std::memory_order_relaxed);
    }

    /**
     * @brief 返回已记录的瞬态检测点总数
     */
    std::uint32_t transientMarkCount() const noexcept {
        return transient_mark_write_.load(std::memory_order_relaxed);
    }
private:
    static constexpr std::size_t kSynthesisRingSize = kFrameSize * 16;
    static constexpr std::size_t kFluxHistorySize = 32;
    static constexpr std::size_t kHeapCapacity = kNumBins * 2;
    static constexpr std::size_t kMaxSuperFluxBands = 256;
    static constexpr std::size_t kMaxSuperFluxWeights = kNumBins * 3;
    // DSPark 通量历史窗口与触发增量（对应 PhaseVocoderEngine 的 kFluxWindow / kFluxDelta）
    static constexpr std::size_t kDsFluxWindow = 12;
    static constexpr float kDsFluxDelta = 0.02f;
    // SuperFlux / DSPark 触发后的冷却帧数（约 64ms @48kHz，覆盖攻击尾巴防止连续重置）
    static constexpr std::size_t kSuperFluxCooldown = 6;
    static constexpr std::size_t kDsParkCooldown = 6;

    struct HeapItem {
        float magnitude = 0.0f;
        std::size_t bin = 0;
        bool is_previous = false;
    };

    void resetPhaseState() {
        previous_analysis_phase_.fill(0.0f);
        previous_synthesis_phase_.fill(0.0f);
        previous_magnitude_.fill(0.0f);
        flux_history_.fill(0.0f);
        flux_history_count_ = 0;
        flux_history_write_ = 0;
        flux_cooldown_ = 0;
        transient_mean_ = 0.0f;
        transient_variance_ = 0.0f;
        transient_frames_ = 0;
        transient_cooldown_ = 0;
        transient_has_mean_ = false;
        superflux_current_.fill(0.0f);
        superflux_previous_max_.fill(0.0f);
        superflux_cooldown_ = 0;
        dsflux_current_.fill(0.0f);
        dsflux_max_prev_.fill(0.0f);
        dsflux_history_.fill(0.0f);
        dsflux_scratch_.fill(0.0f);
        dsflux_pos_ = 0;
        dsflux_filled_ = 0;
        dsflux_cooldown_ = 0;
        first_frame_ = true;
    }

    void recordTransientMark() noexcept {
        const std::uint32_t index = transient_mark_write_.load(std::memory_order_relaxed);
        transient_marks_[index % kMaxTransientMarks].store(input_position_, std::memory_order_relaxed);
        transient_mark_write_.store(index + 1, std::memory_order_relaxed);
    }

    void processFrame(std::size_t analysis_hop) {
        current_ratio_ += 0.25f * (target_ratio_ - current_ratio_);
        constexpr std::size_t kHalfWindow = kWindowSize / 2;

        // PGHI 采用居中的分半窗装入 FFT 输入
        fft_input_.fill(0.0f);
        for (std::size_t i = 0; i < kHalfWindow; ++i) {
            fft_input_[i] = input_ring_[(input_write_ + kHalfWindow + i) % kWindowSize] * window_[kHalfWindow + i];
        }
        const std::size_t pad = kFrameSize - kHalfWindow;
        for (std::size_t i = 0; i < kHalfWindow; ++i) {
            fft_input_[pad + i] = input_ring_[(input_write_ + i) % kWindowSize] * window_[i];
        }
        fft_.FFT(fft_input_.data(), fft_output_.data());

        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            current_re_[bin] = fft_output_[2 * bin];
            current_im_[bin] = fft_output_[2 * bin + 1];
            current_magnitude_[bin] = std::hypot(current_re_[bin], current_im_[bin]);
            current_analysis_phase_[bin] = std::atan2(current_im_[bin], current_re_[bin]);
        }

        bool reset_phase = first_frame_;
        if (!first_frame_) {
            if (transient_mode_ == TransientMode::Flux) {
                reset_phase = detectFluxTransient();
            }
            else if (transient_mode_ == TransientMode::SuperFlux) {
                reset_phase = detectSuperFluxTransient();
            }
            else if (transient_mode_ == TransientMode::Vocoder) {
                reset_phase = detectVocoderTransient();
            }
            else if (transient_mode_ == TransientMode::DsPark) {
                reset_phase = detectDsParkTransient();
            }
        }

        if (reset_phase) {
            current_synthesis_phase_ = current_analysis_phase_;
            if (!first_frame_) {
                // 记录瞬态检测点（帧结束处的输入采样绝对位置）
                recordTransientMark();
            }
        }
        else {
            propagatePghi(analysis_hop);
        }

        synthesizeSpectrum();
        fft_.IFFT(fft_output_.data(), fft_input_.data());
        for (std::size_t i = 0; i < kHalfWindow; ++i) {
            const std::size_t index = static_cast<std::size_t>((synthesis_frame_start_ + i) % kSynthesisRingSize);
            const std::size_t source = kFrameSize - kHalfWindow + i;
            synthesis_ring_[index] += fft_input_[source] * window_[i] * detail::kOlaGain;
        }
        for (std::size_t i = 0; i < kHalfWindow; ++i) {
            const std::size_t index =
                static_cast<std::size_t>((synthesis_frame_start_ + kHalfWindow + i) % kSynthesisRingSize);
            const std::size_t source = i;
            synthesis_ring_[index] += fft_input_[source] * window_[kHalfWindow + i] * detail::kOlaGain;
        }
        // 新合成帧在 hop 边界处的窗值为零，因此该边界样本已完整可读。
        synthesis_available_ = synthesis_frame_start_ + kHopSize + 1;
        synthesis_frame_start_ += kHopSize;

        previous_analysis_phase_ = current_analysis_phase_;
        previous_synthesis_phase_ = current_synthesis_phase_;
        previous_magnitude_ = current_magnitude_;
        first_frame_ = false;
    }

    /**
     * @brief PGHI 相位传播：幅度优先堆遍历，从高幅度 bin 向相邻 bin 扩展相位
     * @param[in] analysis_hop 本帧实际分析 hop
     */
    void propagatePghi(std::size_t analysis_hop) {
        processed_bins_.fill(0);
        heap_size_ = 0;
        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            heapPush({previous_magnitude_[bin], bin, true});
        }

        const float analysis_scale =
            2.0f * std::numbers::pi_v<float> * static_cast<float>(analysis_hop) / static_cast<float>(kFrameSize);
        constexpr float kSynthesisScale =
            2.0f * std::numbers::pi_v<float> * static_cast<float>(kHopSize) / static_cast<float>(kFrameSize);
        const float stretch_ratio = static_cast<float>(kHopSize) / static_cast<float>(analysis_hop);
        std::size_t processed_count = 0;
        while (processed_count < kNumBins && heap_size_ > 0) {
            const HeapItem item = heapPop();
            const std::size_t bin = item.bin;
            if (item.is_previous) {
                if (processed_bins_[bin] != 0) {
                    continue;
                }
                const float expected = analysis_scale * static_cast<float>(bin);
                const float residual =
                    detail::wrapToPi(current_analysis_phase_[bin] - previous_analysis_phase_[bin] - expected);
                current_synthesis_phase_[bin] =
                    detail::wrapToPi(previous_synthesis_phase_[bin] + kSynthesisScale * static_cast<float>(bin)
                                     + stretch_ratio * residual);
                processed_bins_[bin] = 1;
                ++processed_count;
                heapPush({current_magnitude_[bin], bin, false});
                continue;
            }

            if (bin + 1 < kNumBins && processed_bins_[bin + 1] == 0) {
                const float delta = detail::wrapToPi(current_analysis_phase_[bin + 1] - current_analysis_phase_[bin]);
                current_synthesis_phase_[bin + 1] =
                    detail::wrapToPi(current_synthesis_phase_[bin] + stretch_ratio * delta);
                processed_bins_[bin + 1] = 1;
                ++processed_count;
                heapPush({current_magnitude_[bin + 1], bin + 1, false});
            }
            if (bin > 0 && processed_bins_[bin - 1] == 0) {
                const float delta = detail::wrapToPi(current_analysis_phase_[bin - 1] - current_analysis_phase_[bin]);
                current_synthesis_phase_[bin - 1] =
                    detail::wrapToPi(current_synthesis_phase_[bin] + stretch_ratio * delta);
                processed_bins_[bin - 1] = 1;
                ++processed_count;
                heapPush({current_magnitude_[bin - 1], bin - 1, false});
            }
        }
    }

    /**
     * @brief 线性频谱通量瞬态检测
     *
     * 累加各 bin 幅度的正增量作为通量，与最近 32 帧非瞬态通量的中位数
     * 比较，超过 max(4 * baseline, 0.02) 判定为瞬态；触发后冷却 6 帧，
     * 瞬态帧不进入历史，防止异常峰值污染基线。
     */
    bool detectFluxTransient() {
        float positive_difference = 0.0f;
        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            positive_difference += std::max(0.0f, current_magnitude_[bin] - previous_magnitude_[bin]);
        }
        const float flux = positive_difference / static_cast<float>(kNumBins);

        float baseline = 0.0f;
        if (flux_history_count_ > 0) {
            flux_scratch_ = flux_history_;
            const auto middle = flux_scratch_.begin() + static_cast<std::ptrdiff_t>(flux_history_count_ / 2);
            std::nth_element(flux_scratch_.begin(), middle,
                             flux_scratch_.begin() + static_cast<std::ptrdiff_t>(flux_history_count_));
            baseline = *middle;
        }
        const bool transient = flux_cooldown_ == 0 && flux > std::max(4.0f * baseline, 0.02f);
        flux_cooldown_ = transient ? 6 : (flux_cooldown_ > 0 ? flux_cooldown_ - 1 : 0);
        if (!transient) {
            flux_history_[flux_history_write_] = flux;
            flux_history_write_ = (flux_history_write_ + 1) % kFluxHistorySize;
            flux_history_count_ = std::min(flux_history_count_ + 1, kFluxHistorySize);
        }
        return transient;
    }

    /**
     * @brief SuperFlux 瞬态检测
     *
     * 对数频带三角滤波器组（27.5Hz 起、每倍频程 24 带）聚合幅度后做
     * log10 压缩，与上一帧频率方向三邻域最大滤波值比较，只累加正差；
     * 平均正差超过 0.03 判定为瞬态。触发后冷却 6 帧防止攻击尾巴连续重置。
     */
    bool detectSuperFluxTransient() {
        if (superflux_band_count_ == 0) {
            return false;
        }

        float flux = 0.0f;
        constexpr float kReferenceFrameSize = 4096.0f;
        const float scale = kReferenceFrameSize / static_cast<float>(kFrameSize);
        for (std::size_t band = 0; band < superflux_band_count_; ++band) {
            float magnitude = 0.0f;
            for (std::size_t i = 0; i < superflux_band_sizes_[band]; ++i) {
                const std::size_t bin = superflux_band_starts_[band] + i;
                magnitude += current_magnitude_[bin] * superflux_weights_[superflux_band_offsets_[band] + i];
            }
            superflux_current_[band] = std::log10(magnitude * scale + 1.0f);
            if (!first_frame_) {
                flux += std::max(0.0f, superflux_current_[band] - superflux_previous_max_[band]);
            }
        }

        for (std::size_t band = 0; band < superflux_band_count_; ++band) {
            float maximum = superflux_current_[band];
            if (band > 0) {
                maximum = std::max(maximum, superflux_current_[band - 1]);
            }
            if (band + 1 < superflux_band_count_) {
                maximum = std::max(maximum, superflux_current_[band + 1]);
            }
            superflux_previous_max_[band] = maximum;
        }
        // 冷却期内不触发，但上一段的旋转基准保持每帧更新
        const bool transient = superflux_cooldown_ == 0 && !first_frame_
                            && flux / static_cast<float>(superflux_band_count_) > 0.03f;
        superflux_cooldown_ = transient ? kSuperFluxCooldown
                                        : (superflux_cooldown_ > 0 ? superflux_cooldown_ - 1 : 0);
        return transient;
    }

    /**
     * @brief DSPark 相位声码器的谱通量瞬态检测
     *
     * 特征与 SuperFlux 相同（对数频带三角滤波器组 + log10 压缩 + 上一帧
     * 三邻域最大滤波），判定改为：通量相对最近 12 帧通量历史中位数的增量
     * 达到固定门限 0.02 即触发。历史每帧推进；触发后冷却 6 帧防止连续重置。
     */
    bool detectDsParkTransient() {
        if (superflux_band_count_ == 0) {
            return false;
        }

        // 通量：与上一帧三邻域最大滤波值比较，只累加正差，除以频带数
        float flux = 0.0f;
        constexpr float kReferenceFrameSize = 2048.0f;
        const float scale = kReferenceFrameSize / static_cast<float>(kFrameSize);
        for (std::size_t band = 0; band < superflux_band_count_; ++band) {
            float magnitude = 0.0f;
            for (std::size_t i = 0; i < superflux_band_sizes_[band]; ++i) {
                const std::size_t bin = superflux_band_starts_[band] + i;
                magnitude += current_magnitude_[bin] * superflux_weights_[superflux_band_offsets_[band] + i];
            }
            dsflux_current_[band] = std::log10(magnitude * scale + 1.0f);
            if (!first_frame_) {
                flux += std::max(0.0f, dsflux_current_[band] - dsflux_max_prev_[band]);
            }
        }
        flux /= static_cast<float>(std::max<std::size_t>(1, superflux_band_count_));

        // 中位数基线取在当前帧推入历史之前，帧只与先前的材料比较
        bool rise = false;
        if (dsflux_filled_ >= kDsFluxWindow) {
            dsflux_scratch_ = dsflux_history_;
            const auto middle = dsflux_scratch_.begin() + static_cast<std::ptrdiff_t>(kDsFluxWindow / 2);
            std::nth_element(dsflux_scratch_.begin(), middle, dsflux_scratch_.end());
            rise = (flux - *middle) >= kDsFluxDelta;
        }

        // 历史每帧推进，无论是否触发
        dsflux_history_[dsflux_pos_] = flux;
        dsflux_pos_ = (dsflux_pos_ + 1) % kDsFluxWindow;
        dsflux_filled_ = std::min(dsflux_filled_ + 1, kDsFluxWindow);

        // 旋转：当前帧做三邻域最大滤波，作为下一帧的差分基准
        for (std::size_t band = 0; band < superflux_band_count_; ++band) {
            float maximum = dsflux_current_[band];
            if (band > 0) {
                maximum = std::max(maximum, dsflux_current_[band - 1]);
            }
            if (band + 1 < superflux_band_count_) {
                maximum = std::max(maximum, dsflux_current_[band + 1]);
            }
            dsflux_max_prev_[band] = maximum;
        }
        // 冷却期内不触发，历史与旋转基准保持每帧推进
        const bool transient = dsflux_cooldown_ == 0 && !first_frame_ && rise;
        dsflux_cooldown_ = transient ? kDsParkCooldown : (dsflux_cooldown_ > 0 ? dsflux_cooldown_ - 1 : 0);
        return transient;
    }

    /**
     * @brief TransientVocoder 瞬态检测
     *
     * 对 log1p 幅度正增量累加作为通量，除以频率加权的 log1p 能量归一化；
     * 累计超过 4 帧且无冷却时，同时满足
     *   norm_flux > mean + 1.5 * max(0.07, std) 且 norm_flux > mean * 1.35
     * 判定为瞬态。均值/方差用 EMA 更新，瞬态帧 alpha 更大以快速适应。
     */
    bool detectVocoderTransient() {
        float flux = 0.0f;
        float energy = 0.0f;
        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            const float weight = 0.5f + 0.5f * static_cast<float>(bin) / static_cast<float>(kNumBins - 1);
            const float current_log = std::log1p(current_magnitude_[bin]);
            flux += std::max(0.0f, current_log - std::log1p(previous_magnitude_[bin]));
            energy += weight * current_log;
        }
        const float normalized_flux = flux / std::max(energy, 1.0e-9f);
        const float deviation = std::sqrt(transient_variance_);
        const bool transient = transient_frames_ > 4 && transient_cooldown_ == 0
                            && normalized_flux > transient_mean_ + 1.5f * std::max(0.07f, deviation)
                            && normalized_flux > transient_mean_ * 1.35f;

        const float alpha = transient ? 0.3f : 0.12f;
        if (!transient_has_mean_) {
            transient_mean_ = normalized_flux;
            transient_has_mean_ = true;
        }
        else {
            const float delta = normalized_flux - transient_mean_;
            transient_mean_ += alpha * delta;
            transient_variance_ = (1.0f - alpha) * (transient_variance_ + alpha * delta * delta);
        }
        transient_cooldown_ = transient ? 1 : (transient_cooldown_ > 0 ? transient_cooldown_ - 1 : 0);
        ++transient_frames_;
        return transient;
    }

    /**
     * @brief 由幅度与合成相位重建频谱
     *
     * 移调（重采样）后超出奈奎斯特频率的 bin 会被折叠混叠回通带，这里按
     * 当前重采样比例静音掉这些 bin，避免高频混叠。
     */
    void synthesizeSpectrum() {
        fft_output_.fill(0.0f);
        const float shifted_nyquist = static_cast<float>(kNumBins - 1) / std::max(current_ratio_, 1.0e-6f);
        for (std::size_t bin = 0; bin < kNumBins; ++bin) {
            if (static_cast<float>(bin) > shifted_nyquist) {
                break;
            }
            fft_output_[2 * bin] = current_magnitude_[bin] * std::cos(current_synthesis_phase_[bin]);
            fft_output_[2 * bin + 1] = current_magnitude_[bin] * std::sin(current_synthesis_phase_[bin]);
        }
        fft_output_[1] = 0.0f;
        fft_output_[2 * (kNumBins - 1) + 1] = 0.0f;
    }

    /**
     * @brief 构建 SuperFlux 对数频带三角滤波器组
     *
     * 中心频率从 27.5Hz 起每倍频程 24 个，相邻三个中心构成一个三角带。
     */
    void buildSuperFluxFilterBank() {
        constexpr float kMinimumFrequency = 27.5f;
        constexpr float kMaximumFrequency = 16000.0f;
        constexpr int kBandsPerOctave = 24;
        std::array<std::size_t, kMaxSuperFluxBands + 2> centers{};
        std::size_t center_count = 0;
        const float bin_hz = sample_rate_ / static_cast<float>(kFrameSize);
        const float maximum_frequency = std::min(kMaximumFrequency, sample_rate_ * 0.5f * 0.999f);

        for (int index = 0; center_count < centers.size(); ++index) {
            const float frequency =
                kMinimumFrequency * std::exp2(static_cast<float>(index) / static_cast<float>(kBandsPerOctave));
            if (frequency > maximum_frequency) {
                break;
            }
            const std::size_t bin =
                std::min<std::size_t>(static_cast<std::size_t>(std::lround(frequency / bin_hz)), kNumBins - 1);
            if (center_count == 0 || bin > centers[center_count - 1]) {
                centers[center_count++] = bin;
            }
        }

        superflux_band_count_ = 0;
        superflux_weight_count_ = 0;
        for (std::size_t center = 1; center + 1 < center_count; ++center) {
            const std::size_t low = centers[center - 1];
            const std::size_t middle = centers[center];
            const std::size_t high = centers[center + 1];
            const std::size_t size = high - low + 1;
            if (low >= middle || middle >= high || superflux_band_count_ == kMaxSuperFluxBands
                || superflux_weight_count_ + size > kMaxSuperFluxWeights) {
                continue;
            }

            const std::size_t band = superflux_band_count_++;
            superflux_band_starts_[band] = low;
            superflux_band_sizes_[band] = size;
            superflux_band_offsets_[band] = superflux_weight_count_;
            for (std::size_t bin = low; bin <= high; ++bin) {
                const float weight = bin <= middle ? static_cast<float>(bin - low) / static_cast<float>(middle - low)
                                                   : static_cast<float>(high - bin) / static_cast<float>(high - middle);
                superflux_weights_[superflux_weight_count_++] = weight;
            }
        }
    }

    /**
     * @brief 预计算窗函数 sinc 内插核查找表
     */
    void buildSincTable() {
        constexpr float kStep = static_cast<float>(detail::kSincRadius) / static_cast<float>(detail::kSincTableSize);
        for (std::size_t i = 0; i < detail::kSincTableSize + 1; ++i) {
            sinc_table_[i] = detail::sincKernel(static_cast<float>(i) * kStep);
        }
    }

    /**
     * @brief 使用窗函数 sinc 插值重采样合成信号
     *
     * 每个输出采样以重采样位置为中心取 2 * kSincRadius 个抽头，与预计算的
     * Blackman 窗 sinc 核做内积；相比线性插值可显著降低镜频与高频衰减。
     */
    float renderOutputSample() {
        constexpr std::uint64_t kResampleStart = kFrameSize - kHopSize;
        if (!resampler_started_) {
            if (synthesis_available_ <= kResampleStart + 1) {
                return 0.0f;
            }
            resample_position_ = static_cast<double>(kResampleStart);
            while (synthesis_cleared_until_ < kResampleStart) {
                synthesis_ring_[static_cast<std::size_t>(synthesis_cleared_until_ % kSynthesisRingSize)] = 0.0f;
                ++synthesis_cleared_until_;
            }
            resampler_started_ = true;
        }

        const auto lower = static_cast<std::uint64_t>(resample_position_);
        if (lower + detail::kSincRadius >= synthesis_available_) {
            return 0.0f;
        }
        const float fraction = static_cast<float>(resample_position_ - static_cast<double>(lower));
        constexpr float kTableScale =
            static_cast<float>(detail::kSincTableSize) / static_cast<float>(detail::kSincRadius);
        const std::int64_t lower_signed = static_cast<std::int64_t>(lower);

        float output = 0.0f;
        for (std::int64_t tap = -static_cast<std::int64_t>(detail::kSincRadius) + 1;
             tap <= static_cast<std::int64_t>(detail::kSincRadius); ++tap) {
            const float table_pos = std::abs(fraction - static_cast<float>(tap)) * kTableScale;
            const std::size_t index =
                std::min<std::size_t>(detail::kSincTableSize - 1, static_cast<std::size_t>(table_pos));
            const float table_fraction = table_pos - static_cast<float>(index);
            const float kernel =
                sinc_table_[index] + table_fraction * (sinc_table_[index + 1] - sinc_table_[index]);
            const std::uint64_t ring_index = static_cast<std::uint64_t>(lower_signed + tap) % kSynthesisRingSize;
            output += synthesis_ring_[ring_index] * kernel;
        }

        resample_position_ += static_cast<double>(current_ratio_);
        const auto floor_position = static_cast<std::uint64_t>(resample_position_);
        // 保留 kSincRadius 个尾随抽头不清理，供下一次插值读取
        const auto clear_until = floor_position > detail::kSincRadius ? floor_position - detail::kSincRadius : 0;
        while (synthesis_cleared_until_ < clear_until) {
            synthesis_ring_[static_cast<std::size_t>(synthesis_cleared_until_ % kSynthesisRingSize)] = 0.0f;
            ++synthesis_cleared_until_;
        }
        return output;
    }

    void heapPush(HeapItem item) noexcept {
        if (heap_size_ == kHeapCapacity) {
            return;
        }
        std::size_t index = heap_size_++;
        heap_[index] = item;
        while (index > 0) {
            const std::size_t parent = (index - 1) / 2;
            if (heap_[parent].magnitude >= heap_[index].magnitude) {
                break;
            }
            std::swap(heap_[parent], heap_[index]);
            index = parent;
        }
    }

    HeapItem heapPop() noexcept {
        const HeapItem result = heap_[0];
        heap_[0] = heap_[--heap_size_];
        std::size_t index = 0;
        while (true) {
            const std::size_t left = index * 2 + 1;
            if (left >= heap_size_) {
                break;
            }
            const std::size_t right = left + 1;
            const std::size_t largest =
                right < heap_size_ && heap_[right].magnitude > heap_[left].magnitude ? right : left;
            if (heap_[index].magnitude >= heap_[largest].magnitude) {
                break;
            }
            std::swap(heap_[index], heap_[largest]);
            index = largest;
        }
        return result;
    }

    float sample_rate_ = 48000.0f;
    float target_ratio_ = 1.0f;
    float current_ratio_ = 1.0f;
    TransientMode transient_mode_ = TransientMode::Flux;
    bool first_frame_ = true;

    qwqdsp_spectral::RealFFT fft_;
    std::array<float, kWindowSize> window_{};
    std::array<float, kWindowSize> input_ring_{};
    std::array<float, kSynthesisRingSize> synthesis_ring_{};
    std::array<float, kFrameSize> fft_input_{};
    std::array<float, kFrameSize + 2> fft_output_{};
    std::array<float, detail::kSincTableSize + 1> sinc_table_{};

    std::array<float, kNumBins> previous_analysis_phase_{};
    std::array<float, kNumBins> previous_synthesis_phase_{};
    std::array<float, kNumBins> previous_magnitude_{};
    std::array<float, kNumBins> current_re_{};
    std::array<float, kNumBins> current_im_{};
    std::array<float, kNumBins> current_analysis_phase_{};
    std::array<float, kNumBins> current_synthesis_phase_{};
    std::array<float, kNumBins> current_magnitude_{};

    std::array<std::uint8_t, kNumBins> processed_bins_{};
    std::array<HeapItem, kHeapCapacity> heap_{};
    std::size_t heap_size_ = 0;

    std::array<float, kFluxHistorySize> flux_history_{};
    std::array<float, kFluxHistorySize> flux_scratch_{};
    std::size_t flux_history_count_ = 0;
    std::size_t flux_history_write_ = 0;
    std::size_t flux_cooldown_ = 0;
    float transient_mean_ = 0.0f;
    float transient_variance_ = 0.0f;
    std::size_t transient_frames_ = 0;
    std::size_t transient_cooldown_ = 0;
    bool transient_has_mean_ = false;

    std::array<std::size_t, kMaxSuperFluxBands> superflux_band_starts_{};
    std::array<std::size_t, kMaxSuperFluxBands> superflux_band_sizes_{};
    std::array<std::size_t, kMaxSuperFluxBands> superflux_band_offsets_{};
    std::array<float, kMaxSuperFluxWeights> superflux_weights_{};
    std::array<float, kMaxSuperFluxBands> superflux_current_{};
    std::array<float, kMaxSuperFluxBands> superflux_previous_max_{};
    std::size_t superflux_band_count_ = 0;
    std::size_t superflux_weight_count_ = 0;
    std::size_t superflux_cooldown_ = 0;

    // DSPark 谱通量检测器状态（复用 superflux 滤波器组，独立通量状态）
    std::array<float, kMaxSuperFluxBands> dsflux_current_{};
    std::array<float, kMaxSuperFluxBands> dsflux_max_prev_{};
    std::array<float, kDsFluxWindow> dsflux_history_{};
    std::array<float, kDsFluxWindow> dsflux_scratch_{};
    std::size_t dsflux_pos_ = 0;
    std::size_t dsflux_filled_ = 0;
    std::size_t dsflux_cooldown_ = 0;

    // 显示支持：音频线程写入，GUI 线程读取
    std::uint64_t input_position_ = 0;
    std::atomic<std::uint64_t> display_position_{0};
    std::array<float, kDisplaySize> display_ring_{};
    std::size_t display_write_ = 0;
    std::array<std::atomic<std::uint64_t>, kMaxTransientMarks> transient_marks_{};
    std::atomic<std::uint32_t> transient_mark_write_{0};

    std::size_t input_write_ = 0;
    std::size_t collected_samples_ = 0;
    std::size_t samples_since_frame_ = 0;
    float samples_to_frame_ = 0.0f;
    bool analysis_started_ = false;
    std::uint64_t synthesis_frame_start_ = 0;
    std::uint64_t synthesis_available_ = 0;
    std::uint64_t synthesis_cleared_until_ = 0;
    double resample_position_ = 0.0;
    bool resampler_started_ = false;
};

// ------------------------------------------------------------
// PitchShiftDemo
// ------------------------------------------------------------

/**
 * @brief 实时 PGHI 移调演示：瞬态选择、bypass、音高旋钮、波形与瞬态显示
 */
class PitchShiftDemo {
public:
    /**
     * @brief 初始化实时处理器和音高旋钮
     */
    void init() {
        spectral_.init(kSampleRate);
        pitch_knob_.set_bound((kWindowWidth - 92) / 2, 158, 92, 98)
            .set_title("Pitch Shift")
            .set_range(-12.0f, 12.0f, 1.0f, 0.0f)
            .set_value(0.0f)
            .set_sensitivity(2)
            .set_name_font_size(15)
            .set_number_font_size(13)
            .set_fore_color(Color{235, 235, 238, 255})
            .set_bg_color(Color{24, 26, 29, 255});
        pitch_knob_.value_to_text_function = [](float value) { return std::format("{:+.0f} st", value); };
        pitch_knob_.on_value_change = [this](float value) { pitch_shift_.store(value, std::memory_order_relaxed); };
    }

    /**
     * @brief 处理音频设备提供的一段单声道采样
     * @param[in] input 单声道输入缓冲
     * @param[out] output 单声道输出缓冲
     * @param[in] frame_count 缓冲中的采样数
     */
    void processAudio(const float* input, float* output, std::size_t frame_count) {
        const int mode = transient_mode_.load(std::memory_order_relaxed);
        const float semitones = pitch_shift_.load(std::memory_order_relaxed);
        const bool bypass = bypass_.load(std::memory_order_relaxed);
        spectral_.setTransientMode(static_cast<TransientMode>(mode));
        spectral_.setPitchShift(semitones);
        // 无论是否 bypass 都保持 DSP 运行，切换时状态保持热
        spectral_.process(input, output, frame_count);
        if (bypass) {
            // 完全绕过移调器：输出直通输入
            std::copy_n(input, frame_count, output);
        }
    }

    /**
     * @brief 绘制瞬态选择器、bypass 开关、音高旋钮、波形和音频状态
     * @param[in] audio_running 音频设备是否成功启动
     */
    void draw(bool audio_running) {
        static constexpr Color kBackground{18, 20, 23, 255};
        static constexpr Color kPanel{24, 26, 29, 255};
        static constexpr Color kText{232, 233, 235, 255};
        static constexpr Color kMuted{126, 132, 139, 255};
        static constexpr Color kAccent{76, 201, 151, 255};
        static constexpr Color kWarning{225, 172, 74, 255};

        ClearBackground(kBackground);
        DrawRectangle(0, 0, kWindowWidth, 64, kPanel);
        DrawRectangle(0, 63, kWindowWidth, 1, Color{48, 52, 57, 255});
        DrawText("PLAYING", 22, 15, 23, kText);
        DrawText("REAL-TIME / MONO", 23, 42, 10, kMuted);

        const Color status_color = audio_running ? kAccent : kWarning;
        DrawCircle(kWindowWidth - 126, 30, 5.0f, status_color);
        DrawText(audio_running ? "AUDIO ONLINE" : "AUDIO OFFLINE", kWindowWidth - 114, 24, 12, status_color);

        drawTransientSelector();
        drawBypassButton();
        pitch_knob_.display();
        drawWaveform();
        DrawText("right click: reset", kWindowWidth - 118, kWindowHeight - 17, 10, kMuted);
    }

    /**
     * @brief 响应 miniaudio 的全双工单声道数据请求
     * @param[in] device 发起回调的音频设备
     * @param[out] output_buffer 单声道浮点输出缓冲
     * @param[in] input_buffer 单声道浮点输入缓冲
     * @param[in] frame_count 本次请求的音频帧数
     */
    static void audioCallback(ma_device* device, void* output_buffer, const void* input_buffer, ma_uint32 frame_count) {
        auto* output = static_cast<float*>(output_buffer);
        const auto* input = static_cast<const float*>(input_buffer);
        auto* demo = static_cast<PitchShiftDemo*>(device->pUserData);
        if (output == nullptr) {
            return;
        }
        if (input == nullptr || demo == nullptr) {
            std::fill_n(output, frame_count, 0.0f);
            return;
        }
        demo->processAudio(input, output, static_cast<std::size_t>(frame_count));
    }
private:
    /**
     * @brief 绘制瞬态检测模式选择器（None / Flux / SuperFlux / Vocoder / DSPark）
     */
    void drawTransientSelector() {
        static constexpr Color kMuted{126, 132, 139, 255};
        constexpr float kLeft = 20.0f;
        constexpr float kTop = 82.0f;
        constexpr float kGap = 5.0f;
        constexpr float kHeight = 31.0f;
        constexpr float kWidth = (static_cast<float>(kWindowWidth) - 2.0f * kLeft - 4.0f * kGap) / 5.0f;

        DrawText("TRANSIENT DETECTION", static_cast<int>(kLeft), 67, 10, kMuted);

        const Vector2 mouse = GetMousePosition();
        const int selected = transient_mode_.load(std::memory_order_relaxed);

        for (std::size_t i = 0; i < kTransientModeNames.size(); ++i) {
            const Rectangle bounds{kLeft + static_cast<float>(i) * (kWidth + kGap), kTop, kWidth, kHeight};
            const bool hovered = CheckCollisionPointRec(mouse, bounds);
            const bool active = selected == static_cast<int>(i);
            const Color border = hovered ? Color{190, 194, 198, 255} : Color{75, 80, 85, 255};
            const Color text = active ? Color{18, 20, 23, 255} : Color{178, 182, 187, 255};
            if (active) {
                DrawRectangleRec(bounds, Color{232, 233, 235, 255});
            }
            else {
                DrawRectangleLinesEx(bounds, 1.0f, border);
            }
            const int font_size = 12;
            const int text_width = MeasureText(kTransientModeNames[i], font_size);
            DrawText(kTransientModeNames[i], static_cast<int>(bounds.x + (bounds.width - text_width) * 0.5f),
                     static_cast<int>(bounds.y + 9.0f), font_size, text);
            if (hovered && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
                transient_mode_.store(static_cast<int>(i), std::memory_order_relaxed);
            }
        }
    }

    /**
     * @brief 绘制 bypass 开关按钮
     */
    void drawBypassButton() {
        static constexpr Color kText{232, 233, 235, 255};
        static constexpr Color kAccent{76, 201, 151, 255};
        constexpr float kLeft = 20.0f;
        constexpr float kTop = 121.0f;
        constexpr float kWidth = 130.0f;
        constexpr float kHeight = 26.0f;

        const bool bypass = bypass_.load(std::memory_order_relaxed);
        const Rectangle bounds{kLeft, kTop, kWidth, kHeight};
        const bool hovered = CheckCollisionPointRec(GetMousePosition(), bounds);
        if (bypass) {
            DrawRectangleRec(bounds, kAccent);
        }
        else {
            const Color border = hovered ? Color{190, 194, 198, 255} : Color{75, 80, 85, 255};
            DrawRectangleLinesEx(bounds, 1.0f, border);
        }
        const char* label = bypass ? "BYPASS: ON" : "BYPASS: OFF";
        const int font_size = 13;
        const int text_width = MeasureText(label, font_size);
        const Color text_color = bypass ? Color{18, 20, 23, 255} : kText;
        DrawText(label, static_cast<int>(bounds.x + (bounds.width - text_width) * 0.5f),
                 static_cast<int>(bounds.y + 5.0f), font_size, text_color);
        if (hovered && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            bypass_.store(!bypass, std::memory_order_relaxed);
        }
    }

    /**
     * @brief 绘制最近 1 秒的输入波形与瞬态检测点竖线
     */
    void drawWaveform() {
        static constexpr Color kPanel{24, 26, 29, 255};
        static constexpr Color kBorder{48, 52, 57, 255};
        static constexpr Color kText{232, 233, 235, 255};
        static constexpr Color kMuted{126, 132, 139, 255};
        static constexpr Color kWave{76, 201, 151, 255};
        static constexpr Color kMark{225, 172, 74, 255};
        static constexpr Color kCenterLine{60, 64, 69, 255};

        constexpr float kLeft = 20.0f;
        constexpr float kTop = 272.0f;
        constexpr float kWidth = static_cast<float>(kWindowWidth) - 2.0f * kLeft;
        constexpr float kHeight = 180.0f;
        const Rectangle panel{kLeft, kTop, kWidth, kHeight};
        DrawRectangleRec(panel, kPanel);
        DrawRectangleLinesEx(panel, 1.0f, kBorder);
        const std::string title =
            std::format("INPUT WAVEFORM ({:g} s)", RealtimePitchShifter::kDisplayTimeSeconds);
        DrawText(title.c_str(), static_cast<int>(kLeft + 8.0f), static_cast<int>(kTop + 6.0f), 12, kText);

        const int mode = std::clamp(transient_mode_.load(std::memory_order_relaxed), 0, 4);
        DrawText(kTransientModeNames[static_cast<std::size_t>(mode)], static_cast<int>(kLeft + kWidth - 90.0f),
                 static_cast<int>(kTop + 7.0f), 10, mode == 0 ? kMuted : kMark);

        constexpr float kWaveTop = kTop + 26.0f;
        constexpr float kWaveBottom = kTop + kHeight - 8.0f;
        const float mid_y = (kWaveTop + kWaveBottom) * 0.5f;
        const float half_height = (kWaveBottom - kWaveTop) * 0.5f;
        DrawLine(static_cast<int>(kLeft + 1.0f), static_cast<int>(mid_y), static_cast<int>(kLeft + kWidth - 2.0f),
                 static_cast<int>(mid_y), kCenterLine);

        const std::uint64_t position = spectral_.displayPosition();
        const int num_columns = static_cast<int>(kWidth) - 2;
        constexpr float kTotalSamples = static_cast<float>(RealtimePitchShifter::kDisplaySize);
        const float samples_per_pixel = kTotalSamples / static_cast<float>(num_columns);
        const std::int64_t position_signed = static_cast<std::int64_t>(position);

        // 每列绘制该列时间区间内样本的 min/max 包络（最右列为最新样本）
        for (int x = 0; x < num_columns; ++x) {
            const std::int64_t end =
                position_signed
                - static_cast<std::int64_t>(std::floor(static_cast<float>(num_columns - 1 - x) * samples_per_pixel));
            const std::int64_t start = std::max<std::int64_t>(
                0, position_signed - static_cast<std::int64_t>(std::ceil(static_cast<float>(num_columns - x) * samples_per_pixel)));
            float minimum = 1.0f;
            float maximum = -1.0f;
            for (std::int64_t sample = start; sample < end; ++sample) {
                const float value = spectral_.displaySampleAt(static_cast<std::uint64_t>(sample));
                minimum = std::min(minimum, value);
                maximum = std::max(maximum, value);
            }
            minimum = std::clamp(minimum, -1.0f, 1.0f);
            maximum = std::clamp(maximum, -1.0f, 1.0f);
            const float x_pos = kLeft + 1.0f + static_cast<float>(x);
            DrawLine(static_cast<int>(x_pos), static_cast<int>(mid_y - maximum * half_height), static_cast<int>(x_pos),
                     static_cast<int>(mid_y - minimum * half_height), kWave);
        }

        // 瞬态检测点竖线
        if (mode != 0) {
            const std::uint32_t count = spectral_.transientMarkCount();
            const std::uint32_t start = count > static_cast<std::uint32_t>(RealtimePitchShifter::kMaxTransientMarks)
                                            ? count - static_cast<std::uint32_t>(RealtimePitchShifter::kMaxTransientMarks)
                                            : 0;
            for (std::uint32_t i = start; i < count; ++i) {
                const std::uint64_t mark = spectral_.transientMark(i);
                if (mark > position) {
                    continue;
                }
                const std::uint64_t age = position - mark;
                if (age >= static_cast<std::uint64_t>(RealtimePitchShifter::kDisplaySize)) {
                    continue;
                }
                const float x = kLeft + 1.0f + static_cast<float>(num_columns - 1)
                              - static_cast<float>(age) / samples_per_pixel;
                DrawLine(static_cast<int>(x), static_cast<int>(kWaveTop), static_cast<int>(x),
                         static_cast<int>(kWaveBottom), kMark);
            }
        }
    }

    RealtimePitchShifter spectral_;
    Knob pitch_knob_;
    std::atomic<float> pitch_shift_{0.0f};
    std::atomic<int> transient_mode_{1}; // Flux
    std::atomic<bool> bypass_{false};
};

} // namespace playing

/**
 * @brief 启动实时 PGHI 移调演示程序
 * @return 正常退出返回零，音频设备不可用时仍保留界面
 */
int main() {
    using namespace playing;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, kWindowTitle);
    SetTargetFPS(60);

    static PitchShiftDemo demo;
    demo.init();

    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.playback.format = ma_format_f32;
    config.playback.channels = 1;
    config.sampleRate = static_cast<ma_uint32>(kSampleRate);
    config.dataCallback = PitchShiftDemo::audioCallback;
    config.pUserData = &demo;
    config.periodSizeInMilliseconds = 5;

    ma_device device{};
    const ma_result init_result = ma_device_init(nullptr, &config, &device);
    ma_result audio_result = init_result;
    bool audio_running = false;
    if (init_result == MA_SUCCESS) {
        audio_result = ma_device_start(&device);
        audio_running = audio_result == MA_SUCCESS;
    }
    if (!audio_running) {
        std::cout << std::format("miniaudio: device start failed ({})\n", static_cast<int>(audio_result));
    }

    while (!WindowShouldClose()) {
        BeginDrawing();
        demo.draw(audio_running);
        EndDrawing();
    }

    if (init_result == MA_SUCCESS) {
        if (audio_running) {
            ma_device_stop(&device);
        }
        ma_device_uninit(&device);
    }
    CloseWindow();
    return 0;
}
