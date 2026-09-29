#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <vector>

#include "raylib.h"

#include "log_reassign_grid.hpp"

/**
 * @brief 无窗 NC 重分配模式
 *
 * 两个自由度(是否搬频率 / 是否搬时间)在本帧里各由 ``Mode`` 的一个取值控制，
 * 对应 ``labs/nc_reassign`` 的 freq / time / tf 三个变体:
 *   - ``kFreq``:     频率搬到瞬时频率, 时间留在窗口中心(帧时间轴)
 *   - ``kTime``:     频率留在 bin 中心, 时间搬到群延迟
 *   - ``kFreqTime``: 两者都搬(默认; lab 实测线上能量 0.94 的最优变体)
 */
enum class NcReassignMode {
    kFreq,
    kTime,
    kFreqTime,
};

/**
 * @brief 无窗 NC 谱帧 + **重分配**（低频按"至少 k 个周期"取窗）
 *
 * 变体由模板参数 ``Mode`` 选（见 ``NcReassignMode``）。
 *
 * 与 ``WindowlessNcFrame``（同一前端，但**不做重分配**）和 `NcReassignmentFrame`
 * （矩形窗 FFT 的线性 bin 对，也是不重分配）不同：本帧在**对数频率网格**上每行放一个
 * NC bin，逐样本递归滑动 DFT，然后把每个 bin 的能量按 (瞬时频率, 群延迟) 搬进
 * ``LogReassignGrid``。算子与时间参考的推导、Python 参考实现与全部数值验证见
 * ``qwqdsp/labs/nc_reassign``（README + `compare_all.py` 六方法总对比）。
 *
 * 核心公式（论文 arXiv:2410.07982v3）：
 *   (5) 分量频率  f_left/right = f_c ∓ Fs/(2N)
 *   (7) 窗长      N = round( round(2·f_c/W_NC)·Fs/(2·f_c) )
 *   滑动 DFT      X(n) = W·X(n−1) + x[n] − x[n−N]·W^N,   W = e^{−j·2π·f/Fs}
 *                 （相位锚定窗口起点，无需式 8 的相位校正）
 *
 *   NC 幅度       gain   = sqrt(max(0, −(Re_L·Re_R + Im_L·Im_R))) / N
 *   瞬时频率      Y = X_R − X_L，f_inst = −Fs/2π·arg(Y[n]·conj(Y[n−1]))
 *                 Y 恰好等于对同一段样本做**正弦锥**(sin(πm/N))加窗的滑动 DFT
 *                 （e^{−jω_R m} − e^{−jω_L m} = −2j·sin(πm/N)·e^{−jω_c m}）：
 *                 零额外乘法、旁瓣比矩形窗低约 10 dB、滚降快一倍，低频镜像/旁瓣污染
 *                 因此大降（lab 实测 25 Hz 处逐 bin 误差 67.5 → 22.5 cent；原因与
 *                 数据见 labs/nc_reassign README「频率估计器」）。
 *                 本实现的滑动 DFT 对正频率相位随时间递减，故取负号。
 *   群延迟        δ = wrap(arg(X_R·conj(X_L)) + π)
 *                 ⟨m⟩ = (N−1)/2 − N·δ/(2π) = 能量相对**当前样本**的延迟（样本）
 *                 （静止纯音的 ψ 恒为 −π，故以 π 为参考零点；δ = 0 → 落在窗中心）
 *
 *   频率外推      重分配把能量搬到「窗中心 + N·δ/(2π)」处，而上面的 f_inst 是**窗中心**
 *                 处的值；若不同步外推，落点会偏离真值线（lab 实测 0.9459 vs 0.9986）。
 *                 故 kFreqTime 下按 IF 序列斜率一阶外推：
 *                   f_dep = f_inst + rate·(落点相对窗中心的偏移),  rate = Δf_inst/hop
 *                 斜率用**因果后向差分**（实时拿不到未来帧）：lab 实测 0.9977 vs 离线
 *                 中心差分 0.9986，差 0.1%。kFreq 的偏移为 0，自动退化为窗中心频率。
 *
 * 窗长策略（实测见 lab README「窗长 N 的下限」「低频窗长」「超低频用几个周期的窗」）：
 *   - 纯理论窗长（式 7）在低频极长（25 Hz ≈ 8.5 万样本 = 1.8 s）：扫频/瞬态下 `|X|` 由
 *     驻定相位幅度决定、与 N 无关，而 NC 要除以 N → 增益按 1/N 塌缩，且窗中心不再是
 *     有效时间参考 → 谱图变糊（lab 里 `nc-free` 低频段线上能量 0.749）。
 *   - C++ 现有 ``max_window_s = 0.075 s`` 上限：20–30 Hz 只有 1.5–2.25 个周期，略低于
 *     探测/估计质量的最优点。
 *   - **本帧**：在时间上限之外再加一条 **N ≥ minPeriodsFloor·Fs/f_c 的下限**
 *     （默认 4 个周期，20 Hz → 0.2 s），下限可顶开时间上限但不超过式(7)的自然值。
 *     实测（交叉对数 chirp，±50 cent 线上能量）：低频段 0.876 → **0.942**，全带
 *     0.882 → **0.946**；≥53 Hz 的 bin 完全不受影响（下限不生效）。
 *
 * 时间参考：能量落在估计出的**绝对时刻**（窗中心，或重分配后的 窗中心 + 群延迟偏移）上，
 * 所以静止音对所有 bin 都落在同一列；环缓冲按最长窗取 ``ceil((N_max+1.5)/hop)+1`` 个子列，因此整幅显示
 * 固定滞后 N_max 样本——这是长低频窗的固有延迟（0.075 s 上限下 75 ms，4 周期下限下
 * 200 ms）。
 *
 * 幅度：沿用无窗 NC 的 π 归一化（单位幅度纯音 ≈ 0 dB，与 ``WindowlessNcFrame`` 同刻度），
 * 不做 ``LogReassignGrid::SetWindow`` 那套窗相干增益标定；低频行因 NC bin 带宽重叠、
 * 子格内求和会偏亮（与 C++ ``WindowlessNcFrame`` 的低频行为一致）。
 */
template <typename Colormap, NcReassignMode Mode = NcReassignMode::kFreqTime,
          bool EnableFreqInterp = true>
struct WindowlessNcReassignFrame {
    /**
     * @brief 单个 NC bin 的几何与滑动 DFT 状态
     */
    struct Bin {
        float f_center{};                        // 中心频率(Hz, 仅用于显示/映射)
        double f_left{};                         // 左分量参考频率(Hz)
        double f_right{};                        // 右分量参考频率(Hz)
        int N{};                                 // 窗长(样本)
        // 旋转因子与累加器用 double：float 会随时间累积慢速漂移(论文 IV-B 节)
        std::complex<double> Wl{}, Wr{};         // 每样本旋转因子 W
        std::complex<double> WNr_l{}, WNr_r{};   // 跨窗因子 W^N
        std::complex<double> acc_l{}, acc_r{};   // X(n)
        std::complex<double> prev_l{}, prev_r{}; // X(n−1)，用于瞬时频率
        double if_prev{};                        // 上一帧的窗中心瞬时频率(Hz)，用于求斜率
        bool if_prev_valid{};                    // if_prev 是否已有效
    };

    /**
     * @brief 初始化 bin 布局、环缓冲与 log 重分配网格
     *
     * @param sampleRate        采样率(Hz)
     * @param fftSize           帧大小(由上层 SpectrogramColumn 提供)
     * @param hopSize           帧间步进(样本数)
     * @param zeroPad           零填充因子（本方法无零填充，保留占位，与其它帧接口一致）
     * @param outputHeight      输出行数（= 画布高度，每行一个 NC bin）
     * @param freqMin/freqMax   显示频率范围(Hz)
     * @param dbFloor           幅度下限(dB)
     * @param bandwidthScale    NC bin 带宽缩放(默认 1.0)
     * @param minPeriodsFloor   窗长下限(周期数，默认 4.0)：N ≥ k·Fs/f_c，可顶开时间上限
     */
    void Init(int sampleRate, int fftSize, int hopSize, int zeroPad, int outputHeight,
              float freqMin, float freqMax, float dbFloor, float bandwidthScale = 1.0f,
              float minPeriodsFloor = kDefaultMinPeriods) noexcept {
        sampleRate_ = sampleRate;
        fftSize_ = fftSize;
        hopSize_ = hopSize;
        zeroPad_ = zeroPad; // 占位，无零填充
        outputHeight_ = outputHeight;
        freqMin_ = freqMin;
        freqMax_ = freqMax;
        dbFloor_ = dbFloor;
        logMin_ = std::log10(freqMin);
        logMax_ = std::log10(freqMax);

        // 幅度归一化：NC 前向公式 gain = sqrt(ncSum)/N 对单位幅度纯音稳态增益 ≈ A/π
        // （推导见 paper/lab），乘 π 使纯音峰值 ≈ 0 dB，与 WindowlessNcFrame 同刻度。
        gain_norm_ = static_cast<float>(std::numbers::pi_v<double>);

        const int max_window_samples = std::max(8, static_cast<int>(kMaxWindowS * sampleRate));

        // ── 以频率网格为基准：y 轴每个像素对应一个 NC bin（与 WindowlessNcFrame 同布局）──
        const float log_step = (logMax_ - logMin_) / static_cast<float>(outputHeight);
        std::vector<float> centers(static_cast<size_t>(outputHeight));
        for (int y = 0; y < outputHeight; ++y)
            centers[static_cast<size_t>(y)] = std::pow(10.0f, logMax_ - (y + 0.5f) * log_step);

        bins_.clear();
        bins_.reserve(static_cast<size_t>(outputHeight));
        int max_n = 8;
        for (int y = 0; y < outputHeight; ++y) {
            Bin b{};
            b.f_center = centers[static_cast<size_t>(y)];

            // 带宽 = 相邻两 bin 中心频率之差 × 缩放系数（论文 W_NC = f(i+1) − f(i−1)）
            float f_hi = b.f_center; // 更高频侧
            float f_lo = b.f_center; // 更低频侧
            if (y > 0)
                f_hi = centers[static_cast<size_t>(y - 1)];
            if (y + 1 < outputHeight)
                f_lo = centers[static_cast<size_t>(y + 1)];
            float w_nc = std::max(std::abs(f_hi - f_lo) * bandwidthScale, 1e-3f);

            // (7) 自然窗长；再夹到时间上限，最后取周期数下限（下限可顶开上限，但不超自然值）
            const float q = std::round(2.0f * b.f_center / w_nc);
            const int natural_n = static_cast<int>(std::round(q * sampleRate / (2.0f * b.f_center)));
            int n = std::clamp(natural_n, 8, max_window_samples);
            const int floor_n = static_cast<int>(std::round(minPeriodsFloor * sampleRate / b.f_center));
            n = std::min(natural_n, std::max(n, floor_n));
            b.N = std::max(8, n);

            // (5) 左右分量频率（用 double 保证旋转因子精度）
            b.f_left = static_cast<double>(b.f_center) - static_cast<double>(sampleRate) / (2.0 * b.N);
            b.f_right = static_cast<double>(b.f_center) + static_cast<double>(sampleRate) / (2.0 * b.N);
            b.Wl = std::polar(1.0, -2.0 * std::numbers::pi_v<double> * b.f_left / sampleRate);
            b.Wr = std::polar(1.0, -2.0 * std::numbers::pi_v<double> * b.f_right / sampleRate);
            b.WNr_l = std::polar(1.0, -2.0 * std::numbers::pi_v<double> * b.f_left * b.N / sampleRate);
            b.WNr_r = std::polar(1.0, -2.0 * std::numbers::pi_v<double> * b.f_right * b.N / sampleRate);

            max_n = std::max(max_n, b.N);
            bins_.push_back(b);
        }

        // ── 共享样本环回缓冲（长度 = 最长窗长）──
        ring_.assign(static_cast<size_t>(max_n), 0.0);
        ring_size_ = max_n;
        ring_span_ = max_n;
        n_ = 0;

        // ── log 重分配网格：子列数按最长窗取，保证每个落点都在环内 ──
        // 落点子列坐标 = (N_max + 1 − ⟨m⟩)/hop ∈ [0, (N_max+1.5)/hop]，⟨m⟩ ∈ [−0.5, N−0.5]
        const int ring_cols = static_cast<int>(std::ceil((ring_span_ + 1.5) / hopSize_)) + 1;
        grid_.Init(sampleRate, fftSize, ring_cols, outputHeight, freqMin, freqMax, dbFloor);
        column_.resize(static_cast<size_t>(outputHeight));
    }

    /**
     * @brief 处理一帧：用递归滑动 DFT 演进所有 bin，重分配后输出一列
     *
     * window 与 windowed_frame 传入但被忽略（无窗法）。SpectrogramColumn 提供重叠帧：
     * 第一次调用整帧都是新样本；之后每次只有末尾 hopSize_ 个样本是新增的。
     */
    void Process(std::span<const float> raw_frame, std::span<const float> /*window*/,
                 std::span<const float> /*windowed_frame*/) noexcept {
        constexpr double kTwoPi = 2.0 * std::numbers::pi_v<double>;
        const int new_sample_begin =
            (n_ == 0) ? 0 : std::max(0, static_cast<int>(raw_frame.size()) - hopSize_);

        for (int k = new_sample_begin; k < static_cast<int>(raw_frame.size()); ++k) {
            const double s = static_cast<double>(raw_frame[static_cast<size_t>(k)]);
            for (Bin& b : bins_) {
                const double old = (n_ >= b.N) ? ring_[static_cast<size_t>((n_ - b.N) % ring_size_)] : 0.0;
                b.prev_l = b.acc_l;
                b.prev_r = b.acc_r;
                b.acc_l = b.Wl * b.acc_l + s - old * b.WNr_l;
                b.acc_r = b.Wr * b.acc_r + s - old * b.WNr_r;
            }
            ring_[static_cast<size_t>(n_ % ring_size_)] = s;
            ++n_;
        }

        // ── 逐 bin：幅度 + 瞬时频率 + 群延迟 → 落到网格 ──
        const double inv_two_pi_fs = static_cast<double>(sampleRate_) / kTwoPi;
        for (Bin& b : bins_) {
            if (n_ < b.N) // 窗还没填满（开头一段），不落点
                continue;
            const double nc_sum = -(b.acc_l.real() * b.acc_r.real() + b.acc_l.imag() * b.acc_r.imag());
            if (nc_sum <= 0.0)
                continue;
            const float gain =
                static_cast<float>(std::sqrt(nc_sum + kEpsSample) / b.N) * gain_norm_;

            // ── 时间轴：kTime/kFreqTime 搬到 窗中心 + N·δ/(2π)，kFreq 留在窗中心 ──
            // δ = wrap(arg(X_R·conj(X_L)) + π)，⟨m⟩ = (N−1)/2 − N·δ/(2π)
            double mean_delay = 0.5 * (b.N - 1);
            if constexpr (Mode != NcReassignMode::kFreq) {
                const double delta = WrapPi(std::arg(b.acc_r * std::conj(b.acc_l)) + std::numbers::pi_v<double>);
                mean_delay -= b.N * delta / kTwoPi;
            }

            // ── 频率轴：kFreq/kFreqTime 用 Y = X_R − X_L（正弦锥加窗）求**窗中心**处瞬时频率；
            //    kFreqTime 再按 IF 序列斜率（因果后向差分）外推到落点时刻 ──
            float freq_hz = b.f_center;
            if constexpr (Mode != NcReassignMode::kTime) {
                const std::complex<double> y_now = b.acc_r - b.acc_l;
                const std::complex<double> y_prev = b.prev_r - b.prev_l;
                const double if_center = -inv_two_pi_fs * std::arg(y_now * std::conj(y_prev));
                if constexpr (Mode == NcReassignMode::kFreqTime) {
                    const double offset = 0.5 * (b.N - 1) - mean_delay;   // 落点相对窗中心(样本)
                    const double rate = b.if_prev_valid
                                            ? (if_center - b.if_prev) / static_cast<double>(hopSize_)
                                            : 0.0;                        // Hz/样本
                    freq_hz = static_cast<float>(if_center + rate * offset);
                }
                else {
                    freq_hz = static_cast<float>(if_center);
                }
                b.if_prev = if_center;
                b.if_prev_valid = true;
            }
            if (freq_hz < freqMin_ || freq_hz > freqMax_)
                continue;

            const float col = static_cast<float>((ring_span_ + 1.0 - mean_delay) / hopSize_);

            grid_.AddAtColumn(freq_hz, col, gain);
        }
        grid_.Emit(Colormap::kTable, column_);
    }

    std::span<const Color> GetColumn() const noexcept {
        return {column_.data(), static_cast<size_t>(outputHeight_)};
    }

    int ColumnHeight() const noexcept {
        return outputHeight_;
    }

private:
    static constexpr float kMaxWindowS = 0.075f;         // 低音窗长上限(秒, C++ 现口径)
    static constexpr float kDefaultMinPeriods = 4.0f;    // 窗长下限(周期数, lab 实测最优 ~4)
    static constexpr double kEpsSample = 1e-18;          // NC 增益平方根内的保护

    /// @brief 把相位卷绕到 (−π, π]
    static double WrapPi(double x) noexcept {
        return std::atan2(std::sin(x), std::cos(x));
    }

    int sampleRate_{}, fftSize_{}, hopSize_{}, zeroPad_{}, outputHeight_{};
    int ring_size_{}, ring_span_{};
    int n_{};                                            // 已处理样本计数
    float freqMin_{}, freqMax_{}, logMin_{}, logMax_{}, dbFloor_{};
    float gain_norm_{};                                  // 幅度归一化(乘 π ≈ 单位纯音 0 dB)

    std::vector<Bin> bins_;
    std::vector<double> ring_;                           // 共享样本环回缓冲
    std::vector<Color> column_;
    LogReassignGrid<EnableFreqInterp> grid_;
};
