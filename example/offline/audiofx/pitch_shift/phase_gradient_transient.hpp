#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <numbers>
#include <queue>
#include <qwqdsp/spectral/real_fft.hpp>
#include <span>
#include <vector>

namespace qwqdsp_test {

/** 瞬态检测 ODF 类型：Flux（线性幅度谱通量）或 SuperFlux（对数滤波器组） */
enum class OdfType { Flux, SuperFlux };

// ------------------------------------------------------------
// PhaseGradientTransientVocoder
// ------------------------------------------------------------
/**
 * @brief 相位梯度声码器 + flux/superflux 瞬态检测。
 *
 * 结合 PGHI（Phase Gradient Heap Integration）相位传播与瞬态感知：
 * 每帧计算 flux 或 superflux 的 onset 强度（ODF），瞬态帧（ODF 超阈值）
 * 将合成相位重置为分析相位（保留攻击锐度，避免瞬态被相位传播抹平），
 * 非瞬态帧用 PGHI 堆积分传播相位（保持谐波相干，消除 phasiness）。
 *
 * 通过模板参数 Odf 切换检测算法：
 * - OdfType::Flux：线性幅度谱半波整流差分（librosa 谱通量）
 * - OdfType::SuperFlux：对数频率滤波器组 + 频率最大滤波（Boeck 2013）
 *
 * 时间拉伸通过 Aa/As 比率控制，音高移动通过输出重采样实现
 * （与 PhaseGradientVocoder 相同）。
 */
template <OdfType Odf = OdfType::SuperFlux>
class PhaseGradientTransientVocoder {
public:
    static constexpr size_t kDefaultFrameSize = 4096;
    static constexpr size_t kDefaultOverSample = 2;
    /** 默认瞬态阈值：Flux 用 0.07，SuperFlux 用 0.03（ODF 尺度不同） */
    static constexpr float kDefaultFluxThreshold = 0.07f;
    static constexpr float kDefaultSuperFluxThreshold = 0.03f;
    /** Flux 自适应阈值的近期 ODF 历史窗口（帧数） */
    static constexpr size_t kFluxMedianWindow = 32;
    /** Flux 瞬态判定：ODF > 中位数基线 × 该倍数 */
    static constexpr float kFluxMultiplier = 4.0f;
    /** Flux 瞬态后冷却帧数（避免连续重置打断 PGHI） */
    static constexpr size_t kFluxCooldown = 6;
    /** Flux 瞬态判定下限（中位数基线接近 0 时用绝对下限） */
    static constexpr float kFluxMinOdf = 0.02f;

    void SetFrameSize(size_t n) noexcept {
        frame_size_ = n;
    }

    void SetOverSample(size_t m) noexcept {
        over_sample_ = m;
    }

    void SetTimeStretch(float kt) noexcept {
        kt_ = kt;
    }

    void SetPitchShift(float kp) noexcept {
        kp_ = kp;
    }

    /** 设置采样率（SuperFlux 滤波器组依赖，默认 48000） */
    void SetSampleRate(float sr) noexcept {
        sample_rate_ = sr;
    }

    /** 设置瞬态检测阈值（ODF 超过即视为瞬态帧） */
    void SetTransientThreshold(float t) noexcept {
        transient_threshold_ = t;
    }

    std::vector<float> Process(std::span<const float> input) {
        const float kt = kt_;
        const float kp = kp_;
        const size_t Mf = over_sample_;
        const size_t La = frame_size_;
        const size_t Ls = La; // 对称 FFT

        // 总拉伸比 = kt * kp（之后重采样还原 kp 部分）
        const float total_stretch = kt * kp;
        const size_t Aa = static_cast<size_t>(std::round(1024.0f / std::max(total_stretch, 1e-6f)));
        const size_t As = 1024;

        if (Aa == 0 || La == 0)
            return {};

        const size_t fft_size = La * Mf;
        const size_t bins = fft_size / 2 + 1;
        const size_t N = (input.size() - La) / Aa;

        if (N == 0)
            return {};

        const float ra = static_cast<float>(Aa) / static_cast<float>(fft_size);
        const float rs = static_cast<float>(As) / static_cast<float>(fft_size);

        // 窗函数
        window_.resize(La);
        MakeSinSqWindow(window_);

        // FFT
        fft_.Init(fft_size);

        // 缓冲区
        fft_in_.resize(fft_size);
        fft_out_.resize(fft_size + 2);

        prev_re_.resize(bins);
        prev_im_.resize(bins);
        prev_phase_.resize(bins);
        prev_mag_.resize(bins);

        std::vector<float> curr_re(bins);
        std::vector<float> curr_im(bins);
        std::vector<float> curr_phase(bins);
        std::vector<float> curr_mag(bins);

        // 输出
        output_.assign(As * N + Ls, 0.0f);

        // SuperFlux 滤波器组（依赖 sample_rate 的 binHz 映射，首次调用构建）
        if constexpr (Odf == OdfType::SuperFlux) {
            BuildFilterBank(fft_size, bins);
        }

        size_t transient_frames = 0;

        for (size_t i = 0; i < N; ++i) {
            // ---- 分析 ----
            const float* frame_in = input.data() + Aa * i;

            std::fill(fft_in_.begin(), fft_in_.end(), 0.0f);
            for (size_t j = 0; j < La / 2; ++j)
                fft_in_[j] = frame_in[j + La / 2] * window_[j + La / 2];
            const size_t pad = fft_size - La / 2;
            for (size_t j = 0; j < La / 2; ++j)
                fft_in_[pad + j] = frame_in[j] * window_[j];

            fft_.FFT(fft_in_.data(), fft_out_.data());

            for (size_t j = 0; j < bins; ++j) {
                const float re = fft_out_[2 * j];
                const float im = fft_out_[2 * j + 1];
                curr_re[j] = re;
                curr_im[j] = im;
                curr_mag[j] = std::sqrt(re * re + im * im);
            }

            // ---- 瞬态检测 ----
            const float odf = ComputeOdf(curr_mag, fft_size, bins);
            bool is_transient = false;
            if constexpr (Odf == OdfType::Flux) {
                // Flux：自适应阈值（相对近期 ODF 中位数基线）+ 冷却。
                // 线性幅度谱通量对扫频等持续频谱变化也持续响应，固定阈值
                // 会把大多数帧判为瞬态、打断 PGHI；改为只对显著超过
                // 局部基线的帧重置相位。
                if (i > 0 && flux_cooldown_ == 0) {
                    float baseline = 0.0f;
                    if (!flux_hist_.empty()) {
                        std::vector<float> hist(flux_hist_.begin(), flux_hist_.end());
                        std::nth_element(hist.begin(), hist.begin() + hist.size() / 2, hist.end());
                        baseline = hist[hist.size() / 2];
                    }
                    is_transient = odf > std::max(kFluxMultiplier * baseline, kFluxMinOdf);
                }
                if (is_transient) {
                    flux_cooldown_ = kFluxCooldown;
                    ++transient_frames;
                }
                if (flux_cooldown_ > 0)
                    --flux_cooldown_;
                // 更新历史：瞬态帧不推入，避免基线被异常高峰拉高
                if (!is_transient) {
                    flux_hist_.push_back(odf);
                    if (flux_hist_.size() > kFluxMedianWindow)
                        flux_hist_.pop_front();
                }
            }
            else {
                // SuperFlux：固定阈值（对数滤波器组 ODF 尺度稳定）
                is_transient = (i > 0) && (odf > transient_threshold_);
                if (is_transient)
                    ++transient_frames;
            }

            // ---- 相位决策 ----
            if (i == 0 || is_transient) {
                // 首帧或瞬态帧：相位重置为分析相位（保留攻击锐度）
                for (size_t j = 0; j < bins; ++j)
                    curr_phase[j] = std::atan2(curr_im[j], curr_re[j]);
            }
            else {
                // 非瞬态帧：PGHI 堆积分传播
                CalcPGHI(prev_phase_, prev_re_, prev_im_, curr_re, curr_im, curr_phase, ra, rs);
            }

            std::swap(prev_re_, curr_re);
            std::swap(prev_im_, curr_im);
            std::swap(prev_phase_, curr_phase);
            prev_mag_ = curr_mag;

            // ---- 综合 ----
            for (size_t j = 0; j < bins; ++j) {
                const float mag = std::sqrt(prev_re_[j] * prev_re_[j] + prev_im_[j] * prev_im_[j]);
                const float ph = prev_phase_[j];
                fft_out_[2 * j] = mag * std::cos(ph);
                fft_out_[2 * j + 1] = mag * std::sin(ph);
            }

            fft_.IFFT(fft_out_.data(), fft_in_.data());

            // ifftshift + 加窗 + OLA
            for (size_t j = 0; j < La / 2; ++j) {
                const float s = fft_in_[fft_size - La / 2 + j] * window_[j];
                output_[As * i + j] += s;
            }
            for (size_t j = 0; j < La / 2; ++j) {
                const float s = fft_in_[j] * window_[La / 2 + j];
                output_[As * i + La / 2 + j] += s;
            }
        }

        // 归一化
        const float norm = static_cast<float>(As) / static_cast<float>(La);
        for (auto& x : output_)
            x *= norm;

        // ---- 音高移动：重采样 ----
        if (kp != 1.0f) {
            output_ = LinearResample(output_, kp);
        }

        return output_;
    }

    void Reset() noexcept {
        prev_re_.clear();
        prev_im_.clear();
        prev_phase_.clear();
        prev_mag_.clear();
        output_.clear();
        band_cur_.clear();
        band_max_prev_.clear();
        flux_hist_.clear();
        flux_cooldown_ = 0;
    }

private:
    static float WrapToPi(float x) noexcept {
        return x - 2.0f * std::numbers::pi_v<float> * std::round(x / (2.0f * std::numbers::pi_v<float>));
    }

    static void MakeSinSqWindow(std::span<float> w) noexcept {
        const float kScale = std::sqrt(8.0f / 3.0f);
        const size_t N = w.size();
        for (size_t i = 0; i < N; ++i) {
            const float t = std::numbers::pi_v<float> * static_cast<float>(i) / static_cast<float>(N);
            const float s = std::sin(t);
            w[i] = kScale * s * s;
        }
    }

    static std::vector<float> LinearResample(std::span<const float> in, float ratio) {
        const size_t out_len = static_cast<size_t>(std::round(static_cast<float>(in.size()) / ratio));
        if (out_len < 2)
            return {};
        std::vector<float> out(out_len);
        const float step = static_cast<float>(in.size() - 1) / static_cast<float>(out_len - 1);
        for (size_t i = 0; i < out_len; ++i) {
            const float pos = step * static_cast<float>(i);
            const size_t idx = static_cast<size_t>(pos);
            const float frac = pos - static_cast<float>(idx);
            const size_t nxt = std::min(idx + 1, in.size() - 1);
            out[i] = in[idx] + frac * (in[nxt] - in[idx]);
        }
        return out;
    }

    // ------------------------------------------------------------
    // ODF 计算
    // ------------------------------------------------------------
    /**
     * @brief 计算当前帧的 onset 强度（ODF）。
     *
     * @param mag      当前帧幅度（size = bins）
     * @param fft_size FFT 大小
     * @param bins     bin 数
     * @return ODF 值；越大越可能是瞬态帧
     */
    float ComputeOdf(const std::vector<float>& mag, size_t fft_size, size_t bins) {
        if constexpr (Odf == OdfType::Flux) {
            // 线性幅度半波整流差分（librosa 谱通量）
            float flux = 0.0f;
            for (size_t j = 0; j < bins; ++j) {
                const float d = mag[j] - prev_mag_[j];
                if (d > 0.0f)
                    flux += d;
            }
            return flux / static_cast<float>(bins);
        }
        else {
            // SuperFlux：对数滤波器组 + 频率最大滤波（Boeck 2013）
            const size_t num_bands = band_start_.size();
            if (num_bands == 0 || band_cur_.size() != num_bands) {
                band_cur_.assign(num_bands, 0.0f);
                band_max_prev_.assign(num_bands, 0.0f);
            }

            // band 幅度 + log10 压缩
            const float scale = static_cast<float>(kOdfRefFrame / static_cast<double>(fft_size));
            for (size_t b = 0; b < num_bands; ++b) {
                float acc = 0.0f;
                const size_t start = band_start_[b];
                const size_t cnt = band_count_[b];
                for (size_t w = 0; w < cnt; ++w) {
                    const size_t k = start + w;
                    if (k < bins)
                        acc += mag[k] * band_weights_[band_offset_[b] + w];
                }
                band_cur_[b] = std::log10(acc * scale + 1.0f);
            }

            // 半波整流差分（对前帧最大滤波参考）
            float flux = 0.0f;
            for (size_t b = 0; b < num_bands; ++b) {
                const float d = band_cur_[b] - band_max_prev_[b];
                if (d > 0.0f)
                    flux += d;
            }
            flux /= static_cast<float>(num_bands);

            // 旋转：当前 → 前帧，并对前帧做频率 3 邻域最大滤波
            std::vector<float> band_prev = band_cur_;
            for (size_t b = 0; b < num_bands; ++b) {
                float mx = band_prev[b];
                if (b > 0)
                    mx = std::max(mx, band_prev[b - 1]);
                if (b + 1 < num_bands)
                    mx = std::max(mx, band_prev[b + 1]);
                band_max_prev_[b] = mx;
            }

            return flux;
        }
    }

    /** SuperFlux 对数频率三角滤波器组构建（每倍频程 24 band，27.5Hz–16kHz） */
    void BuildFilterBank(size_t fft_size, size_t bins) {
        if (!band_start_.empty())
            return; // 已构建

        const float bin_hz = sample_rate_ / static_cast<float>(fft_size);
        const float f_max = std::min(kFMaxHz, sample_rate_ * 0.5f * 0.999f);

        std::vector<int> centres;
        for (int i = 0;; ++i) {
            const float f = kFMin * std::pow(2.0f, static_cast<float>(i) / static_cast<float>(kBandsPerOctave));
            if (f > f_max)
                break;
            int bin = static_cast<int>(std::lround(f / bin_hz));
            bin = std::clamp(bin, 0, static_cast<int>(bins) - 1);
            if (centres.empty() || bin > centres.back())
                centres.push_back(bin);
        }

        band_start_.clear();
        band_count_.clear();
        band_offset_.clear();
        band_weights_.clear();

        for (size_t j = 1; j + 1 < centres.size(); ++j) {
            const int lo = centres[j - 1];
            const int ce = centres[j];
            const int hi = centres[j + 1];
            if (!(lo < ce && ce < hi))
                continue;
            band_start_.push_back(static_cast<size_t>(lo));
            band_offset_.push_back(band_weights_.size());
            for (int k = lo; k <= hi; ++k) {
                const float w = (k <= ce) ? static_cast<float>(k - lo) / static_cast<float>(ce - lo)
                                          : static_cast<float>(hi - k) / static_cast<float>(hi - ce);
                band_weights_.push_back(w);
            }
            band_count_.push_back(static_cast<size_t>(hi - lo + 1));
        }
        band_cur_.assign(band_start_.size(), 0.0f);
        band_max_prev_.assign(band_start_.size(), 0.0f);
    }

    // ------------------------------------------------------------
    // CalcPGHI — Phase Gradient Heap Integration
    // ------------------------------------------------------------
    /**
     * @brief PGHI 堆积分相位传播（简化版，与 phase_vocoder2.hpp 一致）。
     *
     * 按幅度从大到小处理频谱 bin：时间方向传播（当前帧与上一帧的瞬时
     * 频率偏移）与频率方向传播（相邻 bin 的群延迟），保持谐波相位相干。
     */
    void CalcPGHI(std::span<const float> prev_phase, std::span<const float> prev_re, std::span<const float> prev_im,
                  std::span<const float> curr_re, std::span<const float> curr_im, std::span<float> out_phase, float ra,
                  float rs) const {
        const size_t bins = prev_phase.size();
        const float ratio = rs / ra;
        const float two_pi_ra = 2.0f * std::numbers::pi_v<float> * ra;
        const float two_pi_rs = 2.0f * std::numbers::pi_v<float> * rs;

        // Step 1: max heap for (m, n) tuples
        std::vector<bool> empty(bins, true);
        auto cmp = [](const HeapItem& a, const HeapItem& b) { return a.magnitude < b.magnitude; };
        std::priority_queue<HeapItem, std::vector<HeapItem>, decltype(cmp)> heap(cmp);

        // Step 2: insert (m, n−1) for all m（简化版跳过 set I）
        for (size_t j = 0; j < bins; ++j) {
            const float mag = std::sqrt(prev_re[j] * prev_re[j] + prev_im[j] * prev_im[j]);
            heap.push({mag, j, true});
        }

        // Step 3: while I is not ∅（简化版用 count < bins 替代）
        size_t count = 0;
        while (count < bins) {
            const auto item = heap.top();
            heap.pop();
            const size_t i = item.bin;

            // 时间方向传播（nh = n−1）
            if (item.is_prev) {
                if (!empty[i])
                    continue;

                // Δtϕa：当前帧与上一帧的瞬时频率偏移
                const float cr = curr_re[i] * prev_re[i] + curr_im[i] * prev_im[i];
                const float ci = curr_im[i] * prev_re[i] - curr_re[i] * prev_im[i];
                float dp = std::atan2(ci, cr);
                dp -= two_pi_ra * static_cast<float>(i);

                out_phase[i] = WrapToPi(prev_phase[i] + ratio * WrapToPi(dp) + two_pi_rs * static_cast<float>(i));

                empty[i] = false;
                ++count;
                const float mag = std::sqrt(curr_re[i] * curr_re[i] + curr_im[i] * curr_im[i]);
                heap.push({mag, i, false});
            }
            // 频率方向传播（nh = n）
            else {
                // 前向邻 bin
                if (i + 1 < bins && empty[i + 1]) {
                    const float cr = curr_re[i + 1] * curr_re[i] + curr_im[i + 1] * curr_im[i];
                    const float ci = curr_im[i + 1] * curr_re[i] - curr_re[i + 1] * curr_im[i];
                    const float dp = std::atan2(ci, cr);
                    out_phase[i + 1] = out_phase[i] + ratio * dp;

                    empty[i + 1] = false;
                    ++count;
                    const float mag = std::sqrt(curr_re[i + 1] * curr_re[i + 1] + curr_im[i + 1] * curr_im[i + 1]);
                    heap.push({mag, i + 1, false});
                }
                // 后向邻 bin
                if (i >= 1 && empty[i - 1]) {
                    const float cr = curr_re[i - 1] * curr_re[i] + curr_im[i - 1] * curr_im[i];
                    const float ci = curr_im[i - 1] * curr_re[i] - curr_re[i - 1] * curr_im[i];
                    const float dp = std::atan2(ci, cr);
                    out_phase[i - 1] = out_phase[i] + ratio * dp;

                    empty[i - 1] = false;
                    ++count;
                    const float mag = std::sqrt(curr_re[i - 1] * curr_re[i - 1] + curr_im[i - 1] * curr_im[i - 1]);
                    heap.push({mag, i - 1, false});
                }
            }
        }
    }

    struct HeapItem {
        float magnitude;
        size_t bin;
        bool is_prev;
    };

    // SuperFlux 滤波器组常量
    static constexpr float kFMin = 27.5f;
    static constexpr float kFMaxHz = 16000.0f;
    static constexpr int kBandsPerOctave = 24;
    static constexpr double kOdfRefFrame = 2048.0;

    size_t frame_size_ = kDefaultFrameSize;
    size_t over_sample_ = kDefaultOverSample;
    float kt_ = 1.0f;
    float kp_ = 1.0f;
    float sample_rate_ = 48000.0f;
    float transient_threshold_ = (Odf == OdfType::Flux) ? kDefaultFluxThreshold : kDefaultSuperFluxThreshold;

    qwqdsp_spectral::RealFFT fft_;

    std::vector<float> window_;
    std::vector<float> fft_in_;
    std::vector<float> fft_out_;
    std::vector<float> prev_re_;
    std::vector<float> prev_im_;
    std::vector<float> prev_phase_;
    std::vector<float> prev_mag_;
    std::vector<float> output_;

    // SuperFlux 滤波器组（CSR 式：start bin + 扁平权重 + per-band offset/count）
    std::vector<size_t> band_start_;
    std::vector<float> band_weights_;
    std::vector<size_t> band_offset_;
    std::vector<size_t> band_count_;
    std::vector<float> band_cur_;
    std::vector<float> band_max_prev_;

    // Flux 自适应瞬态状态（仅 OdfType::Flux 使用）
    std::deque<float> flux_hist_; // 近期非瞬态帧 ODF（中位数基线）
    size_t flux_cooldown_ = 0;
};

} // namespace qwqdsp_test
