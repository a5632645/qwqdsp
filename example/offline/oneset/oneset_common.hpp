#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <qwqdsp/spectral/real_fft.hpp>
#include <span>
#include <vector>

namespace qwqdsp_test {
namespace detail {

/** 将相位折叠至 (-pi, pi] */
static inline float PrincArg(float x) noexcept {
    const float two_pi = 2.0f * std::numbers::pi_v<float>;
    return x - two_pi * std::floor(x / two_pi + 0.5f);
}

/** periodic Hann 窗：w[i] = 0.5 * (1 - cos(2πi/N)) */
static inline void MakePeriodicHann(std::span<float> w) noexcept {
    const size_t n = w.size();
    for (size_t i = 0; i < n; ++i) {
        w[i] = 0.5f * (1.0f - std::cos(2.0f * std::numbers::pi_v<float> * static_cast<float>(i) / static_cast<float>(n)));
    }
}

// ------------------------------------------------------------
// StftFrame
// ------------------------------------------------------------
/**
 * @brief 轻量 STFT 帧分析器：Hann 加窗 + RealFFT，输出幅度/相位与历史。
 *
 * 逐帧调用，维护上一帧（prev_mag_）与上上帧（prev_phase2_）历史，
 * 供谱通量 / 复合域 / SuperFlux 等 ODF 使用。状态随帧前进更新。
 */
class StftFrame {
public:
    /** 构造并配置帧长 */
    void Init(size_t fft_size) {
        fft_size_ = fft_size;
        bins_ = fft_size_ / 2 + 1;
        window_.resize(fft_size);
        MakePeriodicHann(window_);
        fft_.Init(fft_size);
        fft_in_.resize(fft_size);
        fft_out_.resize(fft_size + 2);
        mag_.resize(bins_);
        phase_.resize(bins_);
        prev_mag_.assign(bins_, 0.0f);
        prev_phase_.assign(bins_, 0.0f);
        prev_phase2_.assign(bins_, 0.0f);
    }

    /** 分析一帧：@p frame 指向加窗前的帧起点（长度 fft_size_） */
    void Analyze(const float* frame) {
        // 历史滚动须在计算新帧之前：分析前的 mag_/phase_ 是上一帧的值
        prev_phase2_ = prev_phase_;
        prev_phase_ = phase_;
        prev_mag_ = mag_;

        for (size_t j = 0; j < fft_size_; ++j)
            fft_in_[j] = frame[j] * window_[j];

        fft_.FFT(fft_in_.data(), fft_out_.data());

        for (size_t j = 0; j < bins_; ++j) {
            const float re = fft_out_[2 * j];
            const float im = fft_out_[2 * j + 1];
            mag_[j] = std::sqrt(re * re + im * im);
            phase_[j] = std::atan2(im, re);
        }
    }

    size_t fft_size_ = 0;
    size_t bins_ = 0;
    std::vector<float> window_;
    qwqdsp_spectral::RealFFT fft_;
    std::vector<float> fft_in_;
    std::vector<float> fft_out_;
    std::vector<float> mag_;
    std::vector<float> phase_;
    std::vector<float> prev_mag_;
    std::vector<float> prev_phase_;
    std::vector<float> prev_phase2_;
};

// ------------------------------------------------------------
// PeakPick
// ------------------------------------------------------------
/**
 * @brief 自适应峰值拾取（librosa.util.peak_pick / Boeck 2012）。
 *
 * 对 ODF 序列找局部最大值（pre_max 帧窗口内最大），超过前后移动均值
 * + delta，且与上一个 onset 间隔大于 wait 帧时判为 onset。
 * 时间窗默认：pre_max=30ms、post_max=30ms、pre_avg=100ms、
 * post_avg=70ms、wait=30ms（与 librosa/DSPark 一致）。
 *
 * @param odf       onset 强度序列（每帧一个值）
 * @param hop       帧 hop（采样数），用于把 ms 窗口换算成帧
 * @param sample_rate 采样率
 * @param delta     高于移动均值的阈值（越大越少触发）
 * @return onset 帧下标（升序）
 */
static inline std::vector<size_t> PeakPick(std::span<const float> odf, size_t hop, float sample_rate, float delta) {
    const auto ms_to_frames = [&](double ms) {
        return std::max<size_t>(1, static_cast<size_t>(std::lround(ms * 0.001 * static_cast<double>(sample_rate) / static_cast<double>(hop))));
    };
    const size_t pre_max = ms_to_frames(30.0);
    const size_t post_max = ms_to_frames(30.0);
    const size_t pre_avg = ms_to_frames(100.0);
    const size_t post_avg = ms_to_frames(70.0);
    const size_t wait = ms_to_frames(30.0);

    const size_t nf = odf.size();
    std::vector<size_t> onsets;
    if (nf == 0)
        return onsets;

    const int64_t kBig = int64_t(1) << 60;
    int64_t last_frame = -kBig;

    for (size_t f = 0; f < nf; ++f) {
        // 局部最大：前后 post_max/pre_max 帧窗口内最大
        bool is_max = true;
        for (size_t j = f > pre_max ? f - pre_max : 0; j <= std::min(nf - 1, f + post_max); ++j) {
            if (odf[j] > odf[f]) {
                is_max = false;
                break;
            }
        }
        if (!is_max)
            continue;

        // 前后移动均值（含当前帧，librosa 行为）
        float sum = 0.0f;
        size_t cnt = 0;
        const size_t avg_lo = f > pre_avg ? f - pre_avg : 0;
        const size_t avg_hi = std::min(nf - 1, f + post_avg);
        for (size_t j = avg_lo; j <= avg_hi; ++j) {
            sum += odf[j];
            ++cnt;
        }
        const float mean = (cnt > 0) ? sum / static_cast<float>(cnt) : 0.0f;
        if (odf[f] < mean + delta)
            continue;
        if (static_cast<int64_t>(f) - last_frame <= static_cast<int64_t>(wait))
            continue;

        onsets.push_back(f);
        last_frame = static_cast<int64_t>(f);
    }

    return onsets;
}

} // namespace detail

// ------------------------------------------------------------
// OnsetResult
// ------------------------------------------------------------
/**
 * @brief 检测结果：ODF 序列 + 峰值拾取得到的 onset 位置。
 */
struct OnsetResult {
    std::vector<float> odf;         ///< 每帧一个 onset 强度值
    std::vector<size_t> onset_frames; ///< onset 帧下标（升序）
    std::vector<size_t> onset_samples; ///< onset 采样位置（帧下标 × hop）
};

} // namespace qwqdsp_test
