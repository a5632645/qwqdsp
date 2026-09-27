#pragma once
#include <algorithm>
#include <cmath>
#include <numbers>
#include <queue>
#include <qwqdsp/spectral/real_fft.hpp>
#include <span>
#include <vector>

namespace qwqdsp_test {

// ------------------------------------------------------------
// PhaseGradientVocoder
// ------------------------------------------------------------
/**
 * @brief 相位梯度声码器，使用对称 FFT（La == Ls）。
 *
 * 时间拉伸通过 Aa/As 比率控制，音高移动通过输出重采样实现。
 */
class PhaseGradientVocoder {
public:
    static constexpr size_t kDefaultFrameSize = 4096;
    static constexpr size_t kDefaultOverSample = 2;

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

        std::vector<float> curr_re(bins);
        std::vector<float> curr_im(bins);
        std::vector<float> curr_phase(bins);

        // 输出
        output_.assign(As * N + Ls, 0.0f);

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
                curr_re[j] = fft_out_[2 * j];
                curr_im[j] = fft_out_[2 * j + 1];
            }

            // ---- 堆积分 ----
            if (i == 0) {
                for (size_t j = 0; j < bins; ++j)
                    curr_phase[j] = std::atan2(curr_im[j], curr_re[j]);
            }
            else {
                CalcPGHI(prev_phase_, prev_re_, prev_im_, curr_re, curr_im, curr_phase, ra, rs);
            }

            std::swap(prev_re_, curr_re);
            std::swap(prev_im_, curr_im);
            std::swap(prev_phase_, curr_phase);

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
        output_.clear();
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

    struct HeapItem {
        float magnitude;
        size_t bin;
        bool is_prev;
    };

    // ------------------------------------------------------------
    // CalcPGHI — Phase Gradient Heap Integration
    // ------------------------------------------------------------
    /**
     * @brief PGHI 算法核心实现（简化版，不含 abstol 和梯形积分）。
     *
     * 论文完整算法 (Phase Vocoder Done Right 2022):
     *
     * Input:  Phase time derivative (Δtϕa) and magnitude s of frames n and n−1,
     *         phase frequency derivative (Δfϕa) for frame n,
     *         estimated phase ϕs for frame n−1 and relative tolerance tol.
     * 输入：  相位时间导数 (Δtϕa) 以及帧 n 和 n−1 的幅度 s，
     *         帧 n 的相位频率导数 (Δfϕa)，
     *         帧 n−1 的估计相位 ϕs 以及相对容差 tol。
     *
     * Output: Phase estimate ϕs for frame n.
     * 输出：  帧 n 的相位估计 ϕs。
     *
     *   1  abstol ← tol · max(s(m,n) ∪ s(m,n−1))
     *   2  Create set I = { m : s(m,n) > abstol }
     *   3  Assign random values to ϕs(m,n) for m ∉ I
     *   4  Construct a max heap for (m,n) tuples
     *   5  Insert (m,n−1) for m ∈ I into the heap
     *   6  while I is not ∅ do
     *   7      (mh, nh) ← remove the top of the heap
     *   8      if nh = n−1 then
     *   9          if (mh, n) ∈ I then
     *  10              ϕs(mh,n) ← ϕs(mh,n−1) + as/2 · ((Δtϕa)(mh,n−1) + (Δtϕa)(mh,n))
     *  11              Remove (mh,n) from I
     *  12              Insert (mh,n) into the heap
     *  13          end if
     *  14      end if
     *  15      if nh = n then
     *  16          if (mh+1, n) ∈ I then
     *  17              ϕs(mh+1,n) ← ϕs(mh,n) + bs/2 · ((Δfϕa)(mh,n) + (Δfϕa)(mh+1,n))
     *  18              Remove (mh+1,n) from I
     *  19              Insert (mh+1,n) into the heap
     *  20          end if
     *  21          if (mh−1, n) ∈ I then
     *  22              ϕs(mh−1,n) ← ϕs(mh,n) − bs/2 · ((Δfϕa)(mh,n) + (Δfϕa)(mh−1,n))
     *  23              Remove (mh−1,n) from I
     *  24              Insert (mh−1,n) into the heap
     *  25          end if
     *  26      end if
     *  27  end while
     *
     * 注：该简化版跳过了 abstol / set I（步骤 1-3），
     *     对所有 bin 执行 PGHI，且时间/频率传播使用简化公式
     *     （不使用论文里推荐的梯形积分 as/2 · (prev + curr)，实测比简化版效果差）。
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

        // Step 2: insert (m, n−1) for all m (简化版跳过 set I)
        for (size_t j = 0; j < bins; ++j) {
            const float mag = std::sqrt(prev_re[j] * prev_re[j] + prev_im[j] * prev_im[j]);
            heap.push({mag, j, true});
        }

        // Step 3: while I is not ∅ (简化版用 count < bins 替代 I 非空判断)
        size_t count = 0;
        while (count < bins) {
            // Step 4: remove top of heap → (mh, nh)
            const auto item = heap.top();
            heap.pop();
            const size_t i = item.bin;

            // Step 5: nh = n−1 → 时间方向传播
            if (item.is_prev) {
                // Step 6: check (mh, n) ∈ I (简化版只检查是否已处理)
                if (!empty[i])
                    continue;

                // Δtϕa 估计：当前帧与上一帧的瞬时频率偏移
                const float cr = curr_re[i] * prev_re[i] + curr_im[i] * prev_im[i];
                const float ci = curr_im[i] * prev_re[i] - curr_re[i] * prev_im[i];
                float dp = std::atan2(ci, cr);
                dp -= two_pi_ra * static_cast<float>(i);

                // Step 7: ϕs(mh,n) = ϕs(mh,n−1) + as · Δtϕa (简化版，不含梯形积分)
                out_phase[i] = WrapToPi(prev_phase[i] + ratio * WrapToPi(dp) + two_pi_rs * static_cast<float>(i));

                // Step 8: remove (mh,n) from I
                empty[i] = false;
                ++count;
                // Step 9: insert (mh,n) into heap
                const float mag = std::sqrt(curr_re[i] * curr_re[i] + curr_im[i] * curr_im[i]);
                heap.push({mag, i, false});
            }
            // Step 10: nh = n → 频率方向传播
            else {
                // Step 11: forward neighbor (mh+1, n)
                if (i + 1 < bins && empty[i + 1]) {
                    // Δfϕa 估计：相邻 bin 的群延迟
                    const float cr = curr_re[i + 1] * curr_re[i] + curr_im[i + 1] * curr_im[i];
                    const float ci = curr_im[i + 1] * curr_re[i] - curr_re[i + 1] * curr_im[i];
                    const float dp = std::atan2(ci, cr);

                    // Step 12: ϕs(mh+1,n) = ϕs(mh,n) + bs · Δfϕa (简化版)
                    out_phase[i + 1] = out_phase[i] + ratio * dp;

                    // Step 13: remove (mh+1,n) from I
                    empty[i + 1] = false;
                    ++count;
                    // Step 14: insert (mh+1,n) into heap
                    const float mag = std::sqrt(curr_re[i + 1] * curr_re[i + 1] + curr_im[i + 1] * curr_im[i + 1]);
                    heap.push({mag, i + 1, false});
                }
                // Step 15: backward neighbor (mh−1, n)
                if (i >= 1 && empty[i - 1]) {
                    const float cr = curr_re[i - 1] * curr_re[i] + curr_im[i - 1] * curr_im[i];
                    const float ci = curr_im[i - 1] * curr_re[i] - curr_re[i - 1] * curr_im[i];
                    const float dp = std::atan2(ci, cr);

                    // Step 16: ϕs(mh−1,n) = ϕs(mh,n) − bs · Δfϕa (简化版)
                    out_phase[i - 1] = out_phase[i] + ratio * dp;

                    // Step 17: remove (mh−1,n) from I
                    empty[i - 1] = false;
                    ++count;
                    // Step 18: insert (mh−1,n) into heap
                    const float mag = std::sqrt(curr_re[i - 1] * curr_re[i - 1] + curr_im[i - 1] * curr_im[i - 1]);
                    heap.push({mag, i - 1, false});
                }
            }
        }
    }

    size_t frame_size_ = kDefaultFrameSize;
    size_t over_sample_ = kDefaultOverSample;
    float kt_ = 1.0f;
    float kp_ = 1.0f;

    qwqdsp_spectral::RealFFT fft_;

    std::vector<float> window_;
    std::vector<float> fft_in_;
    std::vector<float> fft_out_;
    std::vector<float> prev_re_;
    std::vector<float> prev_im_;
    std::vector<float> prev_phase_;
    std::vector<float> output_;
};

} // namespace qwqdsp_test
