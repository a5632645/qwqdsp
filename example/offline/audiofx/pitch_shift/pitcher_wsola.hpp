#pragma once

#include <qwqdsp/spectral/real_fft.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>

namespace qwqdsp_test {

// ------------------------------------------------------------
// PitcherWsola — WSOLA 音高移位器
// ------------------------------------------------------------
/**
 * @brief 基于环形缓冲的 WSOLA 音高移位器（移植自 qdelay 的 Pitcher）。
 *
 * 原版是 Saike/Joep Vanlier 的 JSFX（WSOLA Pitch shifting library, MIT）
 * 的 C++/JUCE 移植，此处改写为离线单声道版本，FFT 后端换成项目内的
 * qwqdsp_spectral::RealFFT。
 *
 * 原理：
 *   输入逐样本写入环形缓冲，读头以 (1 + speed) 倍速率前进产生音高移动；
 *   读头越过交叉淡化边界时跳回缓冲中部，并用 ACF（自相关）寻找两段波形
 *   的最佳对齐位置，通过余弦交叉淡化平滑拼接，避免跳变咔哒声。
 *
 * 与目录内其他算法的差异：
 *   - 时域算法，输出长度与输入相同（时间-音高耦合，不做时间拉伸）
 *   - 逐样本流式处理，连续调用 Process 会延续内部状态
 */
class PitcherWsola {
public:
    /**
     * @brief 窗口模式（对应原版三种预设）。
     */
    enum class WindowMode {
        kSmall,   ///< 鼓组：窗 128，关闭 ACF
        kMedium,  ///< 通用：窗 1024，启用 ACF
        kLarge,   ///< 人声：窗 2048，启用 ACF
    };

    PitcherWsola() {
        InitInternal(WindowMode::kMedium);
    }

    /**
     * @brief 设置窗口模式并重新初始化内部状态（含 FFT 表，可能分配内存）。
     */
    void SetWindowMode(WindowMode mode) {
        InitInternal(mode);
    }

    /**
     * @brief 以半音设置音高偏移（原版语义：speed = 2^(semitones/12) - 1）。
     *
     * 正半音升调，负半音降调；平滑在逐样本处理中完成。
     */
    void SetSemitones(float semitones) noexcept {
        target_speed_ = std::pow(2.0f, semitones / 12.0f) - 1.0f;
    }

    /**
     * @brief 以比例设置音高偏移（kp>1 升调，kp=1 不变），换算为半音。
     */
    void SetPitchShift(float kp) noexcept {
        SetSemitones(12.0f * std::log2(kp));
    }

    /**
     * @brief 离线处理整段输入，返回与输入等长的输出。
     *
     * 内部按采样流式处理，多次调用会延续状态；独立处理新音频前请先 Reset()。
     *
     * @param input 输入单声道信号
     * @return 音高移位后的输出信号
     */
    std::vector<float> Process(std::span<const float> input) {
        std::vector<float> output;
        output.reserve(input.size());
        for (const float x : input) {
            Update(x);
            output.push_back(out_sample_);
        }
        return output;
    }

    /**
     * @brief 清空内部状态（保留窗口模式与音高设置，不重新分配）。
     */
    void Reset() noexcept {
        ring_.Clear();
        read_head_speed_ = 0.0f;
        last_head_speed_ = read_head_speed_;
        read_head1_ = static_cast<float>(ring_.Size()) * 0.5f;
        read_head2_ = -1.0f;
        fade_ = false;
        fader_ = CosineFade{};
    }

private:
    // ------------------------------------------------------------
    // RingBuffer — 单声道环形缓冲（4 点 Catmull-Rom 插值读取）
    // ------------------------------------------------------------
    class RingBuffer {
    public:
        void Init(int size) {
            size_ = size;
            buffer_.assign(static_cast<size_t>(size), 0.0f);
            while (write_pos_ >= size)
                write_pos_ -= size;
        }

        void Clear() noexcept {
            std::fill(buffer_.begin(), buffer_.end(), 0.0f);
            write_pos_ = 0;
        }

        void Write(float sample) {
            buffer_[static_cast<size_t>(write_pos_)] = sample;
            write_pos_ = (write_pos_ + 1) % size_;
        }

        int Wrap(int idx) const noexcept {
            while (idx < 0)
                idx += size_;
            while (idx >= size_)
                idx -= size_;
            return idx;
        }

        float Read(float dtime) const noexcept {
            float read_pos = static_cast<float>(write_pos_) - dtime;
            while (read_pos < 0.0f)
                read_pos += static_cast<float>(size_);
            while (read_pos >= static_cast<float>(size_))
                read_pos -= static_cast<float>(size_);

            const int i1 = static_cast<int>(std::floor(read_pos));
            const float t = read_pos - static_cast<float>(i1);
            const int i0 = Wrap(i1 - 1);
            const int i2 = Wrap(i1 + 1);
            const int i3 = Wrap(i1 + 2);

            const float t2 = t * t;
            const float t3 = t2 * t;
            const float a0 = -0.5f * t3 + t2 - 0.5f * t;
            const float a1 = 1.5f * t3 - 2.5f * t2 + 1.0f;
            const float a2 = -1.5f * t3 + 2.0f * t2 + 0.5f * t;
            const float a3 = 0.5f * t3 - 0.5f * t2;

            return a0 * buffer_[static_cast<size_t>(i0)]
                 + a1 * buffer_[static_cast<size_t>(i1)]
                 + a2 * buffer_[static_cast<size_t>(i2)]
                 + a3 * buffer_[static_cast<size_t>(i3)];
        }

        void CopyFromBuffer(float* target, int delay, int copy_length) const noexcept {
            int start = write_pos_ - delay;
            while (start < 0)
                start += size_;
            start %= size_;

            int remaining = copy_length;
            int pos = start;
            int out = 0;
            while (remaining > 0) {
                const int chunk = std::min(size_ - pos, remaining);
                std::copy_n(buffer_.begin() + pos, chunk, target + out);
                out += chunk;
                remaining -= chunk;
                pos = (pos + chunk) % size_;
            }
        }

        int Size() const noexcept {
            return size_;
        }

    private:
        std::vector<float> buffer_;
        int write_pos_ = 0;
        int size_ = 0;
    };

    // ------------------------------------------------------------
    // CosineFade — 递归余弦交叉淡化曲线
    // ------------------------------------------------------------
    struct CosineFade {
        float y0 = 0.0f;
        float y1 = 0.0f;
        float y2 = 0.0f;
        float b1 = 0.0f;
        float count = 1.0f;
        float prev_len = 0.0f;  // 上一次淡化长度（Nc）
        float weight = 0.0f;    // 当前淡化权重 [0,1]

        void Prepare(float n) noexcept {
            if (n < 2.0f)
                n = 2.0f;
            prev_len = n;
            const float step = std::numbers::pi_v<float> / (n - 1.0f);
            const float half_pi = std::numbers::pi_v<float> * 0.5f;
            count = n;
            b1 = 2.0f * std::cos(step);
            y1 = std::sin(half_pi - step);
            y2 = std::sin(half_pi - 2.0f * step);
        }

        void Resize(float n) noexcept {
            if (count <= 0.0f || prev_len <= 0.0f)
                return;
            const float clamped = std::max(2.0f, n);
            const float step = std::numbers::pi_v<float> / (clamped - 1.0f);
            b1 = 2.0f * std::cos(step);
            count *= clamped / prev_len;  // 保持淡化进度比例
            prev_len = clamped;
        }

        void Eval() noexcept {
            y0 = b1 * y1 - y2;
            y2 = y1;
            y1 = y0;
            weight = std::clamp(0.5f * (y0 + 1.0f), 0.0f, 1.0f);
            count -= 1.0f;
        }
    };

    // ------------------------------------------------------------
    // 窗口 / FFT 配置
    // ------------------------------------------------------------
    // FFT 长度 = 窗口 * 2
    static constexpr int kOrderSmall = 7;
    static constexpr int kOrderMedium = 10;
    static constexpr int kOrderLarge = 11;

    static constexpr int kWindowSmall = 1 << kOrderSmall;    // 128
    static constexpr int kWindowMedium = 1 << kOrderMedium;  // 1024
    static constexpr int kWindowLarge = 1 << kOrderLarge;    // 2048

    static int WindowSize(WindowMode mode) noexcept {
        switch (mode) {
            case WindowMode::kSmall:
                return kWindowSmall;
            case WindowMode::kMedium:
                return kWindowMedium;
            case WindowMode::kLarge:
                return kWindowLarge;
            default:
                return kWindowMedium;
        }
    }

    static bool AcfFlag(WindowMode mode) noexcept {
        switch (mode) {
            case WindowMode::kSmall:
                return false;
            case WindowMode::kMedium:
            case WindowMode::kLarge:
                return true;
            default:
                return true;
        }
    }

    // ------------------------------------------------------------
    // 内部初始化
    // ------------------------------------------------------------
    void InitInternal(WindowMode mode) {
        mode_ = mode;
        const int win_size = WindowSize(mode);
        ring_.Init(win_size * 2);
        acf_ = AcfFlag(mode);
        // 默认不改音高。原版 init 置 1.0，依赖插件每块调用 setSpeed 覆盖，
        // 离线版本默认应为 kp=1（速度偏移 0）。
        read_head_speed_ = 0.0f;
        last_head_speed_ = read_head_speed_;
        cross_fade_samples_ = win_size / 2;
        fft_.Init(static_cast<size_t>(win_size * 2));
        fftmem1_.assign(static_cast<size_t>(win_size * 2), 0.0f);
        fftmem2_.assign(static_cast<size_t>(win_size * 2), 0.0f);
        spect1_.assign(static_cast<size_t>(win_size * 2) + 2, 0.0f);
        spect2_.assign(static_cast<size_t>(win_size * 2) + 2, 0.0f);
        read_head1_ = static_cast<float>(ring_.Size()) * 0.5f;
        read_head2_ = -1.0f;
        fade_ = false;
        fader_ = CosineFade{};
    }

    // ------------------------------------------------------------
    // ACF 相位对齐
    // ------------------------------------------------------------
    /**
     * @brief 计算两段半窗数据的 ACF 并返回其峰值位置（亚采样精度）。
     *
     * 对 fftmem1_（源区段）与 fftmem2_（目标区段）做频域互相关：
     *   ACF = IFFT(FFT(B) · conj(FFT(A)))
     * 峰值位置即两段波形最佳对齐的偏移量。
     */
    float ComputeMaxAcfPosition() {
        const int fft_size = WindowSize(mode_) * 2;
        const int win_size = WindowSize(mode_);

        fft_.FFT(fftmem2_.data(), spect1_.data());
        fft_.FFT(fftmem1_.data(), spect2_.data());

        // 频域互相关：spect1[k] *= conj(spect2[k])
        // CCS 布局中 DC 与 Nyquist 为实 bin，其余为复数 bin，逐 bin 相乘。
        const size_t num_bins = static_cast<size_t>(fft_size) / 2 + 1;
        for (size_t k = 0; k < num_bins; ++k) {
            const float re1 = spect1_[2 * k];
            const float im1 = spect1_[2 * k + 1];
            const float re2 = spect2_[2 * k];
            const float im2 = spect2_[2 * k + 1];
            // a * conj(b) = (re1*re2 + im1*im2) + i*(im1*re2 - re1*im2)
            spect1_[2 * k] = re1 * re2 + im1 * im2;
            spect1_[2 * k + 1] = im1 * re2 - re1 * im2;
        }

        // 逆变换得到互相关序列（1/L 缩放不影响峰值位置）
        fft_.IFFT(spect1_.data(), fftmem1_.data());

        // 在前四分之一窗口内找峰
        int max_idx = 0;
        float cmax = -1e7f;
        const int search = win_size / 4 - 1;
        for (int idx = 0; idx < search; ++idx) {
            const float current = fftmem1_[static_cast<size_t>(idx)];
            if (current > cmax) {
                cmax = current * 1.05f;
                max_idx = idx;
            }
        }

        // 抛物线亚采样修正
        const float yc = fftmem1_[static_cast<size_t>(max_idx)];
        const float yl = fftmem1_[static_cast<size_t>(std::max(max_idx - 1, 0))];
        const float yr = fftmem1_[static_cast<size_t>(std::min(max_idx + 1, search - 1))];
        const float denom = -2.0f * yc + yl + yr;
        const float correction = (denom == 0.0f) ? 0.0f : (yl - yr) / denom;

        float peak_idx = static_cast<float>(max_idx);
        if (std::fabs(correction) < 4.0f)
            peak_idx += correction;

        return std::max(peak_idx, 0.0f);
    }

    // ------------------------------------------------------------
    // 逐样本处理
    // ------------------------------------------------------------
    void Update(float sample) {
        // 速度平滑（原版 setSpeed 的 EMA）与淡化长度自适应
        read_head_speed_ = read_head_speed_ * 0.999f + 0.001f * target_speed_;
        if (std::fabs(read_head_speed_ - last_head_speed_) > 0.0001f)
            fader_.Resize(static_cast<float>(cross_fade_samples_)
                              / std::max(1.1f, std::fabs(read_head_speed_)) - 16.0f);
        last_head_speed_ = read_head_speed_;

        read_head1_ -= read_head_speed_;

        if (fader_.count <= 0.0f) {
            fade_ = false;
            fader_.count = 1.0f;
            read_head1_ = read_head2_;
        }

        ring_.Write(sample);
        const float l1 = ring_.Read(read_head1_);

        if (fade_) {
            // 交叉淡化中：两路读头混合
            const float l2 = ring_.Read(read_head2_);
            read_head2_ -= read_head_speed_;
            fader_.Eval();
            out_sample_ = l1 * fader_.weight + l2 * (1.0f - fader_.weight);
        }
        else {
            int src;
            int target;
            bool crit;

            if (read_head_speed_ <= 0.0f) {
                // 降调：读头向后越过边界时跳回前部
                src = ring_.Size() - cross_fade_samples_;
                crit = read_head1_ > static_cast<float>(src);
                target = cross_fade_samples_ + cross_fade_samples_;
            }
            else {
                // 升调：读头向前越过边界时跳回中部
                src = cross_fade_samples_;
                crit = read_head1_ <= static_cast<float>(src);
                target = ring_.Size() - cross_fade_samples_;
            }

            if (crit) {
                // 越过拼接边界，用 ACF 对齐后初始化交叉淡化
                float cmax_position = 0.0f;
                if (acf_) {
                    std::fill(fftmem1_.begin(), fftmem1_.end(), 0.0f);
                    std::fill(fftmem2_.begin(), fftmem2_.end(), 0.0f);
                    ring_.CopyFromBuffer(fftmem1_.data(), src, cross_fade_samples_);
                    ring_.CopyFromBuffer(fftmem2_.data(), target, cross_fade_samples_);
                    cmax_position = ComputeMaxAcfPosition();
                }
                fader_.Prepare(std::floor(static_cast<float>(cross_fade_samples_)
                                          / std::max(1.1f, std::fabs(read_head_speed_))));
                read_head2_ = static_cast<float>(target) - cmax_position;
                fade_ = true;
            }

            out_sample_ = l1;
        }
    }

    // ---- 状态 ----
    WindowMode mode_ = WindowMode::kMedium;
    RingBuffer ring_;
    CosineFade fader_;

    qwqdsp_spectral::RealFFT fft_;
    std::vector<float> fftmem1_;  // ACF 输入 A / ACF 输出
    std::vector<float> fftmem2_;  // ACF 输入 B
    std::vector<float> spect1_;   // CCS 频谱缓冲（N+2）
    std::vector<float> spect2_;   // CCS 频谱缓冲（N+2）

    float read_head1_ = 0.0f;
    float read_head2_ = 0.0f;
    float last_head_speed_ = 0.0f;
    float read_head_speed_ = 0.0f;
    float target_speed_ = 0.0f;  // 目标读头速度偏移（由半音换算）
    float out_sample_ = 0.0f;

    bool fade_ = false;
    bool acf_ = false;
    int cross_fade_samples_ = 0;
};

} // namespace qwqdsp_test
