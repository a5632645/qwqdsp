#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <numbers>
#include "miniaudio.h"
#include "raylib.h"
#include "slider.hpp"

// ------------------------------------------------------------
// 常量
// ------------------------------------------------------------
static constexpr int kWindowWidth = 480;
static constexpr int kWindowHeight = 240;
static constexpr float kSampleRate = 48000.0f;

static constexpr float kMinF0 = 75.0f;
static constexpr float kMaxF0 = 500.0f;

static constexpr int kInputBufSize = 8192;   // 须为 2 的幂
static constexpr int kOutputBufSize = 32768;  // 须为 2 的幂
static constexpr int kAlgorithmLatency = 2048;// 合成超前读取的样本数
static constexpr int kGrainMaxLen = 2048;    // 最大 Grain 长度

// MPM 音高检测参数
static constexpr int kMpmMinLag = static_cast<int>(kSampleRate / kMaxF0);
static constexpr int kMpmMaxLag = static_cast<int>(kSampleRate / kMinF0);
static constexpr int kMpmFrameSize = 2048;

// 窗函数 LUT
static constexpr int kWinTableSize = 2048;

// ------------------------------------------------------------
// 工具
// ------------------------------------------------------------
static inline int WrapIndex(int idx, int size) noexcept {
    return idx & (size - 1);
}

// ------------------------------------------------------------
// FormantGrain — 单周期声道共振峰粒度
// ------------------------------------------------------------
struct FormantGrain {
    std::array<float, kGrainMaxLen> samples{};
    int length = 0;
    int t0 = 0;
    bool is_unvoiced = false;
};



// ------------------------------------------------------------
// Hanning 窗 LUT
// ------------------------------------------------------------
static std::array<float, kWinTableSize> s_hann_lut;

static void InitHannLut() noexcept {
    for (int i = 0; i < kWinTableSize; ++i) {
        s_hann_lut[i] = 0.5f * (1.0f - std::cos(2.0f * std::numbers::pi_v<float> * i / (kWinTableSize - 1)));
    }
}

static inline float HannWin(int idx, int len) noexcept {
    const float pos = static_cast<float>(idx) / static_cast<float>(len) * (kWinTableSize - 1);
    const int i0 = static_cast<int>(pos);
    const int i1 = std::min(i0 + 1, kWinTableSize - 1);
    const float f = pos - static_cast<float>(i0);
    return s_hann_lut[i0] + f * (s_hann_lut[i1] - s_hann_lut[i0]);
}

// ------------------------------------------------------------
// 实时 PSOLA — 分析器/合成器解耦架构
// ------------------------------------------------------------
/**
 * @brief 工业级实时 TD-PSOLA，分析器与合成器通过 Triple Buffer 解耦。
 *
 * ┌──────────┐              ┌────────────┐
 * │ Analyzer │── latest ──&gt;│ Synthesizer│
 * │ (GCI)    │   (hold)    │  (synTime) │
 * └──────────┘              └────────────┘
 *
 * 分析器每检测到 GCI 就覆写 latest_grain_，合成器在每个 synTime 位置直读。
 * 无队列/无缓冲，不存在 full/empty 问题。
 *
 * - 分析器：逐样本 GCI + NSDF 音高检测，提取 FormantGrain
 * - 合成器：synTime 按目标周期步进驱动，每步渲染一个 Grain
 * - 浊音 Grain 保持原始时宽（Formant 不变）
 * - 清音 Grain 加相位随机化抖动
 * - formant_shift != 1 时对 Grain 做插值重采样
 */
class RealtimePsola {
public:
    void Init(float sample_rate) noexcept {
        sr_ = sample_rate;
        min_period_ = static_cast<int>(sr_ / kMaxF0);
        max_period_ = static_cast<int>(sr_ / kMinF0);
        def_period_ = (min_period_ + max_period_) / 2;

        std::fill(input_buf_.begin(), input_buf_.end(), 0.0f);
        std::fill(output_buf_.begin(), output_buf_.end(), 0.0f);
        latest_grain_ = FormantGrain{};
        synth_grain_ = FormantGrain{};
        phasor_phase_ = 0.0f;
        syn_rpos_ = 0;
        input_wpos_ = 0;
        prev_epoch_pos_ = -1;
        prev_period_ = static_cast<float>(def_period_);
        pitch_scale_smoothed_ = 1.0f;
        gci_state_ = 0;
        pitch_scale_smoothed_ = 1.0f;
        gci_threshold_ = 0.01f;
        gci_peak_ = 0.0f;
        gci_peak_pos_ = 0;
        gci_since_last_ = 0;
        lpc_xm1_ = 0.0f;
        mpm_prev_period_ = 0.0f;
        mpm_counter_ = 0;
    }

    // ── 原子参数 ──
    std::atomic<float> pitch_shift_semitones_{0.0f};
    std::atomic<float> formant_shift_{1.0f};
    std::atomic<float> pitch_scale_target_{1.0f};

    // --------------------------------------------------------
    // 音频线程主入口
    // --------------------------------------------------------
    void Process(const float* input, float* output, int num_samples) noexcept {
        // ── 参数平滑 ──
        {
            const float target = pitch_scale_target_.load(std::memory_order_relaxed);
            constexpr float kSmooth = 0.005f;
            pitch_scale_smoothed_ += kSmooth * (target - pitch_scale_smoothed_);
        }
        const float pitch_scale = pitch_scale_smoothed_;
        const float formant_scale = formant_shift_.load(std::memory_order_relaxed);

        for (int i = 0; i < num_samples; ++i) {
            const float x = input[i];

            // ── 第 1 步：写入输入缓冲 + 分析器 ──
            input_buf_[input_wpos_] = x;
            RunAnalyzer(input_wpos_, x);
            input_wpos_ = WrapIndex(input_wpos_ + 1, kInputBufSize);

            // ── 第 2 步：Phasor 合成 ──
            // 合成写位置 = 读位置 + 固定超前量（保证 OLA 左半不越界）
            const int syn_wpos = WrapIndex(syn_rpos_ + kAlgorithmLatency, kOutputBufSize);

            // Phasor 步进：频率 = pitch_scale / t0
            const float t0 = (latest_grain_.length > 0)
                ? static_cast<float>(latest_grain_.t0)
                : static_cast<float>(def_period_);
            phasor_phase_ += pitch_scale / t0;

            if (phasor_phase_ >= 1.0f) {
                phasor_phase_ -= 1.0f;

                // 快照 latest_grain_ 到合成器私有缓冲区
                synth_grain_ = latest_grain_;

                if (synth_grain_.length > 0) {
                    // 按 formant_shift 重采样
                    if (std::abs(formant_scale - 1.0f) > 1e-4f && formant_scale > 0.01f) {
                        synth_grain_ = ResampleGrain(synth_grain_, formant_scale);
                    }

                    // OLA 叠加到输出缓冲
                    const int half = synth_grain_.length / 2;
                    for (int j = 0; j < synth_grain_.length; ++j) {
                        const int o = WrapIndex(syn_wpos + j - half + kOutputBufSize, kOutputBufSize);
                        output_buf_[o] += synth_grain_.samples[j];
                    }
                }
            }

            // ── 第 3 步：读取输出 ──
            output[i] = output_buf_[syn_rpos_];
            output_buf_[syn_rpos_] = 0.0f;  // 清零已消费
            syn_rpos_ = WrapIndex(syn_rpos_ + 1, kOutputBufSize);
        }
    }

private:
    // --------------------------------------------------------
    // 分析器
    // --------------------------------------------------------
    void RunAnalyzer(int pos, float x) noexcept {
        if (mpm_counter_++ >= static_cast<int>(0.01f * sr_)) {
            mpm_counter_ = 0;
            RunMpm();
        }
        RunGci(pos, x);
    }

    // --------------------------------------------------------
    // MPM 音高检测 (McLeod Pitch Method)
    // --------------------------------------------------------
    /**
     * @brief MPM 核心：计算 NSDF → 找所有正零交叉峰 →
     *        选第一个超过 0.93×最高峰值的峰（抛物线插值）
     */
    void RunMpm() noexcept {
        const int fs = std::min(kMpmFrameSize, kInputBufSize);
        float frame[kMpmFrameSize];
        for (int i = 0; i < fs; ++i) {
            frame[i] = input_buf_[WrapIndex(input_wpos_ - fs + i + kInputBufSize, kInputBufSize)];
        }

        // ── 计算 NSDF ──
        float nsdf[kMpmMaxLag + 1]{};
        for (int lag = kMpmMinLag; lag <= kMpmMaxLag && lag < fs; ++lag) {
            float sum_diff = 0.0f, sum_sq = 0.0f;
            for (int i = 0; i < fs - lag; ++i) {
                const float d = frame[i] - frame[i + lag];
                sum_diff += d * d;
                sum_sq += frame[i] * frame[i] + frame[i + lag] * frame[i + lag];
            }
            if (sum_sq > 1e-10f)
                nsdf[lag] = 1.0f - sum_diff / (0.5f * sum_sq);
        }

        // ── MPM Peak Picking：找所有正零交叉包围的峰 ──
        // 跳过开头正区
        int pos = kMpmMinLag;
        while (pos <= kMpmMaxLag && nsdf[pos] > 0) ++pos;
        // 跳过负区
        while (pos <= kMpmMaxLag && nsdf[pos] <= 0) ++pos;
        if (pos > kMpmMaxLag) return;

        struct Peak { int idx; float val; };
        Peak peaks[64];
        int n_peaks = 0;
        int cur_max = -1;

        while (pos <= kMpmMaxLag) {
            if (nsdf[pos] > nsdf[pos - 1] && nsdf[pos] >= nsdf[pos + 1]) {
                if (cur_max < 0 || nsdf[pos] > nsdf[cur_max])
                    cur_max = pos;
            }
            ++pos;
            if (pos <= kMpmMaxLag && nsdf[pos] <= 0) {
                if (cur_max >= 0 && n_peaks < 64)
                    peaks[n_peaks++] = {cur_max, nsdf[cur_max]};
                cur_max = -1;
                while (pos <= kMpmMaxLag && nsdf[pos] <= 0) ++pos;
            }
        }
        if (cur_max >= 0 && n_peaks < 64)
            peaks[n_peaks++] = {cur_max, nsdf[cur_max]};

        if (n_peaks == 0) return;

        // ── 找最高峰 ──
        float highest = 0.0f;
        for (int i = 0; i < n_peaks; ++i)
            if (peaks[i].val > highest) highest = peaks[i].val;

        if (highest < 0.15f) return; // 信号太弱

        // ── 选第一个超过 0.93×highest 的峰（抛物线插值）──
        constexpr float kMpmCutoff = 0.93f;
        const float cutoff = kMpmCutoff * highest;
        float best_period = 0.0f;

        for (int i = 0; i < n_peaks; ++i) {
            if (peaks[i].val >= cutoff) {
                const int idx = peaks[i].idx;
                const float ym1 = nsdf[idx - 1];
                const float y0  = nsdf[idx];
                const float yp1 = nsdf[idx + 1];
                const float denom = ym1 - 2.0f * y0 + yp1;
                best_period = static_cast<float>(idx);
                if (std::abs(denom) > 1e-10f)
                    best_period += 0.5f * (ym1 - yp1) / denom;
                break;
            }
        }

        if (best_period > 0.0f) {
            mpm_period_ = best_period;
            mpm_voiced_ = (highest > 0.2f);
            mpm_prev_period_ = best_period;
        }
    }

    // --------------------------------------------------------
    // GCI 逐样本状态机
    // --------------------------------------------------------
    void RunGci(int pos, float x) noexcept {
        const float residual = 0.5f * x - 0.5f * lpc_xm1_;
        lpc_xm1_ = x;
        const float abs_res = std::abs(residual);

        gci_threshold_ = 0.92f * gci_threshold_ + 0.08f * abs_res;
        const float peak_thresh = gci_threshold_ * 2.5f;
        ++gci_since_last_;

        switch (gci_state_) {
        case 0:
            if (abs_res > peak_thresh && gci_since_last_ > min_period_ / 2) {
                gci_state_ = 1;
                gci_peak_ = abs_res;
                gci_peak_pos_ = pos;
            }
            break;
        case 1:
            if (abs_res > gci_peak_) { gci_peak_ = abs_res; gci_peak_pos_ = pos; }
            if (abs_res < gci_peak_ * 0.6f || gci_since_last_ > max_period_) {
                int period = (prev_epoch_pos_ >= 0)
                    ? gci_peak_pos_ - prev_epoch_pos_
                    : static_cast<int>(mpm_period_ > 0 ? mpm_period_ : static_cast<float>(def_period_));
                // 用 MPM 周期修正 GCI 周期
                if (mpm_period_ > 0 && prev_epoch_pos_ >= 0) {
                    const float ratio = static_cast<float>(period) / mpm_period_;
                    if (ratio < 0.5f) period = static_cast<int>(mpm_period_);       // 半周期 → 修正
                    else if (ratio > 1.8f) period = static_cast<int>(mpm_period_);   // 倍周期 → 修正
                }
                if (period >= min_period_ * 0.7f && period <= max_period_ * 1.3f) {
                    PublishGrain(gci_peak_pos_, period);
                    prev_epoch_pos_ = gci_peak_pos_;
                    prev_period_ = static_cast<float>(period);
                } else {
                    // 周期异常时仍然发布一个默认周期 Grain
                    const int c = std::clamp(period, min_period_, max_period_);
                    PublishGrain(gci_peak_pos_, c);
                    prev_epoch_pos_ = gci_peak_pos_;
                }
                gci_state_ = 2;
            }
            break;
        case 2:
            if (gci_since_last_ > prev_period_ * 0.6f) gci_state_ = 0;
            break;
        }
    }

    // --------------------------------------------------------
    // FormantGrain 提取与发布
    // --------------------------------------------------------
    void PublishGrain(int epoch_pos, int t0) noexcept {
        // 直接覆写 latest_grain_，无队列/无缓冲
        const int half = t0;
        const int len = std::min(2 * t0, kGrainMaxLen);

        for (int i = 0; i < len; ++i) {
            const int src = WrapIndex(epoch_pos + i - half + kInputBufSize, kInputBufSize);
            latest_grain_.samples[i] = input_buf_[src] * HannWin(i, len);
        }

        latest_grain_.length = len;
        latest_grain_.t0 = t0;
        latest_grain_.is_unvoiced = !mpm_voiced_;
    }

    // --------------------------------------------------------
    // 合成器 — 颗粒重采样
    // --------------------------------------------------------
    /** 按 ratio 对 Grain 做线性插值重采样 */
    static FormantGrain ResampleGrain(const FormantGrain& src, float ratio) noexcept {
        FormantGrain dst;
        const int rlen = static_cast<int>(static_cast<float>(src.length) / ratio);
        if (rlen < 2 || rlen > kGrainMaxLen) return src;
        dst.length = rlen;
        dst.t0 = static_cast<int>(static_cast<float>(src.t0) / ratio);
        dst.is_unvoiced = src.is_unvoiced;
        const float inv_ratio = 1.0f / ratio;
        for (int i = 0; i < rlen; ++i) {
            const float sp = static_cast<float>(i) * ratio;
            const int i0 = static_cast<int>(sp);
            const int i1 = std::min(i0 + 1, src.length - 1);
            const float f = sp - static_cast<float>(i0);
            dst.samples[i] = src.samples[i0] + f * (src.samples[i1] - src.samples[i0]);
        }
        return dst;
    }

    // --------------------------------------------------------
    // 成员
    // --------------------------------------------------------
    float sr_ = kSampleRate;
    int min_period_{}, max_period_{}, def_period_{};

    std::array<float, kInputBufSize> input_buf_{};
    int input_wpos_ = 0;

    FormantGrain latest_grain_{};   // 分析器覆写
    FormantGrain synth_grain_{};    // 合成器私有拷贝
    float phasor_phase_ = 0.0f;
    int syn_rpos_ = 0;              // 输出读指针
    float pitch_scale_smoothed_ = 1.0f;

    std::array<float, kOutputBufSize> output_buf_{};

    int gci_state_ = 0;
    float gci_threshold_ = 0.0f, gci_peak_ = 0.0f;
    int gci_peak_pos_ = 0, gci_since_last_ = 0;
    int prev_epoch_pos_ = -1;
    float prev_period_ = 0.0f;
    float lpc_xm1_ = 0.0f;

    int mpm_counter_ = 0;
    float mpm_period_ = 0.0f, mpm_prev_period_ = 0.0f;
    bool mpm_voiced_ = false;

};


// 全局实例
static RealtimePsola s_psola;

// ------------------------------------------------------------
//  miniaudio 回调
// ------------------------------------------------------------
extern "C" void MaCallback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount) {
    (void)pDevice;
    auto* input = static_cast<const float*>(pInput);
    auto* output = static_cast<float*>(pOutput);
    if (!input || !output) return;
    s_psola.Process(input, output, static_cast<int>(frameCount));
}

// ------------------------------------------------------------
//  旋钮工具
// ------------------------------------------------------------
static Knob MakeKnob(int& kx, int ky, int w, int h, const char* title,
                     float vmin, float vmax, float vstep, float vdef) {
    Knob knob;
    knob.set_bound(kx, ky, w, h)
        .set_title(title)
        .set_range(vmin, vmax, vstep, vdef)
        .set_value(vdef)
        .set_name_font_size(9)
        .set_number_font_size(9)
        .set_fore_color(Color{200, 200, 210, 255})
        .set_bg_color(Color{30, 30, 35, 255});
    kx += w + 14;
    return knob;
}

// ------------------------------------------------------------
//  main
// ------------------------------------------------------------
int main(void) {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWindowWidth, kWindowHeight, "Realtime PSOLA — Async Analysis/Synth");
    SetTargetFPS(60);

    InitHannLut();
    s_psola.Init(kSampleRate);
    s_psola.pitch_scale_target_.store(1.0f, std::memory_order_relaxed);
    s_psola.pitch_shift_semitones_.store(0.0f, std::memory_order_relaxed);
    s_psola.formant_shift_.store(1.0f, std::memory_order_relaxed);

    ma_device_config config = ma_device_config_init(ma_device_type_duplex);
    config.capture.format = ma_format_f32;
    config.capture.channels = 1;
    config.playback.format = ma_format_f32;
    config.playback.channels = 1;
    config.sampleRate = (ma_uint32)kSampleRate;
    config.dataCallback = MaCallback;
    config.pUserData = nullptr;
    config.periodSizeInMilliseconds = 10;

    ma_device device;
    bool audio_ok = (ma_device_init(nullptr, &config, &device) == MA_SUCCESS);
    if (audio_ok) ma_device_start(&device);

    constexpr int kKnobW = 80, kKnobH = 82;
    constexpr int kGap = 30;
    int kx = (kWindowWidth - (2 * kKnobW + kGap)) / 2;
    int ky = (kWindowHeight - kKnobH) / 2 - 10;

    Knob knob_pitch = MakeKnob(kx, ky, kKnobW, kKnobH, "Pitch", -12.0f, 12.0f, 0.5f, 0.0f);
    Knob knob_formant = MakeKnob(kx, ky, kKnobW, kKnobH, "Formant", 0.5f, 2.0f, 0.05f, 1.0f);

    knob_pitch.on_value_change = [](float v) {
        s_psola.pitch_scale_target_.store(std::pow(2.0f, v / 12.0f), std::memory_order_relaxed);
        s_psola.pitch_shift_semitones_.store(v, std::memory_order_relaxed);
    };
    knob_pitch.value_to_text_function = [](float v) -> std::string {
        return (v >= 0) ? TextFormat("+%.1f st", v) : TextFormat("%.1f st", v);
    };
    knob_formant.on_value_change = [](float v) {
        s_psola.formant_shift_.store(v, std::memory_order_relaxed);
    };
    knob_formant.value_to_text_function = [](float v) -> std::string {
        return TextFormat("%.2fx", v);
    };

    {
        const float pv = knob_pitch.get_value();
        s_psola.pitch_scale_target_.store(std::pow(2.0f, pv / 12.0f), std::memory_order_relaxed);
        s_psola.pitch_shift_semitones_.store(pv, std::memory_order_relaxed);
        s_psola.formant_shift_.store(knob_formant.get_value(), std::memory_order_relaxed);
    }

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(Color{10, 10, 12, 255});
        knob_pitch.display();
        knob_formant.display();
        DrawText("Realtime PSOLA — Phasor-driven Synthesis", 10, kWindowHeight - 24, 10, Color{80, 80, 90, 255});
        EndDrawing();
    }

    if (audio_ok) { ma_device_stop(&device); ma_device_uninit(&device); }
    CloseWindow();
    return 0;
}
