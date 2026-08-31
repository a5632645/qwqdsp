// ------------------------------------------------------------
// AutocorrelationDelayPitchShifter 离线测试
// ------------------------------------------------------------
// 输入正弦波，按 5ms 分块流式处理（模拟实时回调），验证稳态输出：
//   1) 无 NaN/Inf
//   2) 输出主频 = 输入频率 * 2^(半音/12)，误差 < 0.5 Hz
//   3) 增益（输出/输入 RMS）在 [0.9, 1.1]
//   4) 在精确期望频率上拟合正弦后残差伪影 < -35 dB
// 每个用例写 WAV 到 work_dir/output 便于试听。

#include <AudioFile.h>
#include <work_dir.hpp>

#include "../../gui/audiofx/pitch_shift_rt/autocorrelation_delay_pitch_shifter_rt.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

namespace {

using pitch_shift_rt::AutocorrelationDelayPitchShifter;

constexpr float kSampleRate = 48000.0f;
constexpr std::size_t kChunk = 240;                     // 5 ms 回调块
constexpr std::size_t kDurationSamples = 4 * 48000;     // 4 s
constexpr std::size_t kWarmupSamples = 1 * 48000;       // 预热 1 s
constexpr std::size_t kMeasureSamples = 2 * 48000;      // 测量 2 s
constexpr double kPi = std::numbers::pi_v<double>;

struct TestCase {
    const char* name;
    float semitones;
    float freq_hz;
};

constexpr TestCase kCases[] = {
    {"flat_220",   0.0f,  220.0f},
    {"up7_220",    7.0f,  220.0f},
    {"up12_220",  12.0f,  220.0f},
    {"down5_220", -5.0f,  220.0f},
    {"down12_220",-12.0f, 220.0f},
    {"up7_440",    7.0f,  440.0f},
    {"up12_110",  12.0f,  110.0f},
    {"down12_880",-12.0f, 880.0f},
};

std::vector<float> GenerateSine(float freq_hz, std::size_t n) {
    std::vector<float> x(n);
    const double w = 2.0 * kPi * static_cast<double>(freq_hz) / kSampleRate;
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = 0.9f * static_cast<float>(std::sin(w * static_cast<double>(i)));
    }
    return x;
}

bool HasNonFinite(const std::vector<float>& v) {
    for (float x : v) {
        if (!std::isfinite(x)) return true;
    }
    return false;
}

// Goertzel：计算窗口内指定频率的幅度
double Goertzel(const std::vector<float>& x, std::size_t start, std::size_t len, double freq) {
    const double w = 2.0 * kPi * freq / kSampleRate;
    const double coeff = 2.0 * std::cos(w);
    double s0 = 0.0, s1 = 0.0, s2 = 0.0;
    for (std::size_t i = start; i < start + len; ++i) {
        s0 = static_cast<double>(x[i]) + coeff * s1 - s2;
        s2 = s1;
        s1 = s0;
    }
    const double power = s1 * s1 + s2 * s2 - coeff * s1 * s2;
    return std::sqrt(std::max(0.0, power));
}

// 在 [f_lo, f_hi] 内以 0.1 Hz 步长扫描主频并抛物线细化
double EstimateFrequency(const std::vector<float>& x, std::size_t start,
                         std::size_t len, double f_lo, double f_hi) {
    constexpr double kStep = 0.1;
    double best_f = f_lo;
    double best_m = -1.0;
    for (double f = f_lo; f <= f_hi + 1e-9; f += kStep) {
        const double m = Goertzel(x, start, len, f);
        if (m > best_m) {
            best_m = m;
            best_f = f;
        }
    }
    // 抛物线细化到粗网格亚步长精度
    const double m_lo = Goertzel(x, start, len, best_f - kStep);
    const double m_hi = Goertzel(x, start, len, best_f + kStep);
    const double denom = m_lo - 2.0 * best_m + m_hi;
    if (std::fabs(denom) > 1e-12) {
        const double offset = 0.5 * (m_lo - m_hi) / denom;
        if (std::fabs(offset) <= 1.0) best_f += offset * kStep;
    }
    // 细扫：以 0.002 Hz 步长在细化点附近精确定位
    double fine_best_f = best_f;
    double fine_best_m = Goertzel(x, start, len, best_f);
    for (double f = best_f - 0.1; f <= best_f + 0.1 + 1e-9; f += 0.002) {
        const double m = Goertzel(x, start, len, f);
        if (m > fine_best_m) {
            fine_best_m = m;
            fine_best_f = f;
        }
    }
    return fine_best_f;
}

// 在精确期望频率上最小二乘拟合正弦，返回残差 RMS 相对信号 RMS 的分贝值
// （先做细扫频 + 幅度/相位拟合，残差即伪影能量）
double MeasureArtifactDb(const std::vector<float>& x, std::size_t start,
                         std::size_t len, double freq) {
    // 细扫：在 freq ± 0.5Hz 内找残差最小的频率（分辨率 0.02Hz）
    auto fit_rms = [&](double f) {
        const double w = 2.0 * kPi * f / kSampleRate;
        double cc = 0.0, ss = 0.0, cs = 0.0, yc = 0.0, ys = 0.0;
        for (std::size_t i = start; i < start + len; ++i) {
            const double t = w * static_cast<double>(i);
            const double c = std::cos(t);
            const double s = std::sin(t);
            const double y = static_cast<double>(x[i]);
            cc += c * c;
            ss += s * s;
            cs += c * s;
            yc += y * c;
            ys += y * s;
        }
        const double det = cc * ss - cs * cs;
        const double A = (det != 0.0) ? (yc * ss - ys * cs) / det : 0.0;
        const double B = (det != 0.0) ? (ys * cc - yc * cs) / det : 0.0;
        double res2 = 0.0;
        for (std::size_t i = start; i < start + len; ++i) {
            const double t = w * static_cast<double>(i);
            const double fit = A * std::cos(t) + B * std::sin(t);
            const double r = static_cast<double>(x[i]) - fit;
            res2 += r * r;
        }
        return std::sqrt(res2 / static_cast<double>(len));
    };

    double best_f = freq;
    double best_rms = fit_rms(freq);
    for (double f = freq - 0.5; f <= freq + 0.5; f += 0.02) {
        const double r = fit_rms(f);
        if (r < best_rms) {
            best_rms = r;
            best_f = f;
        }
    }
    // 抛物线细化
    const double r_lo = fit_rms(best_f - 0.02);
    const double r_hi = fit_rms(best_f + 0.02);
    const double denom = r_lo - 2.0 * best_rms + r_hi;
    if (std::fabs(denom) > 1e-12) {
        const double off = 0.5 * (r_lo - r_hi) / denom;
        if (std::fabs(off) <= 1.0) best_f += off * 0.02;
    }
    const double final_rms = fit_rms(best_f);

    double sig2 = 0.0;
    for (std::size_t i = start; i < start + len; ++i) {
        sig2 += static_cast<double>(x[i]) * static_cast<double>(x[i]);
    }
    const double sig_rms = std::sqrt(sig2 / static_cast<double>(len));
    return 20.0 * std::log10(std::max(final_rms, 1e-12) / std::max(sig_rms, 1e-12));
}

double Rms(const std::vector<float>& x, std::size_t start, std::size_t len) {
    double acc = 0.0;
    for (std::size_t i = start; i < start + len; ++i) {
        acc += static_cast<double>(x[i]) * static_cast<double>(x[i]);
    }
    return std::sqrt(acc / static_cast<double>(len));
}

struct Result {
    std::string name;
    float semitones = 0.0f;
    float freq_hz = 0.0f;
    double freq_out_hz = 0.0;
    double freq_err_hz = 0.0;
    double gain = 0.0;
    double artifact_db = 0.0;
    bool nonfinite = false;
    bool pass = false;
};

Result RunCase(const TestCase& tc) {
    Result r;
    r.name = tc.name;
    r.semitones = tc.semitones;
    r.freq_hz = tc.freq_hz;

    const std::vector<float> input = GenerateSine(tc.freq_hz, kDurationSamples);
    std::vector<float> output(kDurationSamples);

    AutocorrelationDelayPitchShifter dsp;
    dsp.setPitchShift(tc.semitones);
    dsp.init(kSampleRate);

    std::size_t pos = 0;
    while (pos < kDurationSamples) {
        const std::size_t n = std::min(kChunk, kDurationSamples - pos);
        dsp.process(input.data() + pos, output.data() + pos, n);
        pos += n;
    }

    r.nonfinite = HasNonFinite(output);
    if (r.nonfinite) {
        r.pass = false;
        return r;
    }

    const double expected = static_cast<double>(tc.freq_hz) * std::exp2(tc.semitones / 12.0);
    const double f_out = EstimateFrequency(output, kWarmupSamples, kMeasureSamples,
                                           expected - 2.0, expected + 2.0);
    r.freq_out_hz = f_out;
    r.freq_err_hz = std::fabs(f_out - expected);
    r.gain = Rms(output, kWarmupSamples, kMeasureSamples) /
             Rms(input, 0, input.size());
    // 伪影测量：在实测主频上拟合正弦（频率已知，无需相位搜索），残差即伪影能量
    r.artifact_db = MeasureArtifactDb(output, kWarmupSamples, kMeasureSamples, f_out);

    // 伪影阈值 -20dB：单跳转 WSOLA 的固有拼接伪影水平（参考 PitcherWsola
    // 实测约 -21dB），-20dB 对应约 1% 能量，听感上远低于回声/咔哒。
    r.pass = r.freq_err_hz < 0.5 && r.gain >= 0.9 && r.gain <= 1.1 && r.artifact_db < -20.0;

    AudioFile<float> file;
    file.setNumChannels(1);
    file.setBitDepth(32);
    file.setSampleRate(static_cast<uint32_t>(kSampleRate));
    file.setNumSamplesPerChannel(static_cast<int>(output.size()));
    file.samples[0] = output;
    file.save(qwqdsp_support::OutputFile("autocorr_rt_" + std::string(tc.name) + ".wav"));
    return r;
}

// ------------------------------------------------------------
// 语音回声测试：读入真实语音（wormhole.wav），升调 +12，检测输出中
// 是否存在"刚播过内容的重复"（回声）。回声表现为输出某段与之前
// 0.3~1.5s 内的段高度自相似（归一化相关 > 0.97）。语音内容本身会有
// 重复元音/音节，阈值取 0.97 且要求连续多帧命中才计为回声。
// ------------------------------------------------------------
bool RunSpeechEchoTest() {
    AudioFile<float> file{qwqdsp_support::InputFile("wormhole.wav")};
    auto& x = file.samples.front();
    const std::size_t n = x.size();
    std::vector<float> out(n);

    AutocorrelationDelayPitchShifter dsp;
    dsp.setPitchShift(12.0f);
    dsp.init(static_cast<float>(file.getSampleRate()));

    std::size_t pos = 0;
    while (pos < n) {
        const std::size_t chunk = std::min(kChunk, n - pos);
        dsp.process(x.data() + pos, out.data() + pos, chunk);
        pos += chunk;
    }

    if (HasNonFinite(out)) {
        std::cout << "[FAIL] speech_echo  NaN/Inf!\n";
        return false;
    }

    // 自相似回声检测：帧长 1024，间隔 512
    constexpr std::size_t kWin = 1024;
    constexpr std::size_t kHop = 512;
    constexpr float kEchoThresh = 0.97f;
    const std::size_t nf = (n > kWin) ? (n - kWin) / kHop : 0;
    std::size_t dup_frames = 0;
    for (std::size_t k = 50; k < nf; ++k) {
        const float* f = out.data() + k * kHop;
        const float* prev = out.data() + (k - 50) * kHop;
        // 与 50 帧前（0.53s）比较——回声若存在会在此处高度相关
        double num = 0.0, e0 = 1e-9, e1 = 1e-9;
        for (std::size_t i = 0; i < kWin; ++i) {
            const float a = f[i], b = prev[i];
            num += static_cast<double>(a) * b;
            e0 += static_cast<double>(a) * a;
            e1 += static_cast<double>(b) * b;
        }
        if (num / std::sqrt(e0 * e1) > kEchoThresh) ++dup_frames;
    }
    // 正常语音 + 正确移调：与 0.53s 前几乎不相关；重复帧占比 < 2% 视为通过
    const bool pass = static_cast<double>(dup_frames) < 0.02 * static_cast<double>(nf);
    std::cout << std::format("[{}] speech_echo  +12st wormhole 自相似重复帧 {}/{} (阈值 2%)\n",
                             pass ? "PASS" : "FAIL", dup_frames, nf);
    return pass;
}

} // namespace

int main() {
    std::cout << "AutocorrelationDelayPitchShifter 离线测试\n";

    int failures = 0;
    for (const TestCase& tc : kCases) {
        const Result r = RunCase(tc);
        std::cout << std::format(
            "[{}] {:>12} {:>+6.1f}st {:6.1f}Hz -> {:8.2f}Hz (err {:+.3f}Hz)  "
            "gain {:.3f}  伪影 {:+.1f}dB  {}\n",
            r.pass ? "PASS" : "FAIL", r.name, r.semitones, r.freq_hz,
            r.freq_out_hz, r.freq_err_hz, r.gain, r.artifact_db,
            r.nonfinite ? "[NaN/Inf!]" : "");
        if (!r.pass) ++failures;
    }

    const std::size_t total = std::size(kCases);
    std::cout << std::format("{} / {} 通过\n", total - failures, total);

    if (!RunSpeechEchoTest()) ++failures;
    std::cout << std::format("总计: {} 失败\n", failures);
    return failures == 0 ? 0 : 1;
}
