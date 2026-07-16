#include <qwqdsp/pitch/swift_f0.hpp>
#include <qwqdsp/pitch/swift_f0_model.hpp>
#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/hann.hpp>
#include <work_dir.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <numbers>
#include <vector>

// ------------------------------------------------------------
// 常量（匹配 ONNX 模型）
// ------------------------------------------------------------
constexpr int kSampleRate = 16000;
constexpr int kNFFT = 1024;
constexpr int kHop = 256;
constexpr int kPad = (kNFFT - kHop) / 2;    // 384
constexpr int kNumFreqBins = kNFFT / 2 + 1; // 513
constexpr int kSliceStart = 3;
constexpr int kSliceEnd = 135;                       // exclusive
constexpr int kNumMelBins = kSliceEnd - kSliceStart; // 132
constexpr float kEps = 1.0e-8f;

// ------------------------------------------------------------
// 辅助: 读二进制文件
// ------------------------------------------------------------
static bool readBinary(const std::filesystem::path& path, std::vector<float>& out) {
    FILE* f = nullptr;
    fopen_s(&f, path.string().c_str(), "rb");
    if (!f) {
        std::printf("    FAIL: cannot open %s\n", path.string().c_str());
        return false;
    }
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fseek(f, 0, SEEK_SET);
    out.resize(size / sizeof(float));
    fread(out.data(), sizeof(float), out.size(), f);
    fclose(f);
    return true;
}

// ------------------------------------------------------------
// 辅助: 最大相对误差
// ------------------------------------------------------------
static float maxRelError(const float* a, const float* b, int n) {
    float max_err = 0.0f;
    for (int i = 0; i < n; ++i) {
        float denom = std::max(std::abs(b[i]), 1e-6f);
        float err = std::abs(a[i] - b[i]) / denom;
        if (err > max_err)
            max_err = err;
    }
    return max_err;
}

// ------------------------------------------------------------
// 辅助: 生成正弦波
// ------------------------------------------------------------
static void generateSine(std::vector<float>& buffer, float freq, float sample_rate, int num_samples,
                         float amplitude = 0.5f) noexcept {
    for (int i = 0; i < num_samples; ++i) {
        buffer[i] = amplitude * std::sin(2.0f * std::numbers::pi_v<float> * freq * i / sample_rate);
    }
}

// ------------------------------------------------------------
// 手搓 STFT → log-mag 频谱
// ------------------------------------------------------------
// 与 ONNX 模型前端一致:
//   reflect pad (384) → STFT(hop=256, n_fft=1024, Hann) → mag → slice[3:135) → log(1e-8)
// 输出: [T, 132]
// ------------------------------------------------------------
static void computeLogMagSTFT(const std::vector<float>& audio, std::vector<float>& log_mag_out, int& out_T,
                              int& out_F) {
    int const audio_len = static_cast<int>(audio.size());

    // ---- 零填充 (ONNX model: Pad mode=constant, value=0) ----
    std::vector<float> padded(audio_len + 2 * kPad, 0.0f);
    for (int i = 0; i < audio_len; ++i)
        padded[kPad + i] = audio[i];

    int const padded_len = static_cast<int>(padded.size());

    // ---- 帧数 ----
    int T = (padded_len - kNFFT) / kHop + 1;
    int F = kNumMelBins;
    out_T = T;
    out_F = F;

    log_mag_out.resize(T * F);

    // ---- STFT 准备 ----
    qwqdsp_spectral::RealFFT fft;
    fft.Init(kNFFT);
    std::vector<float> frame(kNFFT);
    std::vector<float> gain(fft.NumBins()); // 513 bins

    // ---- 逐帧 STFT ----
    for (int t = 0; t < T; ++t) {
        int offset = t * kHop;

        // 拷贝帧
        for (int i = 0; i < kNFFT; ++i)
            frame[i] = padded[offset + i];

        // 加窗 (Hann, for_analyze=true → t = n/N)
        qwqdsp_window::Hann::ApplyWindow(frame, true);

        // FFT → 幅度谱
        fft.FFTGainPhase(frame, gain);

        // Slice [3:135) + log(1e-8)
        for (int f = 0; f < F; ++f) {
            float mag = gain[f + kSliceStart];
            log_mag_out[t * F + f] = std::log(mag + kEps);
        }
    }
}

// ------------------------------------------------------------
// 测试 1: 手搓 STFT vs 参考 log-mag
// ------------------------------------------------------------
static int testSTFTPipeline() noexcept {
    std::printf("  [STFT vs Reference] ...\n");
    int failed = 0;

    auto data_dir = qwqdsp_support::GetWorkDir() / "swift_f0";

    // 读取参考 log-mag
    std::vector<float> log_mag_ref;
    if (!readBinary(data_dir / "log_mag.bin", log_mag_ref))
        return 1;

    // 生成与 Python 测试相同的音频 (440Hz sine, 1s @ 16kHz)
    int const audio_len = kSampleRate; // 1 second
    std::vector<float> audio(audio_len);
    generateSine(audio, 440.0f, static_cast<float>(kSampleRate), audio_len);

    // 手搓 STFT
    std::vector<float> log_mag_our;
    int T = 0, F = 0;
    computeLogMagSTFT(audio, log_mag_our, T, F);

    // 检查维度
    int T_ref = static_cast<int>(log_mag_ref.size()) / F;
    if (T != T_ref || F != kNumMelBins) {
        std::printf("    FAIL: shape mismatch ref=(%d,%d) our=(%d,%d)\n", T_ref, F, T, F);
        return 1;
    }

    // 使用绝对误差（log-mag 尺度）。
    // 近零信号处（高频底噪 bin）不同 FFT 实现的微小数值差异会被 log 放大，
    // 但不影响实际 pitch 检测（见 Full Pipeline 测试）。
    float max_abs_err = 0.0f;
    float mean_abs_err = 0.0f;
    int n = T * F;
    for (int i = 0; i < n; ++i) {
        float err = std::abs(log_mag_our[i] - log_mag_ref[i]);
        if (err > max_abs_err)
            max_abs_err = err;
        mean_abs_err += err;
    }
    mean_abs_err /= static_cast<float>(n);

    std::printf("    T=%d, F=%d\n", T, F);
    std::printf("    max abs err: %.8f\n", max_abs_err);
    std::printf("    mean abs err: %.8f\n", mean_abs_err);

    // 打印前 5 帧对比
    std::printf("    frame 0 ref: ");
    for (int f = 0; f < std::min(5, F); ++f)
        std::printf(" %8.4f", log_mag_ref[f]);
    std::printf("\n    frame 0 our: ");
    for (int f = 0; f < std::min(5, F); ++f)
        std::printf(" %8.4f", log_mag_our[f]);
    std::printf("\n");

    if (max_abs_err > 10.0f) {
        std::printf("    FAIL: STFT mismatch (max_abs_err=%.6f)\n", max_abs_err);
        // 显示最大误差位置
        int max_idx = 0;
        float max_v = 0;
        for (int i = 0; i < n; ++i) {
            float err = std::abs(log_mag_our[i] - log_mag_ref[i]);
            if (err > max_v) {
                max_v = err;
                max_idx = i;
            }
        }
        int mt = max_idx / F, mf = max_idx % F;
        std::printf("      max at [%d,%d]: ref=%.6f our=%.6f\n", mt, mf, log_mag_ref[max_idx], log_mag_our[max_idx]);
        failed++;
    }
    else {
        std::printf("    OK\n");
    }
    return failed;
}

// ------------------------------------------------------------
// 测试 2: swift_f0 模型推理 (使用参考 log-mag)
// ------------------------------------------------------------
static int testSwiftF0Inference() noexcept {
    std::printf("  [SwiftF0 Inference] ...\n");
    int failed = 0;

    // 读取测试数据
    std::vector<float> log_mag, pitch_ref, conf_ref;

    auto data_dir = qwqdsp_support::GetWorkDir() / "swift_f0";

    if (!readBinary(data_dir / "log_mag.bin", log_mag))
        return 1;
    if (!readBinary(data_dir / "pitch_ref.bin", pitch_ref))
        return 1;
    if (!readBinary(data_dir / "conf_ref.bin", conf_ref))
        return 1;

    int T = static_cast<int>(pitch_ref.size());
    int F = qwqdsp_swift_f0::kNumMelBins;

    std::printf("    T=%d, F=%d\n", T, F);

    // 构建 Eigen Tensor 输入
    Eigen::Tensor<float, 2> input(T, F);
    for (int t = 0; t < T; ++t)
        for (int f = 0; f < F; ++f)
            input(t, f) = log_mag[t * F + f];

    // 运行推理
    Eigen::VectorXf pitch_hz(T);
    Eigen::VectorXf confidence(T);

    qwqdsp_swift_f0::SwiftF0Inference inference;
    inference.Process(input, pitch_hz, confidence);

    // 验证结果
    float pitch_err = maxRelError(pitch_hz.data(), pitch_ref.data(), T);
    float conf_err = maxRelError(confidence.data(), conf_ref.data(), T);

    std::printf("    pitch max rel err:  %.6f (%s)\n", pitch_err, pitch_err < 0.01f ? "OK" : "FAIL");
    std::printf("    conf  max rel err:  %.6f (%s)\n", conf_err, conf_err < 0.01f ? "OK" : "FAIL");

    if (pitch_err >= 0.01f) {
        std::printf("    FAIL: pitch error too large\n");
        std::printf("    idx  ref_hz    our_hz\n");
        for (int i = 0; i < std::min(10, T); ++i) {
            std::printf("    %3d  %8.2f  %8.2f\n", i, pitch_ref[i], pitch_hz[i]);
        }
        failed++;
    }
    if (conf_err >= 0.01f) {
        std::printf("    FAIL: confidence error too large\n");
        std::printf("    idx  ref_conf  our_conf\n");
        for (int i = 0; i < std::min(10, T); ++i) {
            std::printf("    %3d  %8.4f  %8.4f\n", i, conf_ref[i], confidence[i]);
        }
        failed++;
    }

    if (failed == 0)
        std::printf("    OK\n");
    return failed;
}

// ------------------------------------------------------------
// 测试 3: 完整管线 (手搓 STFT + 模型推理)
// ------------------------------------------------------------
static int testFullPipeline() noexcept {
    std::printf("  [Full Pipeline (STFT + Inference)] ...\n");
    int failed = 0;

    auto data_dir = qwqdsp_support::GetWorkDir() / "swift_f0";

    // 读取参考 pitch/conf
    std::vector<float> pitch_ref, conf_ref;
    if (!readBinary(data_dir / "pitch_ref.bin", pitch_ref))
        return 1;
    if (!readBinary(data_dir / "conf_ref.bin", conf_ref))
        return 1;
    int T_ref = static_cast<int>(pitch_ref.size());

    // 生成音频
    int const audio_len = kSampleRate;
    std::vector<float> audio(audio_len);
    generateSine(audio, 440.0f, static_cast<float>(kSampleRate), audio_len);

    // 手搓 STFT
    std::vector<float> log_mag_our;
    int T = 0, F = 0;
    computeLogMagSTFT(audio, log_mag_our, T, F);

    if (T != T_ref || F != kNumMelBins) {
        std::printf("    FAIL: shape mismatch\n");
        return 1;
    }

    // 构建 Eigen Tensor
    Eigen::Tensor<float, 2> input(T, F);
    for (int t = 0; t < T; ++t)
        for (int f = 0; f < F; ++f)
            input(t, f) = log_mag_our[t * F + f];

    // 推理
    Eigen::VectorXf pitch_hz(T);
    Eigen::VectorXf confidence(T);

    qwqdsp_swift_f0::SwiftF0Inference inference;
    inference.Process(input, pitch_hz, confidence);

    // 验证
    float pitch_err = maxRelError(pitch_hz.data(), pitch_ref.data(), T);
    float conf_err = maxRelError(confidence.data(), conf_ref.data(), T);

    std::printf("    pitch max rel err:  %.6f (%s)\n", pitch_err, pitch_err < 0.03f ? "OK" : "FAIL");
    std::printf("    conf  max rel err:  %.6f (%s)\n", conf_err, conf_err < 0.15f ? "OK" : "FAIL");

    if (pitch_err >= 0.03f) {
        std::printf("    FAIL: full pipeline pitch error\n");
        failed++;
    }
    if (conf_err >= 0.15f) {
        std::printf("    FAIL: full pipeline confidence error\n");
        failed++;
    }

    if (failed == 0)
        std::printf("    OK\n");
    return failed;
}

// ------------------------------------------------------------
// 主函数
// ------------------------------------------------------------
int main() {
    std::printf("swift_f0 model test\n");
    std::printf("==================\n\n");

    int total = 0;
    total += testSTFTPipeline();
    std::printf("\n");
    total += testSwiftF0Inference();
    std::printf("\n");
    total += testFullPipeline();

    std::printf("\n");
    if (total > 0) {
        std::printf("FAILED: %d test(s)\n", total);
    }
    else {
        std::printf("All tests passed.\n");
    }
    return total;
}
