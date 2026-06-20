#include "AudioFile.h"
#include "work_dir.hpp"
#include <cassert>
#include <iostream>
#include <span>
#include <string>
#include <vector>

#if 0

/**
 * @brief 多相抽取器 / Polyphase downsampler (decimator)
 *
 * 将多相分解与转置 FIR 结合，在降低采样率的同时提供抗混叠滤波。
 *
 * 使用流程:
 * @code
 *   PolyphaseDownsampler d;
 *   d.Init(ir, D);          // D: 下采样倍数
 *   d.Reset();
 *   auto out = d.Process(in);
 * @endcode
 *
 * @ref https://ccrma.stanford.edu/~jos/sasp/ Polyphase implementation
 */
class PolyphaseDownsampler {
public:
    /**
     * @brief 初始化多相滤波器组
     * @param coeff   原型低通 FIR 系数（长度任意）
     * @param downsample  下采样倍数 D
     *
     * 将原型滤波器系数按多相分解重新排列为 D 个相位子滤波器，
     * 每个子滤波器长度为 ceil(len(coeff) / D)。
     */
    void Init(std::span<float> coeff, int downsample) {
        downsample_ = downsample;

        int len_coeff = static_cast<int>(coeff.size());
        each_phase_len_ = len_coeff / downsample;
        if (len_coeff % downsample != 0) {
            ++each_phase_len_;
        }
        
        phase_coeffs_.resize(downsample * each_phase_len_);
        for (int i = 0; i < downsample; ++i) {
            auto dst_it = phase_coeffs_.begin() + i * each_phase_len_;
            for (int j = 0; j < each_phase_len_; ++j) {
                int idx = j * downsample + i;
                if (idx < coeff.size()) {
                    *dst_it = coeff[idx];
                }
                else {
                    *dst_it = 0;
                }
                ++dst_it;
            }
        }

        phase_state_.resize(phase_coeffs_.size());
    }

    /**
     * @brief 重置所有内部状态（相位累加器、延迟线）
     *
     * 在开始处理新的一段信号前调用，或在切换输入流时调用。
     */
    void Reset() noexcept {
        std::fill(phase_state_.begin(), phase_state_.end(), 0);
        current_phase_ = 0;
        acc_ = 0;
    }

    /**
     * @brief 逐块处理（块大小 == D），返回一个输出样本
     * @param x  长度为 D 的输入块
     * @return 下采样后的单个输出样本
     *
     * 调用方需保证 x.size() == downsample_。
     * 适合在已知输入已按 D 分块的场景下使用。
     */
    float Tick(std::span<float> x) noexcept {
        assert(x.size() == downsample_);

        float y = 0;
        for (int i = 0; i < downsample_; ++i) {
            auto* lag_it = phase_state_.data() + i * each_phase_len_;
            auto* coeff_it = phase_coeffs_.data() + i * each_phase_len_;
            y += TransposeFir(lag_it, coeff_it, x[i], each_phase_len_);
        }
        return y;
    }

    /**
     * @brief 逐样本处理（样本级接口），自动累积到满 D 个样本后输出
     * @param x  单个输入样本
     * @param y  输出指针的引用；当累积满 D 个输入样本时写入 *y 并前移 y
     *
     * 用法示例:
     * @code
     *   float out_buf[1024];
     *   float* y = out_buf;
     *   for (float s : input_samples) {
     *       downsampler.Tick(s, y);
     *   }
     *   size_t out_len = y - out_buf;
     * @endcode
     */
    void Tick(float x, float*& y) noexcept {
        auto* lag_it = phase_state_.data() + current_phase_ * each_phase_len_;
        auto* coeff_it = phase_coeffs_.data() + current_phase_ * each_phase_len_;
        acc_ += TransposeFir(lag_it, coeff_it, x, each_phase_len_);

        ++current_phase_;
        if (current_phase_ >= downsample_) {
            current_phase_ = 0;
            *y = acc_;
            ++y;
            acc_ = 0;
        }
    }

    /**
     * @brief 批量处理整个输入缓冲区，一次完成下采样
     * @param x  输入样本序列
     * @return 下采样后的输出序列，长度为 floor(x.size() / D) 或 ceil(x.size() / D)
     *
     * 内部自动按 D 分块，依次调用 Tick(span) 完成处理。
     * 最后不足 D 个的尾部样本同样会产出一个输出样本（标量模式）。
     *
     * @note 调用前需要先 Init() 并 Reset()。
     */
    std::vector<float> Process(std::span<float> x) {
        std::vector<float> out;

        int full_loop_count = x.size() / downsample_;
        int scalar_count = x.size() - full_loop_count * downsample_;
        out.reserve(scalar_count == 0 ? full_loop_count : full_loop_count + 1);

        for (int block = 0; block < full_loop_count; ++block) {
            int const block_start = block * downsample_;
            float y = 0;
            for (int j = 0; j < downsample_; ++j) {
                auto* lag_it = phase_state_.data() + j * each_phase_len_;
                auto* coeff_it = phase_coeffs_.data() + j * each_phase_len_;
                y += TransposeFir(lag_it, coeff_it, x[block_start + j], each_phase_len_);
            }
            out.push_back(y);
        }

        if (scalar_count != 0) {
            int const block_start = full_loop_count * downsample_;
            float y = 0;
            for (int j = 0; j < scalar_count; ++j) {
                auto* lag_it = phase_state_.data() + j * each_phase_len_;
                auto* coeff_it = phase_coeffs_.data() + j * each_phase_len_;
                y += TransposeFir(lag_it, coeff_it, x[block_start + j], each_phase_len_);
            }
            for (int j = scalar_count; j < downsample_; ++j) {
                auto* lag_it = phase_state_.data() + j * each_phase_len_;
                auto* coeff_it = phase_coeffs_.data() + j * each_phase_len_;
                y += TransposeFir(lag_it, coeff_it, 0, each_phase_len_);
            }
            out.push_back(y);
        }

        return out;
    }
private:
    /**
     * @brief 转置形式 FIR 滤波器（Transposed Direct-Form II）
     * @param lag    延迟线状态（长度 len）
     * @param coeff  FIR 系数（长度 len）
     * @param x      当前输入样本
     * @param len    滤波器阶数
     * @return 滤波后的输出样本
     *
     * 相比标准 FIR，转置形式在一个时钟周期内即可产出一个输出，
     * 适合流水线 / SIMD 友好的实现。
     */
    static float TransposeFir(float* lag, float const* coeff, float x, int len) noexcept {
        float const y = coeff[0] * x + lag[0];
        for (size_t i = 0; i < len - 1; ++i) {
            lag[i] = lag[i + 1] + coeff[i + 1] * x;
        }
        lag[len - 1] = coeff[len - 1] * x;
        return y;
    }

    /** 下采样倍数 D */
    int downsample_{};
    /** 每个相位子滤波器的长度 = ceil(len(coeff) / D) */
    int each_phase_len_{};
    /** 多相滤波器系数，按相位连续存储：phase_coeffs_[phase * each_phase_len_ + tap] */
    std::vector<float> phase_coeffs_;
    /** 各相位对应的延迟线状态 */
    std::vector<float> phase_state_;

    /** 当前相位索引 (0 .. D-1)，样本级 Tick 使用 */
    int current_phase_{};
    /** 跨相位累积的输出值，样本级 Tick 使用 */
    float acc_{};
};

#include <qwqdsp/filter/window_fir.hpp>
#include <qwqdsp/oscillator/vic_sine_osc.hpp>
#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/hann.hpp>


// int main() {
//     float x[1024]{};
//     float y[1024]{};

//     qwqdsp_oscillator::VicSineOsc osc;
//     osc.SetFreq(3, 1024);
//     for (int i = 0; i < 1024; ++i) {
//         x[i] = osc.Tick();
//     }

//     float ir[63];
//     qwqdsp_filter::WindowFIR::Lowpass(ir, std::numbers::pi_v<float> / 8);
//     qwqdsp_window::Hann::ApplyWindow(ir, false);

//     PolyphaseDownsampler downsample;
//     downsample.Init(ir, 8);
//     downsample.Reset();

//     int j = 0;
//     for (int i = 0; i < 1024; i+=8) {
//         y[j++] = downsample.Tick({&x[i], 8});
//     }

//     float g[1024/2 + 1];
//     qwqdsp_spectral::RealFFT fft;
//     fft.Init(1024);
//     fft.FFTGainPhase(y, g);
// }

int main() {
    AudioFile<float> file;
    file.load(qwqdsp_support::WormholeWav());

    int const downsample_factor = 8;
    int const filter_len = 63;

    // 设计抗混叠低通 FIR (截止频率 = pi / downsample_factor)
    std::vector<float> ir(filter_len);
    qwqdsp_filter::WindowFIR::Lowpass(ir, std::numbers::pi_v<float> / downsample_factor);
    qwqdsp_window::Hann::ApplyWindow(ir, false);

    // 初始化多相下采样器
    PolyphaseDownsampler downsampler;
    downsampler.Init(ir, downsample_factor);

    // 提取第一通道
    auto& left_channel = file.samples[0];
    float original_sample_rate = static_cast<float>(file.getSampleRate());

    // 下采样
    downsampler.Reset();
    std::vector<float> downsampled = downsampler.Process(left_channel);

    // 构建输出文件 (单声道)
    AudioFile<float> out_file;
    out_file.setNumChannels(1);
    out_file.setNumSamplesPerChannel(static_cast<int>(downsampled.size()));
    out_file.setSampleRate(static_cast<uint32_t>(original_sample_rate / downsample_factor));
    out_file.setBitDepth(file.getBitDepth());
    out_file.samples[0] = std::move(downsampled);

    // 保存
    std::string out_path = qwqdsp_support::OutputFile("wormhole_downsampled.wav");
    out_file.save(out_path);

    std::cout << "Saved " << out_path << " ("
              << original_sample_rate << " Hz -> "
              << out_file.getSampleRate() << " Hz, "
              << out_file.getNumSamplesPerChannel() << " samples)"
              << std::endl;

    return 0;
}

#else

#include "AudioFile.h"
#include <cassert>
#include <iostream>
#include <span>
#include <string>
#include <vector>

/**
 * @brief 多相插值器 / Polyphase upsampler (interpolator)
 *
 * 先对输入信号插零（零填充），再通过多相分解的 FIR 滤波器组
 * 完成抗镜像滤波，从而提升采样率。
 *
 * 相比直接先插零再滤波，多相结构将计算量降低到约 1/U。
 *
 * 使用流程:
 * @code
 *   PolyphaseUpsampler u;
 *   u.Init(ir, U);          // U: 上采样倍数
 *   u.Reset();
 *   auto out = u.Process(in);
 * @endcode
 *
 * @ref https://ccrma.stanford.edu/~jos/sasp/ Polyphase implementation
 */
class PolyphaseUpsampler {
public:
    /**
     * @brief 初始化多相滤波器组
     * @param coeff   原型低通 FIR 系数（长度任意）
     * @param upsample  上采样倍数 U
     *
     * 将原型滤波器系数按多相分解重新排列为 U 个相位子滤波器，
     * 每个子滤波器长度为 ceil(len(coeff) / U)。
     *
     * @note 系数自动乘以 U 以补偿插零带来的增益损失。
     */
    void Init(std::span<float> coeff, int upsample) {
        upsample_ = upsample;

        int len_coeff = static_cast<int>(coeff.size());
        each_phase_len_ = len_coeff / upsample;
        if (len_coeff % upsample != 0) {
            ++each_phase_len_;
        }

        phase_coeffs_.resize(upsample * each_phase_len_);
        for (int i = 0; i < upsample; ++i) {
            auto dst_it = phase_coeffs_.begin() + i * each_phase_len_;
            for (int j = 0; j < each_phase_len_; ++j) {
                int idx = j * upsample + i;
                if (idx < coeff.size()) {
                    *dst_it = coeff[idx] * upsample;
                }
                else {
                    *dst_it = 0;
                }
                ++dst_it;
            }
        }

        phase_state_.resize(each_phase_len_);
    }

    /**
     * @brief 重置延迟线状态和相位计数器
     *
     * 处理新信号段前调用，避免跨段交界处的瞬态。
     */
    void Reset() noexcept {
        std::fill(phase_state_.begin(), phase_state_.end(), 0);
        current_phase_ = 0;
    }

    /**
     * @brief 逐样本处理（块输出接口），一个输入样本产生 U 个输出
     * @param x  单个输入样本
     * @param y  输出 span，长度需为 U；写入后 y[i] = 第 i 个相位子滤波器的输出
     *
     * 内部先更新延迟线 (PushFirState)，再依次计算所有 U 个相位的输出。
     */
    void Tick(float x, std::span<float> y) noexcept {
        assert(y.size() == upsample_);

        PushFirState(x);
        for (int i = 0; i < upsample_; ++i) {
            auto* coeff_it = phase_coeffs_.data() + i * each_phase_len_;
            y[i] = Fir(coeff_it);
        }
    }

    /**
     * @brief 逐样本处理（样本级输出接口），每次调用返回一个输出样本
     * @param x  输入指针的引用；当消耗完一个输入样本时自动前移
     * @return 当前相位子滤波器输出的单个样本
     *
     * 用法示例:
     * @code
     *   float const* in = input.data();
     *   for (size_t i = 0; i < output.size(); ++i) {
     *       output[i] = upsampler.Tick(in);
     *   }
     * @endcode
     */
    float Tick(float*& x) noexcept {
        if (current_phase_ == 0) {
            PushFirState(*x);
        }

        auto* coeff_it = phase_coeffs_.data() + current_phase_ * each_phase_len_;
        float y = Fir(coeff_it);

        ++current_phase_;
        if (current_phase_ >= upsample_) {
            current_phase_ = 0;
            ++x;
        }

        return y;
    }

    /**
     * @brief 批量处理整个输入缓冲区，一次完成上采样
     * @param x  输入样本序列
     * @return 上采样后的输出序列，长度为 x.size() * U
     *
     * 内部对每个输入样本依次调用 PushFirState + 遍历 U 个相位。
     *
     * @note 调用前需要先 Init() 并 Reset()。
     */
    std::vector<float> Process(std::span<float> x) {
        std::vector<float> out;
        out.reserve(x.size() * upsample_);

        for (auto s : x) {
            PushFirState(s);
            for (int i = 0; i < upsample_; ++i) {
                auto* coeff_it = phase_coeffs_.data() + i * each_phase_len_;
                out.push_back(Fir(coeff_it));
            }
        }

        return out;
    }
private:
    /**
     * @brief 直接形式 FIR 内积
     * @param coeff  当前相位的滤波器系数
     * @return 滤波后的输出样本
     */
    float Fir(float const* coeff) noexcept {
        float sum = 0;
        for (int i = 0; i < each_phase_len_; ++i) {
            sum += phase_state_[i] * coeff[i];
        }
        return sum;
    }

    /**
     * @brief 将新样本推入延迟线（标准 Direct-Form 移位）
     * @param x  输入样本
     *
     * 每次收到新的输入样本时调用一次（无论上采样倍数是多少）。
     */
    void PushFirState(float x) noexcept {
        for (int i = each_phase_len_ - 1; i > 0; --i) {
            phase_state_[i] = phase_state_[i - 1];
        }
        phase_state_[0] = x;
    }

    /** 上采样倍数 U */
    int upsample_{};
    /** 每个相位子滤波器的长度 = ceil(len(coeff) / U) */
    int each_phase_len_{};
    /** 多相滤波器系数，按相位连续存储：phase_coeffs_[phase * each_phase_len_ + tap] */
    std::vector<float> phase_coeffs_;
    /** 直接形式 FIR 延迟线（长度 = each_phase_len_） */
    std::vector<float> phase_state_;

    /** 当前相位索引 (0 .. U-1)，样本级 Tick 使用 */
    int current_phase_{};
};

#include <qwqdsp/filter/window_fir.hpp>
#include <qwqdsp/oscillator/vic_sine_osc.hpp>
#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/hann.hpp>

// int main() {
//     float x[1024 / 8]{};
//     float y[1024]{};

//     qwqdsp_oscillator::VicSineOsc osc;
//     osc.SetFreq(3, 1024 / 8);
//     for (int i = 0; i < 1024 / 8; ++i) {
//         x[i] = osc.Tick();
//     }

//     float ir[63];
//     qwqdsp_filter::WindowFIR::Lowpass(ir, std::numbers::pi_v<float> / 8);
//     qwqdsp_window::Hann::ApplyWindow(ir, false);

//     PolyphaseUpsampler downsample;
//     downsample.Init(ir, 8);
//     downsample.Reset();

//     for (int i = 0; i < 1024 / 8; ++i) {
//         downsample.Tick(x[i], {&y[i * 8], 8});
//     }

//     float g[1024/2 + 1];
//     qwqdsp_spectral::RealFFT fft;
//     fft.Init(1024);
//     fft.FFTGainPhase(y, g);
// }

int main() {
    AudioFile<float> file;
    file.load(qwqdsp_support::InputFile("wormhole_downsampled.wav"));

    int const upsample_factor = 8;
    int const filter_len = 63;

    // 设计抗混叠低通 FIR (截止频率 = pi / upsample_factor)
    std::vector<float> ir(filter_len);
    qwqdsp_filter::WindowFIR::Lowpass(ir, std::numbers::pi_v<float> / upsample_factor);
    qwqdsp_window::Hann::ApplyWindow(ir, false);

    // 初始化多相上采样器
    PolyphaseUpsampler upsampler;
    upsampler.Init(ir, upsample_factor);

    // 提取第一通道
    auto& left_channel = file.samples[0];
    float original_sample_rate = static_cast<float>(file.getSampleRate());

    // 上采样
    upsampler.Reset();
    std::vector<float> upsampled = upsampler.Process(left_channel);

    // 构建输出文件 (单声道)
    AudioFile<float> out_file;
    out_file.setNumChannels(1);
    out_file.setNumSamplesPerChannel(static_cast<int>(upsampled.size()));
    out_file.setSampleRate(static_cast<uint32_t>(original_sample_rate * upsample_factor));
    out_file.setBitDepth(file.getBitDepth());
    out_file.samples[0] = std::move(upsampled);

    // 保存
    std::string out_path = qwqdsp_support::OutputFile("wormhole_upsampled.wav");
    out_file.save(out_path);

    std::cout << "Saved " << out_path << " (" << original_sample_rate << " Hz -> " << out_file.getSampleRate()
              << " Hz, " << out_file.getNumSamplesPerChannel() << " samples)" << std::endl;

    return 0;
}

#endif
