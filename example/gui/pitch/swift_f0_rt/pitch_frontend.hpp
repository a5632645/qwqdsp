#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <span>
#include <thread>
#include <vector>

#include <qwqdsp/pitch/hide/swift_f0_model.hpp>
#include <qwqdsp/spectral/real_fft_adv.hpp>
#include <qwqdsp/window/hann.hpp>

#include "inference_fast.hpp"

namespace swift_f0_rt {

// ------------------------------------------------------------
// 常量
// ------------------------------------------------------------

/// 模型帧长 / 跳长 / 每帧 mel bin 数
inline constexpr int kFrameLen = qwqdsp_swift_f0::kNFFT;
inline constexpr int kHop = qwqdsp_swift_f0::kHopLength;
inline constexpr int kNumMelBins = qwqdsp_swift_f0::kNumMelBins;
/// 模型前端的起始零填充（与训练时一致）
inline constexpr int kStftPad = qwqdsp_swift_f0::kSTFTPadding;
/// 卷积感受野半径（5 层 5x5 SAME），输出帧需左右各 kRecep 帧上下文
inline constexpr int kRecep = 10;
/// 每次推理的输出帧数（块处理摊薄推理开销）
inline constexpr int kBlock = 4;

// ------------------------------------------------------------
// LogMagStft
// ------------------------------------------------------------
// 16kHz → STFT(hop=256, n_fft=1024, Hann) → mag → slice[3:135) → log(mag+1e-8)
// 复刻模型前端，起始 384 零填充。逐块流式，每凑满一帧回调一次。
// ------------------------------------------------------------

class LogMagStft {
public:
    void Init() noexcept {
        fft_.Init(kFrameLen);
        buf_.assign(kFrameLen, 0.0f);
        count_ = kStftPad; // 起始零填充已“在缓冲里”
        consumed_ = 0;
        gain_.resize(fft_.NumBins());
        frame_out_.resize(kNumMelBins);
    }

    void Reset() noexcept {
        std::fill(buf_.begin(), buf_.end(), 0.0f);
        count_ = kStftPad;
        consumed_ = 0;
    }

    /**
     * @brief 输入 16kHz 样本，每产出一帧调用 on_frame(mel_bins, center_16k)
     * @param center_16k 该帧分析窗中心对应的 16kHz 原始时间轴样本位置（可为负）
     */
    template <typename OnFrame>
    void Process(std::span<const float> in, OnFrame&& on_frame) {
        size_t i = 0;
        while (i < in.size()) {
            size_t const need = static_cast<size_t>(kFrameLen) - static_cast<size_t>(count_);
            size_t const take = std::min(need, in.size() - i);
            std::copy_n(in.begin() + static_cast<std::ptrdiff_t>(i), take,
                        buf_.begin() + static_cast<std::ptrdiff_t>(count_));
            count_ += static_cast<int>(take);
            consumed_ += static_cast<int64_t>(take);
            i += take;

            if (count_ == kFrameLen) {
                // 加窗
                for (int n = 0; n < kFrameLen; ++n) {
                    work_[n] = buf_[n] * Window(n);
                }
                fft_.FFTGainPhase(work_, gain_);
                for (int f = 0; f < kNumMelBins; ++f) {
                    frame_out_[f] = std::log(gain_[f + qwqdsp_swift_f0::kSliceStart] + qwqdsp_swift_f0::kEps);
                }
                int64_t const center_16k = consumed_ - kFrameLen / 2;
                on_frame(std::span<const float>{frame_out_.data(), frame_out_.size()}, center_16k);

                // 左移 hop
                std::move(buf_.begin() + kHop, buf_.end(), buf_.begin());
                count_ -= kHop;
            }
        }
    }

private:
    /// Hann 窗（周期性，与分析一致）
    static float Window(int n) noexcept {
        constexpr float kTwoPi = 6.28318530717958647692f;
        float const t = static_cast<float>(n) / static_cast<float>(kFrameLen);
        return 0.5f * (1.0f - std::cos(kTwoPi * t));
    }

    qwqdsp_spectral::RealFftAdv fft_;
    std::vector<float> buf_;
    std::array<float, kFrameLen> work_{};
    std::vector<float> gain_;
    std::vector<float> frame_out_;
    int count_{};
    int64_t consumed_{};
};

// ------------------------------------------------------------
// FrameRing — 单生产者单消费者，传递 log-mag 帧
// ------------------------------------------------------------

struct Frame {
    std::array<float, kNumMelBins> lm{};
    int64_t center_sample{}; ///< 分析窗中心在原始采样率时间轴上的样本位置
    uint64_t index{};        ///< 全局帧序号
};

class FrameRing {
public:
    static constexpr int kCapacity = 1024;

    /// @brief 音频线程调用
    void Push(const Frame& f) noexcept {
        uint64_t const w = write_.load(std::memory_order_relaxed);
        slots_[w % kCapacity] = f;
        write_.store(w + 1, std::memory_order_release);
    }

    /// @brief 工作线程调用；成功取出则返回 true
    bool Pop(Frame& f) noexcept {
        uint64_t const r = read_.load(std::memory_order_relaxed);
        if (r == write_.load(std::memory_order_acquire)) {
            return false;
        }
        f = slots_[r % kCapacity];
        read_.store(r + 1, std::memory_order_release);
        return true;
    }

    int Depth() const noexcept {
        uint64_t const w = write_.load(std::memory_order_acquire);
        uint64_t const r = read_.load(std::memory_order_relaxed);
        return static_cast<int>(w - r);
    }

private:
    std::array<Frame, kCapacity> slots_{};
    std::atomic<uint64_t> write_{0};
    std::atomic<uint64_t> read_{0};
};

// ------------------------------------------------------------
// ResultRing — 工作线程 → GUI 线程，传递逐帧基频结果
// ------------------------------------------------------------

struct PitchResult {
    int64_t center_sample{}; ///< 分析窗中心在原始采样率时间轴上的样本位置
    float pitch_hz{};
    float confidence{};
};

class ResultRing {
public:
    static constexpr int kCapacity = 4096;

    void Push(const PitchResult& r) noexcept {
        uint64_t const w = write_.load(std::memory_order_relaxed);
        slots_[w % kCapacity] = r;
        write_.store(w + 1, std::memory_order_release);
    }

    bool Pop(PitchResult& r) noexcept {
        uint64_t const rd = read_.load(std::memory_order_relaxed);
        if (rd == write_.load(std::memory_order_acquire)) {
            return false;
        }
        r = slots_[rd % kCapacity];
        read_.store(rd + 1, std::memory_order_release);
        return true;
    }

private:
    std::array<PitchResult, kCapacity> slots_{};
    std::atomic<uint64_t> write_{0};
    std::atomic<uint64_t> read_{0};
};

// ------------------------------------------------------------
// SwiftF0Worker
// ------------------------------------------------------------
// 在独立线程消费 log-mag 帧，维护滑动历史，按块（kBlock 帧）调用推理，
// 每个输出帧使用左右各 kRecep 帧上下文，与整段批推理数值等价。
// 推理开销 ~14ms/块，块周期 kBlock*16ms，留出充足余量。
// ------------------------------------------------------------

class SwiftF0Worker {
public:
    void Start(FrameRing* in, ResultRing* out) {
        in_ = in;
        out_ = out;
        inference_.Init();
        running_.store(true, std::memory_order_release);
        thread_ = std::thread([this] { Run(); });
    }

    void Stop() noexcept {
        running_.store(false, std::memory_order_release);
        if (thread_.joinable()) {
            thread_.join();
        }
    }

    /// 最近一次块推理耗时（毫秒），供 GUI 显示
    float LastInferenceMs() const noexcept {
        return last_ms_.load(std::memory_order_relaxed);
    }

    /// 当前已发布到的最新帧序号
    uint64_t OutputFrame() const noexcept {
        return out_frame_.load(std::memory_order_relaxed);
    }

private:
    void Run() {
        std::vector<float> window;   // 连续窗口缓冲 [W * 132]
        std::vector<float> pitch, conf;

        while (running_.load(std::memory_order_acquire)) {
            // ---- 拉取所有可用帧 ----
            Frame f;
            while (in_->Pop(f)) {
                if (hist_.empty()) {
                    out_idx_ = f.index; // 对齐到首个可用帧
                }
                else if (f.index != hist_.back().index + 1) {
                    // 环形缓冲被覆盖导致序号跳变，重置后重新对齐
                    hist_.clear();
                    out_idx_ = f.index;
                }
                hist_.push_back(f);
            }

            // ---- 攒够一个块 + 右侧上下文就推理 ----
            bool did_work = false;
            while (true) {
                uint64_t const need_back = out_idx_ + kBlock - 1 + kRecep;
                if (hist_.empty() || hist_.back().index < need_back) {
                    break;
                }
                uint64_t lo = (out_idx_ > kRecep) ? (out_idx_ - kRecep) : 0;
                // 历史起点可能因环形缓冲跳变而晚于理论左边界
                lo = std::max(lo, hist_.front().index);
                uint64_t const hi = need_back;
                int const w = static_cast<int>(hi - lo + 1);

                // 拷到连续缓冲（hist_ 起始全局帧号）
                int const offset = static_cast<int>(lo - hist_.front().index);
                int const out_offset = static_cast<int>(out_idx_ - lo);
                window.resize(static_cast<size_t>(w) * kNumMelBins);
                for (int i = 0; i < w; ++i) {
                    auto const& src = hist_[static_cast<size_t>(offset + i)].lm;
                    std::copy(src.begin(), src.end(), window.begin() + static_cast<std::ptrdiff_t>(i) * kNumMelBins);
                }

                pitch.assign(static_cast<size_t>(w), 0.0f);
                conf.assign(static_cast<size_t>(w), 0.0f);
                auto const t0 = std::chrono::steady_clock::now();
                inference_.Process(window.data(), w, pitch.data(), conf.data());
                auto const t1 = std::chrono::steady_clock::now();
                last_ms_.store(static_cast<float>(std::chrono::duration<double, std::milli>(t1 - t0).count()),
                               std::memory_order_relaxed);

                // 发布本块 kBlock 帧
                for (int i = 0; i < kBlock; ++i) {
                    auto const& src = hist_[static_cast<size_t>(offset + out_offset + i)];
                    PitchResult r;
                    r.center_sample = src.center_sample;
                    r.pitch_hz = pitch[static_cast<size_t>(out_offset + i)];
                    r.confidence = conf[static_cast<size_t>(out_offset + i)];
                    out_->Push(r);
                }
                out_idx_ += kBlock;
                out_frame_.store(out_idx_, std::memory_order_relaxed);
                did_work = true;

                // 丢弃不再需要的过去帧
                while (!hist_.empty() && hist_.front().index + kRecep < out_idx_) {
                    hist_.pop_front();
                }
            }

            if (!did_work) {
                std::this_thread::sleep_for(std::chrono::milliseconds(2));
            }
        }
    }

    FrameRing* in_{};
    ResultRing* out_{};
    FastInference inference_;
    std::deque<Frame> hist_;
    uint64_t out_idx_{0};
    std::atomic<bool> running_{false};
    std::atomic<float> last_ms_{0.0f};
    std::atomic<uint64_t> out_frame_{0};
    std::thread thread_;
};

} // namespace swift_f0_rt
