#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

#include "../basicfx/maximizer_dsp.hpp"
#include "qwqdsp/filter/biquad.hpp"
#include "qwqdsp/filter/rbj.hpp"

namespace erb_multiband_maximizer_demo {

constexpr size_t kNumBands = 12;
constexpr size_t kNumCrossovers = kNumBands - 1;
constexpr float kErbLowHz = 20.0f;
constexpr float kErbHighHz = 20000.0f;

/// @brief ERB-rate: 21.4 * log10(1 + 0.00437 f)
static inline float HzToErb(float hz) noexcept {
    return 21.4f * std::log10(1.0f + 0.00437f * hz);
}

static inline float ErbToHz(float erb) noexcept {
    return (std::pow(10.0f, erb / 21.4f) - 1.0f) / 0.00437f;
}

// ------------------------------------------------------------
//  Highpass4 — 四阶 Butterworth 高通 (两个二阶节)
// ------------------------------------------------------------
class Highpass4 {
public:
    void SetCutoff(float f, float fs) noexcept {
        qwqdsp_filter::RBJ design;
        design.Highpass(qwqdsp::convert::Freq2W(f, fs), std::numbers::sqrt2_v<float> * 0.5f);
        filter1_.Set(design.b0, design.b1, design.b2, design.a1, design.a2);
        filter2_.Set(design.b0, design.b1, design.b2, design.a1, design.a2);
    }

    void Reset() noexcept {
        filter1_.Reset();
        filter2_.Reset();
    }

    float Tick(float x) noexcept {
        return filter2_.Tick(filter1_.Tick(x));
    }
private:
    qwqdsp_filter::Biquad filter1_;
    qwqdsp_filter::Biquad filter2_;
};

// ------------------------------------------------------------
//  ErbBandSplitter — 12 段 ERB 频程分频器
//
//  11 个高通分频点等间距落在 ERB 尺度上 (20 Hz - 20 kHz)。
//  相邻高通输出相减得到各频段: band_i = HP_i - HP_{i+1}
//  (HP_0 = 输入, HP_12 = 0)，因此 12 段相加逐样本等于原信号，
//  没有串联 LR 分频在密集频段上的相位抵消。
// ------------------------------------------------------------
class ErbBandSplitter {
public:
    using Edges = std::array<float, kNumBands + 1>;

    void Init(float sample_rate) noexcept {
        float const erb_low = HzToErb(kErbLowHz);
        float const erb_high = HzToErb(kErbHighHz);
        edges_[0] = kErbLowHz;
        edges_[kNumBands] = kErbHighHz;
        for (size_t i = 0; i < kNumCrossovers; ++i) {
            float const t = static_cast<float>(i + 1) / static_cast<float>(kNumBands);
            edges_[i + 1] = ErbToHz(erb_low + (erb_high - erb_low) * t);
            filter_l_[i].SetCutoff(edges_[i + 1], sample_rate);
            filter_r_[i].SetCutoff(edges_[i + 1], sample_rate);
        }
        Reset();
    }

    Edges const& GetEdges() const noexcept {
        return edges_;
    }

    void Reset() noexcept {
        for (auto& filter : filter_l_) {
            filter.Reset();
        }
        for (auto& filter : filter_r_) {
            filter.Reset();
        }
    }

    void Process(float l, float r, float* band_l, float* band_r) noexcept {
        float previous_l = l;
        float previous_r = r;
        for (size_t i = 0; i < kNumCrossovers; ++i) {
            float const highpass_l = filter_l_[i].Tick(l);
            float const highpass_r = filter_r_[i].Tick(r);
            band_l[i] = previous_l - highpass_l;
            band_r[i] = previous_r - highpass_r;
            previous_l = highpass_l;
            previous_r = highpass_r;
        }
        band_l[kNumBands - 1] = previous_l;
        band_r[kNumBands - 1] = previous_r;
    }
private:
    std::array<Highpass4, kNumCrossovers> filter_l_;
    std::array<Highpass4, kNumCrossovers> filter_r_;
    Edges edges_{};
};

// ------------------------------------------------------------
//  ErbMultibandMaximizer — 12 段独立最大化器 + 末端最大化器
//
//  每段: Drive -> look-ahead limiter -> 段 Ceiling
//  12 段相加后再过一个最终 limiter，保证总输出不超过 ceiling。
// ------------------------------------------------------------
class ErbMultibandMaximizer {
public:
    static constexpr size_t kChunkFrames = 512;

    void Init(float sample_rate) noexcept {
        splitter_.Init(sample_rate);
        for (auto& maximizer : band_maximizers_) {
            maximizer.Init(sample_rate);
        }
        final_maximizer_.Init(sample_rate);
        Reset();
    }

    void Reset() noexcept {
        splitter_.Reset();
        for (auto& maximizer : band_maximizers_) {
            maximizer.Reset();
        }
        final_maximizer_.Reset();
        band_gain_reduction_db_.fill(0.0f);
        final_gain_reduction_db_ = 0.0f;
    }

    void SetParameters(std::array<float, kNumBands> const& drives_db, float ceiling_db, float release_ms,
                       float lookahead_ms) noexcept {
        for (size_t band = 0; band < kNumBands; ++band) {
            band_maximizers_[band].SetParameters(drives_db[band], ceiling_db, release_ms, lookahead_ms);
        }
        final_maximizer_.SetParameters(0.0f, ceiling_db, release_ms, lookahead_ms);
    }

    void Process(float* buffer, size_t frame_count) noexcept {
        size_t offset = 0;
        while (offset < frame_count) {
            size_t const frames = std::min(kChunkFrames, frame_count - offset);
            ProcessChunk(buffer + offset * 2, frames);
            offset += frames;
        }
    }

    float GetBandGainReductionDb(size_t band) const noexcept {
        return band_gain_reduction_db_[band];
    }

    float GetFinalGainReductionDb() const noexcept {
        return final_gain_reduction_db_;
    }

    ErbBandSplitter::Edges const& GetBandEdges() const noexcept {
        return splitter_.GetEdges();
    }
private:
    void ProcessChunk(float* buffer, size_t frame_count) noexcept {
        for (size_t i = 0; i < frame_count; ++i) {
            float band_l[kNumBands];
            float band_r[kNumBands];
            splitter_.Process(buffer[2 * i], buffer[2 * i + 1], band_l, band_r);

            for (size_t band = 0; band < kNumBands; ++band) {
                band_buffer_[band][2 * i] = band_l[band];
                band_buffer_[band][2 * i + 1] = band_r[band];
            }
        }

        for (size_t band = 0; band < kNumBands; ++band) {
            band_maximizers_[band].Process(band_buffer_[band].data(), frame_count);
        }

        for (size_t i = 0; i < frame_count; ++i) {
            float sum_l = 0.0f;
            float sum_r = 0.0f;
            for (size_t band = 0; band < kNumBands; ++band) {
                sum_l += band_buffer_[band][2 * i];
                sum_r += band_buffer_[band][2 * i + 1];
            }
            buffer[2 * i] = sum_l;
            buffer[2 * i + 1] = sum_r;
        }

        final_maximizer_.Process(buffer, frame_count);

        for (size_t band = 0; band < kNumBands; ++band) {
            band_gain_reduction_db_[band] = band_maximizers_[band].GetGainReductionDb();
        }
        final_gain_reduction_db_ = final_maximizer_.GetGainReductionDb();
    }

    ErbBandSplitter splitter_;
    std::array<maximizer_demo::StereoMaximizer, kNumBands> band_maximizers_;
    maximizer_demo::StereoMaximizer final_maximizer_;
    std::array<std::array<float, kChunkFrames * 2>, kNumBands> band_buffer_{};
    std::array<float, kNumBands> band_gain_reduction_db_{};
    float final_gain_reduction_db_{};
};

} // namespace erb_multiband_maximizer_demo
