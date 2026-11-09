#pragma once

#include <algorithm>
#include <array>
#include <cstddef>

#include "maximizer_dsp.hpp"
#include "qwqdsp/filter/linkwitz_riley.hpp"

namespace multiband_maximizer_demo {

constexpr size_t kNumBands = 3;

// 经典三分频：120 Hz / 2.5 kHz
constexpr std::array<float, kNumBands - 1> kCrossoverHz{
    120.0f,
    2500.0f,
};

// ------------------------------------------------------------
//  ThreeBandSplitter — 立体声三段 LR4 分频
// ------------------------------------------------------------
class ThreeBandSplitter {
public:
    void Init(float sample_rate) noexcept {
        for (size_t i = 0; i < filter_l_.size(); ++i) {
            filter_l_[i].SetCutoff(kCrossoverHz[i], sample_rate);
            filter_r_[i].SetCutoff(kCrossoverHz[i], sample_rate);
        }
        Reset();
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
        auto const [low_l, high_l] = filter_l_[0].Tick(l);
        auto const [low_r, high_r] = filter_r_[0].Tick(r);
        band_l[0] = low_l;
        band_r[0] = low_r;

        auto const [mid_l, high2_l] = filter_l_[1].Tick(high_l);
        auto const [mid_r, high2_r] = filter_r_[1].Tick(high_r);
        band_l[1] = mid_l;
        band_r[1] = mid_r;
        band_l[2] = high2_l;
        band_r[2] = high2_r;
    }

private:
    std::array<qwqdsp_filter::LinkwitzRiley4, kNumBands - 1> filter_l_;
    std::array<qwqdsp_filter::LinkwitzRiley4, kNumBands - 1> filter_r_;
};

// ------------------------------------------------------------
//  MultibandMaximizer — 三段独立 maximizer + 末端 maximizer
//
//  每段：Drive -> look-ahead limiter -> 段 Ceiling
//  三段相加后，再经过一个最终 look-ahead limiter 保证总输出
//  不超过 Output Ceiling。
//
//  每个 maximizer 都带 look-ahead，所以总延迟约为 2 x look-ahead。
// ------------------------------------------------------------
class MultibandMaximizer {
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

    void SetParameters(std::array<float, kNumBands> const& drives,
                       std::array<float, kNumBands> const& ceilings,
                       float release_ms, float lookahead_ms,
                       float output_ceiling_db) noexcept {
        for (size_t band = 0; band < kNumBands; ++band) {
            band_maximizers_[band].SetParameters(
                drives[band], ceilings[band], release_ms, lookahead_ms);
        }
        final_maximizer_.SetParameters(0.0f, output_ceiling_db, release_ms, lookahead_ms);
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
            buffer[2 * i] = band_buffer_[0][2 * i] + band_buffer_[1][2 * i] + band_buffer_[2][2 * i];
            buffer[2 * i + 1] =
                band_buffer_[0][2 * i + 1] + band_buffer_[1][2 * i + 1] + band_buffer_[2][2 * i + 1];
        }

        final_maximizer_.Process(buffer, frame_count);

        for (size_t band = 0; band < kNumBands; ++band) {
            band_gain_reduction_db_[band] = band_maximizers_[band].GetGainReductionDb();
        }
        final_gain_reduction_db_ = final_maximizer_.GetGainReductionDb();
    }

    ThreeBandSplitter splitter_;
    std::array<maximizer_demo::StereoMaximizer, kNumBands> band_maximizers_;
    maximizer_demo::StereoMaximizer final_maximizer_;
    std::array<std::array<float, kChunkFrames * 2>, kNumBands> band_buffer_{};
    std::array<float, kNumBands> band_gain_reduction_db_{};
    float final_gain_reduction_db_{};
};

} // namespace multiband_maximizer_demo
