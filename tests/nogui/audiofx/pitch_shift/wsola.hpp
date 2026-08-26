#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <qwqdsp/segement/mono_reader.hpp>
#include <span>
#include <vector>

namespace qwqdsp_test {
namespace detail {

// ------------------------------------------------------------
// OlaBuffer
// ------------------------------------------------------------
/**
 * @brief 环形重叠相加（OLA）缓冲区，用于 WSOLA 合成输出。
 *
 * 缓冲区大小为 2 的幂，通过掩码实现环形环绕。
 */
class OlaBuffer {
public:
    void Init(size_t size) {
        size_ = size;
        mask_ = size - 1;
        buffer_.resize(size);
    }

    void Add(float* buffer, float* window, size_t count, size_t hop) {
        size_t pos = ola_add_pos_;
        for (size_t i = 0; i < count; ++i) {
            buffer_[pos++] += buffer[i] * window[i];
            pos &= mask_;
        }
        ola_add_end_ += hop;
        ola_add_pos_ += hop;
        size_ += hop;
        ola_add_end_ &= mask_;
        ola_add_pos_ &= mask_;
    }

    void Read(float* buffer, size_t count) {
        size_t read = std::min(size_, count);
        for (size_t i = 0; i < read; ++i) {
            buffer[i] = buffer_[rpos_];
            buffer_[rpos_] = 0;
            rpos_ = (rpos_ + 1) & mask_;
        }
        size_ -= read;
        std::fill_n(buffer + read, count - read, 0.0f);
    }

    void ReadPushback(std::vector<float>& buffer) {
        for (size_t i = 0; i < size_; ++i) {
            buffer.push_back(buffer_[rpos_]);
            buffer_[rpos_] = 0;
            rpos_ = (rpos_ + 1) & mask_;
        }
        size_ = 0;
    }

private:
    std::vector<float> buffer_;
    size_t ola_add_pos_{};
    size_t ola_add_end_{};
    size_t rpos_{};
    size_t size_{};
    size_t mask_{};
};

} // namespace detail

// ------------------------------------------------------------
// RunWsola
// ------------------------------------------------------------
/**
 * @brief WSOLA 时间拉伸（Waveform Similarity Overlap-Add）。
 *
 * 通过相关性匹配选择最相似的重叠位置，实现无相位失真的时间拉伸。
 *
 * @param input         输入单声道信号
 * @param stretch_ratio 时间拉伸比（>1 变慢变长，<1 变快变短）
 * @return 拉伸后的输出信号
 */
static inline std::vector<float> RunWsola(std::span<const float> input, float stretch_ratio = 2.0f) {
    constexpr size_t synthsis_block = 1024;
    constexpr size_t overlap_length = 256;
    constexpr size_t synthsis_hop = synthsis_block - overlap_length;
    const size_t analyze_hop = static_cast<size_t>(static_cast<float>(synthsis_hop) / stretch_ratio);
    constexpr int search_range = overlap_length * 2;

    std::array<float, synthsis_block> hann_window;
    hann_window.fill(1.0f);
    for (size_t i = 0; i < overlap_length; ++i) {
        float x = static_cast<float>(i) / static_cast<float>(overlap_length);
        hann_window[i] = std::sin(x * std::numbers::pi_v<float> * 0.5f);
        hann_window[i + (synthsis_block - overlap_length)] = std::cos(x * std::numbers::pi_v<float> * 0.5f);
    }

    std::vector<float> output;
    std::vector<float> copy(input.begin(), input.end()); // MonoReader 需要可变 span
    qwqdsp_segement::MonoReader slice{copy};
    detail::OlaBuffer ola;
    ola.Init(synthsis_block * 4);

    // init
    std::array<float, synthsis_block> fill_block;
    slice.Read(synthsis_block, fill_block);
    slice.Step(analyze_hop);
    ola.Add(fill_block.data(), hann_window.data(), synthsis_block, synthsis_hop);
    ola.ReadPushback(output);

    while (!slice.IsEnd()) {
        // step
        size_t this_rpos = slice.GetReadpos();

        // find match block
        constexpr size_t search_block_size = search_range + overlap_length;
        std::array<float, search_block_size> search_block;
        int search_begin = this_rpos - search_range / 2;
        if (search_begin < 0)
            search_begin = 0;
        slice.ReadAbsolute(search_begin, search_block_size, search_block.data());
        std::array<float, overlap_length> match_block;
        std::copy_n(fill_block.data() + (synthsis_block - overlap_length), overlap_length, match_block.data());

        // correlation
        float max_corr = 99999999999.0f;
        int best_offset = 0;
        for (size_t i = 0; i < search_range; ++i) {
            float sum = 0;
            for (size_t j = 0; j < overlap_length; ++j) {
                float a = match_block[j];
                float b = search_block[j + i];
                sum += (a - b) * (a - b);
            }
            if (sum < max_corr) {
                max_corr = sum;
                best_offset = i;
            }
        }

        // get
        slice.ReadAbsolute(search_begin + best_offset, synthsis_block, fill_block.data());
        ola.Add(fill_block.data(), hann_window.data(), synthsis_block, synthsis_hop);
        ola.ReadPushback(output);
        slice.Step(analyze_hop);
    }

    return output;
}

} // namespace qwqdsp_test
