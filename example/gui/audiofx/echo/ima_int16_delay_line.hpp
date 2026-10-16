#pragma once

// ------------------------------------------------------------
//  ImaInt16DelayLine — int16 域的 4bit IMA ADPCM 块压缩延迟线
//
//  每 32 个样点压成一块，块内布局共 19 字节：
//    偏移 0..1  int16（小端）= 块开始前的预测值
//    偏移 2     int8         = 块开始前的步长索引（0..88）
//    偏移 3..18 16 字节       = 32 个 4bit 码字，低半字节存放块内靠前的样本
//
//  块头只承载跨块的编解码状态，块内 32 个样点全部被编码；编解码都在 int16 域内
//  做饱和加，步长表用标准 IMA/DVI4 原值。与标准 IMA WAV 的差别只有 nibble 顺序。
// ------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace echo_demo {

/// @brief int16 域 IMA ADPCM 编解码器状态
struct ImaInt16State {
    int32_t predictor = 0;
    int32_t index = 0;
};

// ------------------------------------------------------------
//  ImaInt16DelayLine — 以 19 字节块为存储单元的整数延迟线
//
//  写侧按 32 个样点聚成一块编码，读侧按需解出目标所在块。
// ------------------------------------------------------------
class ImaInt16DelayLine {
public:
    static constexpr size_t kSamplesPerBlock = 32;
    static constexpr size_t kBlockBytes = 19;
    static constexpr size_t kIndexOffset = 2;
    static constexpr size_t kDataOffset = 3;
    static constexpr size_t kStepIndexCount = 89;
    static constexpr int32_t kStepIndexMax = static_cast<int32_t>(kStepIndexCount) - 1;

    /// @brief 最小可用延迟样本数：整块写完才能解码
    static constexpr size_t kMinDelaySamples = kSamplesPerBlock;

    /// @brief 把一块样点编码进 19 字节，并把块开始前的状态写进块头
    /// @param[in,out] state 传入块开始前的状态，传出块结束后的状态
    /// @param[in] samples 块内 kSamplesPerBlock 个 int16 样点
    /// @param[out] block 收 kBlockBytes 字节的块
    static void EncodeBlock(ImaInt16State& state, const int16_t* samples, uint8_t* block) noexcept {
        int32_t predictor = state.predictor;
        int32_t index = state.index;
        StoreLe16(block, static_cast<int16_t>(predictor));
        block[kIndexOffset] = static_cast<uint8_t>(static_cast<int8_t>(index));

        uint8_t* p = block + kDataOffset;
        for (size_t i = 0; i < kSamplesPerBlock / 2; ++i) {
            uint8_t byte = 0;
            for (size_t k = 0; k < 2; ++k) {
                int32_t const sample = samples[2 * i + k];
                int32_t step = kStepTable[index];
                int32_t diff = sample - predictor;
                int32_t code = 0;
                int32_t delta = step >> 3;
                if (diff < 0) {
                    diff = -diff;
                    code = 8;
                }
                if (diff >= step) {
                    code |= 4;
                    diff -= step;
                    delta += step;
                }
                step >>= 1;
                if (diff >= step) {
                    code |= 2;
                    diff -= step;
                    delta += step;
                }
                step >>= 1;
                if (diff >= step) {
                    code |= 1;
                    delta += step;
                }
                predictor = Sat16(predictor + (((code & 8) != 0) ? -delta : delta));
                index = std::clamp(index + kIndexTable[code & 7], 0, kStepIndexMax);
                byte = static_cast<uint8_t>(byte | ((k == 0) ? code : (code << 4)));
            }
            *p++ = byte;
        }

        state.predictor = predictor;
        state.index = index;
    }

    /// @brief 从 19 字节块解出 kSamplesPerBlock 个 int16 样点，起始状态取自块头
    static void DecodeBlock(const uint8_t* block, int16_t* samples) noexcept {
        int32_t predictor = LoadLe16(block);
        int32_t index = static_cast<int8_t>(block[kIndexOffset]);
        const uint8_t* p = block + kDataOffset;
        for (size_t i = 0; i < kSamplesPerBlock / 2; ++i) {
            uint8_t const byte = *p++;
            int32_t code = byte & 0x0F;
            for (size_t k = 0; k < 2; ++k) {
                int32_t const step = kStepTable[index];
                int32_t delta = step >> 3;
                if ((code & 4) != 0) {
                    delta += step;
                }
                if ((code & 2) != 0) {
                    delta += step >> 1;
                }
                if ((code & 1) != 0) {
                    delta += step >> 2;
                }
                if ((code & 8) != 0) {
                    delta = -delta;
                }
                predictor = Sat16(predictor + delta);
                samples[2 * i + k] = static_cast<int16_t>(predictor);
                index = std::clamp(index + kIndexTable[code & 7], 0, kStepIndexMax);
                code = byte >> 4;
            }
        }
    }

    /// @brief 按容量分配块环，块数按 max_samples 精确算出（不向上取 2 的幂）
    void Init(size_t max_samples) {
        block_count_ = std::max<size_t>((max_samples + kSamplesPerBlock - 1) / kSamplesPerBlock, 2);
        buffer_.assign(block_count_ * kBlockBytes, 0);
        pending_.fill(0);
        scratch_.fill(0);
        state_ = ImaInt16State{};
        write_index_ = 0;
        cached_block_ = kNoBlock;
    }

    /// @brief 编码并写入一个样本，凑满一块时整块落盘
    void Push(float x) noexcept {
        size_t const in_block = write_index_ % kSamplesPerBlock;
        pending_[in_block] = ToInt16(x);
        ++write_index_;
        if (in_block + 1 == kSamplesPerBlock) {
            EncodeBlock(state_, pending_.data(), &buffer_[BlockOffset(write_index_ / kSamplesPerBlock - 1)]);
        }
    }

    /// @brief 读取 delay_samples 次 Push 之前写入的样本
    /// @param[in] delay_samples 延迟样本数，内部钳到 [kMinDelaySamples, Capacity()]
    float GetBeforePush(size_t delay_samples) noexcept {
        size_t const delay = std::clamp(delay_samples, kMinDelaySamples, Capacity());
        if (write_index_ < delay) {
            return 0.0f;
        }
        size_t const target = write_index_ - delay;
        size_t const block = target / kSamplesPerBlock;
        if (block != cached_block_) {
            DecodeBlock(&buffer_[BlockOffset(block)], scratch_.data());
            cached_block_ = block;
        }
        return ToFloat(scratch_[target % kSamplesPerBlock]);
    }

    /// @brief 可用延迟样本数上限（块对齐，留一块给正在写入的块）
    size_t Capacity() const noexcept { return (block_count_ - 1) * kSamplesPerBlock; }

    /// @brief 分配出的块数
    size_t BlockCount() const noexcept { return block_count_; }
private:
    static constexpr size_t kNoBlock = static_cast<size_t>(-1);
    static constexpr float kScale = 32768.0f;
    static constexpr float kInvScale = 1.0f / kScale;

    static constexpr int8_t kIndexTable[8] = {-1, -1, -1, -1, 2, 4, 6, 8};
    static constexpr int16_t kStepTable[89] = {
        7,     8,     9,     10,    11,    12,    13,    14,    16,    17,    19,    21,    23,    25,    28,
        31,    34,    37,    41,    45,    50,    55,    60,    66,    73,    80,    88,    97,    107,   118,
        130,   143,   157,   173,   190,   209,   230,   253,   279,   307,   337,   371,   408,   449,   494,
        544,   598,   658,   724,   796,   876,   963,   1060,  1166,  1282,  1411,  1552,  1707,  1878,  2066,
        2272,  2499,  2749,  3024,  3327,  3660,  4026,  4428,  4871,  5358,  5894,  6484,  7132,  7845,  8630,
        9493,  10442, 11487, 12635, 13899, 15289, 16818, 18500, 20350, 22385, 24623, 27086, 29794, 32767};

    static int32_t Sat16(int32_t v) noexcept { return std::clamp(v, -32768, 32767); }

    static int16_t ToInt16(float x) noexcept {
        return static_cast<int16_t>(std::lrintf(std::clamp(x * kScale, -kScale, kScale - 1.0f)));
    }

    static float ToFloat(int16_t v) noexcept { return static_cast<float>(v) * kInvScale; }

    static void StoreLe16(uint8_t* p, int16_t v) noexcept {
        uint16_t const u = static_cast<uint16_t>(v);
        p[0] = static_cast<uint8_t>(u);
        p[1] = static_cast<uint8_t>(u >> 8);
    }

    static int16_t LoadLe16(const uint8_t* p) noexcept {
        return static_cast<int16_t>(static_cast<uint16_t>(p[0]) | (static_cast<uint16_t>(p[1]) << 8));
    }

    size_t BlockOffset(size_t block) const noexcept { return (block % block_count_) * kBlockBytes; }

    std::vector<uint8_t> buffer_;
    std::array<int16_t, kSamplesPerBlock> pending_{};
    std::array<int16_t, kSamplesPerBlock> scratch_{};
    ImaInt16State state_;
    size_t write_index_ = 0;
    size_t block_count_ = 0;
    size_t cached_block_ = kNoBlock;
};

// ------------------------------------------------------------
//  ImaInt16Storage — int16 域的 4bit IMA ADPCM 块压缩格式
// ------------------------------------------------------------
struct ImaInt16Storage {
    using Line = ImaInt16DelayLine;

    static constexpr const char* kName = "IMA16";
    static constexpr const char* kDescription = "4 bit IMA, int16 domain, 32 samples/19 B";
};

} // namespace echo_demo
