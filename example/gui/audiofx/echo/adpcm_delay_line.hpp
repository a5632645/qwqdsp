#pragma once

// ------------------------------------------------------------
//  AdpcmDelayLine — Q8.23 域的 IMA ADPCM 块压缩延迟线
//
//  每 64 个样本压成一块，块内布局共 kBlockBytes 字节：
//    偏移 0..3         int32（Q8.23，小端）= 块开始前的预测值
//    偏移 4            int8               = 块开始前的步长索引
//    偏移 5..kBlockBytes-1                 = 64 个码字，字节内低位存放靠前的样本
//
//  码字宽度是模板参数：kMagnitudeBits 个幅度位 + 1 个符号位。幅度位为 3 时就是标准
//  4bit IMA（每字节两个码字，低半字节在先）；为 7 时是 8bit 码字（每样本一字节）。
//  重建步距与 IMA 一致：delta = (2·m + 1)·step / 2^kMagnitudeBits。
// ------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "quantized_delay_line.hpp"

namespace echo_demo {

/// @brief 编解码器状态（Q8.23 预测值 + 步长索引）
struct AdpcmState {
    int32_t predictor = 0;
    int32_t step_index = 0;
};

// ------------------------------------------------------------
//  AdpcmDelayLine — 直接以压缩块为存储单元的整数延迟线
//
//  写侧按 64 个样本聚成一块编码，读侧按需解出目标所在块。
// ------------------------------------------------------------
template <int kMagnitudeBits = 7>
class AdpcmDelayLine {
    static_assert(kMagnitudeBits >= 3 && kMagnitudeBits <= 7, "幅度位数需落在 [3, 7]");
public:
    using Value = int32_t;

    static constexpr int kCodeBits = kMagnitudeBits + 1;
    static constexpr size_t kSamplesPerBlock = 64;
    static constexpr size_t kCodesPerByte = 8 / kCodeBits;
    static constexpr size_t kPredictorOffset = 0;
    static constexpr size_t kIndexOffset = 4;
    static constexpr size_t kDataOffset = 5;
    static constexpr size_t kDataBytes = kSamplesPerBlock / kCodesPerByte;
    static constexpr size_t kBlockBytes = kDataOffset + kDataBytes;
    static constexpr int32_t kStepIndexMax = 88;
    static constexpr uint8_t kMagnitudeMask = static_cast<uint8_t>((1u << kMagnitudeBits) - 1);
    static constexpr uint8_t kSignMask = static_cast<uint8_t>(1u << kMagnitudeBits);
    static constexpr uint8_t kCodeMask = static_cast<uint8_t>(kMagnitudeMask | kSignMask);

    /// @brief 最小可用延迟样本数：整块写完才能解码
    static constexpr size_t kMinDelaySamples = kSamplesPerBlock;

    /// @brief 把一块样本编码进 kBlockBytes 字节，并把块开始前的状态写进块头
    /// @param[in,out] state 传入块开始前的状态，传出块结束后的状态
    /// @param[in] samples 块内 kSamplesPerBlock 个 Q8.23 样本
    /// @param[out] block 收 kBlockBytes 字节的块
    static void EncodeBlock(AdpcmState& state, const int32_t* samples, uint8_t* block) noexcept {
        StorePredictor(block, state.predictor);
        block[kIndexOffset] = static_cast<uint8_t>(static_cast<int8_t>(state.step_index));
        for (size_t i = 0; i < kSamplesPerBlock; ++i) {
            uint8_t const code = Encode(state, samples[i]);
            size_t const shift = (i % kCodesPerByte) * kCodeBits;
            uint8_t& byte = block[kDataOffset + i / kCodesPerByte];
            if (shift == 0) {
                byte = 0;
            }
            byte = static_cast<uint8_t>(byte | (code << shift));
        }
    }

    /// @brief 从 kBlockBytes 字节块解出 kSamplesPerBlock 个 Q8.23 样本，起始状态取自块头
    static void DecodeBlock(const uint8_t* block, int32_t* samples) noexcept {
        // 钳到预测值上限：块头若被写坏也不会让差分运算溢出
        AdpcmState state{std::clamp(LoadPredictor(block), -kPredictorLimit, kPredictorLimit),
                         static_cast<int32_t>(static_cast<int8_t>(block[kIndexOffset]))};
        for (size_t i = 0; i < kSamplesPerBlock; ++i) {
            size_t const shift = (i % kCodesPerByte) * kCodeBits;
            uint8_t const code =
                static_cast<uint8_t>((block[kDataOffset + i / kCodesPerByte] >> shift) & kCodeMask);
            samples[i] = Advance(state, code);
        }
    }

    /// @brief 按容量分配块环，容量向上取到 2 的幂
    void Init(size_t max_samples) {
        size_t blocks = (max_samples + kSamplesPerBlock - 1) / kSamplesPerBlock;
        size_t size = 2;
        while (size < blocks) {
            size <<= 1;
        }
        block_count_ = size;
        block_mask_ = size - 1;
        buffer_.assign(block_count_ * kBlockBytes, 0);
        pending_.fill(0);
        scratch_.fill(0);
        encoder_state_ = AdpcmState{};
        write_index_ = 0;
        cached_block_ = kNoBlock;
    }

    /// @brief 编码并写入一个样本，凑满一块时整块落盘
    void Push(Value x) noexcept {
        size_t const in_block = write_index_ % kSamplesPerBlock;
        pending_[in_block] = x;
        ++write_index_;
        if (in_block + 1 == kSamplesPerBlock) {
            EncodeBlock(encoder_state_, pending_.data(), &buffer_[BlockOffset(write_index_ / kSamplesPerBlock - 1)]);
        }
    }

    /// @brief 读取 delay_samples 次 Push 之前写入的样本
    /// @param[in] delay_samples 延迟样本数，内部钳到 [kMinDelaySamples, Capacity()]
    Value GetBeforePush(size_t delay_samples) noexcept {
        size_t const delay = std::clamp(delay_samples, kMinDelaySamples, Capacity());
        if (write_index_ < delay) {
            return 0;
        }
        size_t const target = write_index_ - delay;
        size_t const block = target / kSamplesPerBlock;
        if (block != cached_block_) {
            DecodeBlock(&buffer_[BlockOffset(block)], scratch_.data());
            cached_block_ = block;
        }
        return scratch_[target % kSamplesPerBlock];
    }

    /// @brief 可用延迟样本数上限（块对齐，留一块给正在写入的块）
    size_t Capacity() const noexcept { return (block_count_ - 1) * kSamplesPerBlock; }

    /// @brief 分配出的块数
    size_t BlockCount() const noexcept { return block_count_; }
private:
    static constexpr size_t kNoBlock = static_cast<size_t>(-1);
    /// @brief 预测值钳位上限（Q8.23 的整数位范围 ±128），保证差分不出 int32
    static constexpr int32_t kPredictorLimit = 1 << 30;

    static constexpr int8_t kIndexAdjust[8] = {-1, -1, -1, -1, 2, 4, 6, 8};
    static constexpr int16_t kStepTable[89] = {
        7,     8,     9,     10,    11,    12,    13,    14,    16,    17,    19,    21,    23,    25,    28,
        31,    34,    37,    41,    45,    50,    55,    60,    66,    73,    80,    88,    97,    107,   118,
        130,   143,   157,   173,   190,   209,   230,   253,   279,   307,   337,   371,   408,   449,   494,
        544,   598,   658,   724,   796,   876,   963,   1060,  1166,  1282,  1411,  1552,  1707,  1878,  2066,
        2272,  2499,  2749,  3024,  3327,  3660,  4026,  4428,  4871,  5358,  5894,  6484,  7132,  7845,  8630,
        9493,  10442, 11487, 12635, 13899, 15289, 16818, 18500, 20350, 22385, 24623, 27086, 29794, 32767};

    /// @brief 第 index 档步长的 Q8.23 值
    static int32_t Step(int32_t index) noexcept { return static_cast<int32_t>(kStepTable[index]) << 8; }

    /// @brief 按码字推进状态并返回重建样本的 Q8.23 值（编码与解码共用，保证两端状态一致）
    static int32_t Advance(AdpcmState& state, uint8_t code) noexcept {
        int32_t const step = Step(state.step_index);
        // 幅度 m 的重建量 = (2m + 1)·step / 2^kMagnitudeBits：先加上半个最小步距
        int32_t delta = step >> kMagnitudeBits;
        for (int k = 0; k < kMagnitudeBits; ++k) {
            if ((code & (1u << (kMagnitudeBits - 1 - k))) != 0) {
                delta += step >> k;
            }
        }
        state.predictor = ((code & kSignMask) != 0) ? state.predictor - delta : state.predictor + delta;
        state.predictor = std::clamp(state.predictor, -kPredictorLimit, kPredictorLimit);
        // 步长索引只看幅度的粗档位（8 档），与 4bit 版保持同一套自适应动态
        int32_t const magnitude = code & kMagnitudeMask;
        state.step_index = std::clamp(state.step_index + kIndexAdjust[magnitude >> (kMagnitudeBits - 3)], 0,
                                      kStepIndexMax);
        return state.predictor;
    }

    /// @brief 把一个样本编码成码字并推进状态
    static uint8_t Encode(AdpcmState& state, int32_t sample) noexcept {
        int32_t const step = Step(state.step_index);
        int64_t diff = static_cast<int64_t>(sample) - state.predictor;
        uint8_t code = 0;
        if (diff < 0) {
            code = kSignMask;
            diff = -diff;
        }
        // 逐位与阈值比较等价于取最近的幅度档：阈值依次为 step/2^k
        for (int k = 0; k < kMagnitudeBits; ++k) {
            int32_t const threshold = step >> k;
            if (diff >= threshold) {
                code = static_cast<uint8_t>(code | (1u << (kMagnitudeBits - 1 - k)));
                diff -= threshold;
            }
        }
        Advance(state, code);
        return code;
    }

    static void StorePredictor(uint8_t* block, int32_t value) noexcept {
        uint32_t const u = static_cast<uint32_t>(value);
        block[0] = static_cast<uint8_t>(u);
        block[1] = static_cast<uint8_t>(u >> 8);
        block[2] = static_cast<uint8_t>(u >> 16);
        block[3] = static_cast<uint8_t>(u >> 24);
    }

    static int32_t LoadPredictor(const uint8_t* block) noexcept {
        return static_cast<int32_t>(static_cast<uint32_t>(block[0]) | (static_cast<uint32_t>(block[1]) << 8)
                                    | (static_cast<uint32_t>(block[2]) << 16)
                                    | (static_cast<uint32_t>(block[3]) << 24));
    }

    size_t BlockOffset(size_t block) const noexcept { return (block & block_mask_) * kBlockBytes; }

    std::vector<uint8_t> buffer_;
    std::array<int32_t, kSamplesPerBlock> pending_{};
    std::array<int32_t, kSamplesPerBlock> scratch_{};
    AdpcmState encoder_state_;
    size_t write_index_ = 0;
    size_t block_count_ = 0;
    size_t block_mask_ = 0;
    size_t cached_block_ = kNoBlock;
};

// ------------------------------------------------------------
//  Adpcm8Storage — 8bit IMA 块压缩存储，运算域沿用 Q8.23（int32）
// ------------------------------------------------------------
struct Adpcm8Storage : Q8_23Storage {
    using Line = AdpcmDelayLine<7>;

    static constexpr const char* kName = "IMA23/8";
    static constexpr const char* kDescription = "8 bit IMA storage + Q8.23 (int32) arithmetic";
};

} // namespace echo_demo
