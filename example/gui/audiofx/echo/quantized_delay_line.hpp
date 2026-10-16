#pragma once

// ------------------------------------------------------------
//  整数延迟线（无插值）与四种样本存储格式
// ------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace echo_demo {

/// @brief 24 位有符号定点数，按小端 3 字节紧凑存放
struct Int24 {
    uint8_t bytes[3];
};

/// @brief 把浮点样本按 2^-kFracBits 定点量化，向最近取整并饱和到格式表示范围
/// @param[in] x 待量化样本，超出表示范围时饱和
/// @return kTotalBits 位有符号整数形式，kTotalBits 不得超过 32
template <int kFracBits, int kTotalBits>
static int32_t Quantize(float x) noexcept {
    constexpr double kScale = static_cast<double>(uint64_t{1} << kFracBits);
    constexpr double kMin = -static_cast<double>(int64_t{1} << (kTotalBits - 1));
    constexpr double kMax = static_cast<double>((int64_t{1} << (kTotalBits - 1)) - 1);
    // 在 double 域内先饱和再取整，避免 llrint 收到溢出参数
    return static_cast<int32_t>(std::llrint(std::clamp(static_cast<double>(x) * kScale, kMin, kMax)));
}

// ------------------------------------------------------------
//  存储策略 — 每种格式给出量化、反量化、延迟线类型与界面显示用的名称
// ------------------------------------------------------------

template <typename Storage>
class QuantizedDelayLine;

/// @brief 单精度浮点存储，不做量化
struct FloatStorage {
    using Value = float;
    using Line = QuantizedDelayLine<FloatStorage>;
    static constexpr const char* kName = "f32";
    static constexpr const char* kDescription = "32 bit float, no quantization";
    static constexpr int kFracBits = 0;
    static constexpr int kTotalBits = 32;

    static Value FromFloat(float x) noexcept { return x; }
    static float ToFloat(Value v) noexcept { return v; }
};

/// @brief Q8.23 定点存储：23 位小数、上限 ±128
struct Q8_23Storage {
    using Value = int32_t;
    using Line = QuantizedDelayLine<Q8_23Storage>;
    static constexpr const char* kName = "Q8.23";
    static constexpr const char* kDescription = "23 frac bits, 32 bit/sample, range +-128";
    static constexpr int kFracBits = 23;
    static constexpr int kTotalBits = 32;

    static Value FromFloat(float x) noexcept { return Quantize<kFracBits, kTotalBits>(x); }
    static float ToFloat(Value v) noexcept { return static_cast<float>(v) * kInvScale; }
private:
    static constexpr float kInvScale = 1.0f / static_cast<float>(uint32_t{1} << kFracBits);
};

/// @brief Q3.20 定点存储：20 位小数、上限 ±8，每样本占 3 字节
struct Q3_20Storage {
    using Value = Int24;
    using Line = QuantizedDelayLine<Q3_20Storage>;
    static constexpr const char* kName = "Q3.20";
    static constexpr const char* kDescription = "20 frac bits, 24 bit/sample, range +-8";
    static constexpr int kFracBits = 20;
    static constexpr int kTotalBits = 24;

    static Value FromFloat(float x) noexcept {
        int32_t const v = Quantize<kFracBits, kTotalBits>(x);
        return Value{{static_cast<uint8_t>(v), static_cast<uint8_t>(v >> 8), static_cast<uint8_t>(v >> 16)}};
    }

    static float ToFloat(Value v) noexcept {
        int32_t i = static_cast<int32_t>(v.bytes[0]) | (static_cast<int32_t>(v.bytes[1]) << 8)
                  | (static_cast<int32_t>(v.bytes[2]) << 16);
        if ((i & 0x800000) != 0) {
            i -= 0x1000000;
        }
        return static_cast<float>(i) * kInvScale;
    }
private:
    static constexpr float kInvScale = 1.0f / static_cast<float>(uint32_t{1} << kFracBits);
};

/// @brief Q3.12 定点存储：12 位小数、上限 ±8
struct Q3_12Storage {
    using Value = int16_t;
    using Line = QuantizedDelayLine<Q3_12Storage>;
    static constexpr const char* kName = "Q3.12";
    static constexpr const char* kDescription = "12 frac bits, 16 bit/sample, range +-8";
    static constexpr int kFracBits = 12;
    static constexpr int kTotalBits = 16;

    static Value FromFloat(float x) noexcept {
        return static_cast<Value>(Quantize<kFracBits, kTotalBits>(x));
    }
    static float ToFloat(Value v) noexcept { return static_cast<float>(v) * kInvScale; }
private:
    static constexpr float kInvScale = 1.0f / static_cast<float>(uint32_t{1} << kFracBits);
};

// ------------------------------------------------------------
//  QuantizedDelayLine — 整数延迟、无插值的单声道延迟线
//
//  容量向上取到 2 的幂，读写都经过 Storage 的量化与反量化。
// ------------------------------------------------------------
template <typename Storage>
class QuantizedDelayLine {
public:
    /// @brief 最小可用延迟样本数：每个样本独立存放，可读到刚写入的前一个
    static constexpr size_t kMinDelaySamples = 1;

    /// @brief 按容量分配缓冲
    /// @param[in] max_samples 可用延迟样本数上限，向上取到 2 的幂
    void Init(size_t max_samples) {
        size_ = 1;
        while (size_ < max_samples) {
            size_ <<= 1;
        }
        mask_ = size_ - 1;
        wpos_ = 0;
        buffer_.assign(size_, Storage::FromFloat(0.0f));
    }

    /// @brief 写入一个样本
    void Push(float x) noexcept {
        buffer_[wpos_] = Storage::FromFloat(x);
        wpos_ = (wpos_ + 1) & mask_;
    }

    /// @brief 读取 delay_samples 次 Push 之前写入的样本
    /// @param[in] delay_samples 延迟样本数，取值须落在 [1, Capacity()]
    /// @return 该位置存储值的反量化结果
    float GetBeforePush(size_t delay_samples) const noexcept {
        return Storage::ToFloat(buffer_[(wpos_ - delay_samples) & mask_]);
    }

    /// @brief 可用延迟样本数上限
    size_t Capacity() const noexcept { return size_; }
private:
    std::vector<typename Storage::Value> buffer_;
    size_t wpos_ = 0;
    size_t mask_ = 0;
    size_t size_ = 0;
};

} // namespace echo_demo
