#pragma once

// ------------------------------------------------------------
//  延迟线的标量运算域与整数延迟线
//
//  每个存储策略同时定义该格式的**运算域**：Value 类型、浮点换算、饱和加、四舍五入
//  定点乘；整个 echo（读取、反馈乘、混音、写回）都在这套运算里完成。16 位格式额外
//  提供一阶误差反馈量化，把量化余量高通整形。
// ------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace echo_demo {

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

/// @brief 24 位有符号定点数，按小端 3 字节紧凑存放
struct Int24 {
    uint8_t bytes[3];

    /// @brief 取出符号扩展后的整数值
    int32_t ToInt() const noexcept {
        int32_t v = static_cast<int32_t>(bytes[0]) | (static_cast<int32_t>(bytes[1]) << 8)
                  | (static_cast<int32_t>(bytes[2]) << 16);
        if ((v & 0x800000) != 0) {
            v -= 0x1000000;
        }
        return v;
    }

    /// @brief 由整数值构造，超出 24 位范围时饱和
    static Int24 FromInt(int32_t value) noexcept {
        int32_t const v = std::clamp(value, -(1 << 23), (1 << 23) - 1);
        uint32_t const u = static_cast<uint32_t>(v);
        return Int24{{static_cast<uint8_t>(u), static_cast<uint8_t>(u >> 8), static_cast<uint8_t>(u >> 16)}};
    }
};

namespace detail {

/// @brief 一阶误差反馈量化：把细网格输入（1/256 LSB 单位）量化到 16 位定点网格
/// @param[in] v_fine 细网格输入值，不含误差项
/// @param[in,out] error 传入上一轮余量，传出本轮余量；两者都保持在半个 LSB 以内
/// @return 量化后的 16 位定点值
inline int16_t QuantizeInt16Fine(int64_t v_fine, int64_t& error) noexcept {
    constexpr int64_t kMinFine = -(int64_t{1} << 15) * 256;
    constexpr int64_t kMaxFine = ((int64_t{1} << 15) - 1) * 256;
    v_fine = std::clamp(v_fine + error, kMinFine, kMaxFine);
    int64_t const q = ((v_fine + 128) >> 8) << 8;
    error = v_fine - q;
    return static_cast<int16_t>(q >> 8);
}

} // namespace detail

template <typename Storage>
class QuantizedDelayLine;

// ------------------------------------------------------------
//  FloatStorage — 纯浮点运算域，不做任何量化
// ------------------------------------------------------------
struct FloatStorage {
    using Value = float;
    using Line = QuantizedDelayLine<FloatStorage>;

    static constexpr const char* kName = "f32";
    static constexpr const char* kDescription = "f32 arithmetic, no quantization";
    static constexpr bool kUseErrorFeedback = false;

    static Value FromFloat(float x) noexcept { return x; }
    static float ToFloat(Value v) noexcept { return v; }
    static Value Add(Value a, Value b) noexcept { return a + b; }
    static Value Mul(Value a, Value b) noexcept { return a * b; }
};

// ------------------------------------------------------------
//  Q8_23Storage — int32 上的 Q8.23 运算域
// ------------------------------------------------------------
struct Q8_23Storage {
    using Value = int32_t;
    using Line = QuantizedDelayLine<Q8_23Storage>;

    static constexpr const char* kName = "Q8.23";
    static constexpr const char* kDescription = "Q8.23 arithmetic (int32), all in-format";
    static constexpr int kFracBits = 23;
    static constexpr bool kUseErrorFeedback = false;

    static Value FromFloat(float x) noexcept { return Quantize<kFracBits, 32>(x); }
    static float ToFloat(Value v) noexcept { return static_cast<float>(v) * kInvScale; }

    /// @brief 饱和加
    static Value Add(Value a, Value b) noexcept {
        return static_cast<Value>(std::clamp<int64_t>(static_cast<int64_t>(a) + b, kMin, kMax));
    }

    /// @brief 四舍五入的定点乘，结果饱和
    static Value Mul(Value a, Value b) noexcept {
        int64_t const product = static_cast<int64_t>(a) * b;
        return static_cast<Value>(std::clamp<int64_t>((product + kHalf) >> kFracBits, kMin, kMax));
    }
private:
    static constexpr int64_t kMin = -(int64_t{1} << 31);
    static constexpr int64_t kMax = (int64_t{1} << 31) - 1;
    static constexpr int64_t kHalf = int64_t{1} << (kFracBits - 1);
    static constexpr float kInvScale = 1.0f / static_cast<float>(uint32_t{1} << kFracBits);
};

// ------------------------------------------------------------
//  Q3_20Storage — int24 上的 Q3.20 运算域
// ------------------------------------------------------------
struct Q3_20Storage {
    using Value = Int24;
    using Line = QuantizedDelayLine<Q3_20Storage>;

    static constexpr const char* kName = "Q3.20";
    static constexpr const char* kDescription = "Q3.20 arithmetic (int24), all in-format";
    static constexpr int kFracBits = 20;
    static constexpr bool kUseErrorFeedback = false;

    static Value FromFloat(float x) noexcept { return Value::FromInt(Quantize<kFracBits, 24>(x)); }
    static float ToFloat(Value v) noexcept { return static_cast<float>(v.ToInt()) * kInvScale; }

    /// @brief 饱和加
    static Value Add(Value a, Value b) noexcept {
        return Value::FromInt(static_cast<int32_t>(std::clamp<int64_t>(
            static_cast<int64_t>(a.ToInt()) + b.ToInt(), kMin, kMax)));
    }

    /// @brief 四舍五入的定点乘，结果饱和
    static Value Mul(Value a, Value b) noexcept {
        int64_t const product = static_cast<int64_t>(a.ToInt()) * b.ToInt();
        return Value::FromInt(static_cast<int32_t>(
            std::clamp<int64_t>((product + kHalf) >> kFracBits, kMin, kMax)));
    }
private:
    static constexpr int64_t kMin = -(int64_t{1} << 23);
    static constexpr int64_t kMax = (int64_t{1} << 23) - 1;
    static constexpr int64_t kHalf = int64_t{1} << (kFracBits - 1);
    static constexpr float kInvScale = 1.0f / static_cast<float>(uint32_t{1} << kFracBits);
};

// ------------------------------------------------------------
//  Q16Arithmetic — int16 上的定点运算域（含一阶误差反馈）
//
//  量化余量保存在 1/256 LSB 的细网格上，所以 12 位小数与 15 位小数两档共用同一套实现。
// ------------------------------------------------------------
template <int kFrac>
struct Q16Arithmetic {
    using Value = int16_t;

    static constexpr bool kUseErrorFeedback = true;
    static constexpr int kFracBits = kFrac;

    static Value FromFloat(float x) noexcept { return static_cast<Value>(Quantize<kFrac, 16>(x)); }
    static float ToFloat(Value v) noexcept { return static_cast<float>(v) * kInvScale; }

    /// @brief 饱和加
    static Value Add(Value a, Value b) noexcept {
        return static_cast<Value>(std::clamp<int32_t>(static_cast<int32_t>(a) + b, -32768, 32767));
    }

    /// @brief 四舍五入的定点乘，结果饱和
    static Value Mul(Value a, Value b) noexcept {
        return static_cast<Value>(RoundShift(static_cast<int64_t>(a) * b));
    }

    /// @brief 带一阶误差反馈的输入量化
    static Value FromFloatShaped(float x, int64_t& error) noexcept {
        double const scale = std::exp2(static_cast<double>(kFrac) + 8.0);
        return detail::QuantizeInt16Fine(std::llrint(static_cast<double>(x) * scale), error);
    }

    /// @brief 带一阶误差反馈的定点乘
    static Value MulShaped(Value a, Value b, int64_t& error) noexcept {
        int64_t const product = static_cast<int64_t>(a) * b;
        int64_t const fine = (product + (int64_t{1} << (kFrac - 9))) >> (kFrac - 8);
        return detail::QuantizeInt16Fine(fine, error);
    }
private:
    static int64_t RoundShift(int64_t product) noexcept {
        return std::clamp<int64_t>((product + (int64_t{1} << (kFrac - 1))) >> kFrac, -32768, 32767);
    }

    static constexpr float kInvScale = 1.0f / static_cast<float>(uint32_t{1} << kFrac);
};

// ------------------------------------------------------------
//  Q24Arithmetic — 24bit PCM 混音域（值右对齐放在 int32，1.0 = 2^23）
//
//  对应厂商链里「混音与反馈仍走 24bit」的那一段；饱和到 24bit（±0x7FFFFF）。
// ------------------------------------------------------------
struct Q24Arithmetic {
    using Value = int32_t;

    static constexpr bool kUseErrorFeedback = false;
    static constexpr int kFracBits = 23;

    static Value FromFloat(float x) noexcept { return Quantize<kFracBits, 24>(x); }
    static float ToFloat(Value v) noexcept { return static_cast<float>(v) * kInvScale; }

    /// @brief 饱和加
    static Value Add(Value a, Value b) noexcept {
        return static_cast<Value>(std::clamp<int64_t>(static_cast<int64_t>(a) + b, kMin, kMax));
    }

    /// @brief 四舍五入的定点乘，结果饱和
    static Value Mul(Value a, Value b) noexcept {
        int64_t const product = static_cast<int64_t>(a) * b;
        return static_cast<Value>(std::clamp<int64_t>((product + kHalf) >> kFracBits, kMin, kMax));
    }
private:
    static constexpr int64_t kMin = -(int64_t{1} << 23);
    static constexpr int64_t kMax = (int64_t{1} << 23) - 1;
    static constexpr int64_t kHalf = int64_t{1} << (kFracBits - 1);
    static constexpr float kInvScale = 1.0f / static_cast<float>(uint32_t{1} << kFracBits);
};

// ------------------------------------------------------------
//  Q3_12Storage — int16 上的 Q3.12 运算域（带一阶误差反馈）
// ------------------------------------------------------------
struct Q3_12Storage : Q16Arithmetic<12> {
    using Line = QuantizedDelayLine<Q3_12Storage>;

    static constexpr const char* kName = "Q3.12";
    static constexpr const char* kDescription = "Q3.12 arithmetic (int16) + 1st-order EF";
};

// ------------------------------------------------------------
//  QuantizedDelayLine — 整数延迟、无插值的单声道延迟线
//
//  容量向上取到 2 的幂；收发都是运算域里的原始值，不做任何换算。
// ------------------------------------------------------------
template <typename Storage>
class QuantizedDelayLine {
public:
    using Value = typename Storage::Value;

    /// @brief 最小可用延迟样本数
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
        buffer_.assign(size_, Value{});
    }

    /// @brief 写入一个样本
    void Push(Value x) noexcept {
        buffer_[wpos_] = x;
        wpos_ = (wpos_ + 1) & mask_;
    }

    /// @brief 读取 delay_samples 次 Push 之前写入的样本
    /// @param[in] delay_samples 延迟样本数，取值须落在 [1, Capacity()]
    Value GetBeforePush(size_t delay_samples) const noexcept {
        return buffer_[(wpos_ - delay_samples) & mask_];
    }

    /// @brief 可用延迟样本数上限
    size_t Capacity() const noexcept { return size_; }
private:
    std::vector<Value> buffer_;
    size_t wpos_ = 0;
    size_t mask_ = 0;
    size_t size_ = 0;
};

} // namespace echo_demo
