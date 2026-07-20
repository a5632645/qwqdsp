#pragma once

#include <array>
#include <cstdint>
#include <span>

// ------------------------------------------------------------
// ScaleHelper — 调性辅助工具
// ------------------------------------------------------------
/**
 * @brief 提供音阶级别上的调性判断与 MIDI 音符查找。
 *
 * 音级 (note class) 映射: 0=C, 1=C#, 2=D, ..., 11=B。
 * 通过位掩码 (uint16_t, 低 12 位) 表示哪些音级被允许。
 */
struct ScaleHelper {
    enum class Type : uint8_t {
        kMajor,         ///< 大调音阶
        kMinorNatural,  ///< 自然小调
        kMinorHarmonic, ///< 和声小调
        kMinorMelodic,  ///< 旋律小调（上行）
        kChromatic,     ///< 半音阶（全部 12 个音）
    };

    // 从根音开始的半音间隔
    static constexpr std::array<int, 7> kMajorIntervals{0, 2, 4, 5, 7, 9, 11};
    static constexpr std::array<int, 7> kMinorNaturalIntervals{0, 2, 3, 5, 7, 8, 10};
    static constexpr std::array<int, 7> kMinorHarmonicIntervals{0, 2, 3, 5, 7, 8, 11};
    static constexpr std::array<int, 7> kMinorMelodicIntervals{0, 2, 3, 5, 7, 9, 11};
    static constexpr std::array<int, 12> kChromaticIntervals{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    /**
     * @brief 生成指定调性的音级掩码
     * @param root  根音音级 (0=C, 1=C#, …, 11=B)
     * @param type  音阶类型
     * @return 低 12 位表示允许的音级
     */
    static uint16_t MakeMask(int root, Type type) noexcept {
        uint16_t mask = 0;
        std::span<const int> intervals;
        switch (type) {
        case Type::kMajor:
            intervals = kMajorIntervals;
            break;
        case Type::kMinorNatural:
            intervals = kMinorNaturalIntervals;
            break;
        case Type::kMinorHarmonic:
            intervals = kMinorHarmonicIntervals;
            break;
        case Type::kMinorMelodic:
            intervals = kMinorMelodicIntervals;
            break;
        case Type::kChromatic:
            intervals = kChromaticIntervals;
            break;
        }
        for (int iv : intervals) {
            mask |= static_cast<uint16_t>(1u << ((root + iv) % 12));
        }
        return mask;
    }

    /**
     * @brief 判断 MIDI 音级是否在允许集合中
     */
    static bool IsAllowed(int note_class, uint16_t mask) noexcept {
        return (mask >> static_cast<uint16_t>(note_class % 12)) & 1u;
    }

    /**
     * @brief 调性描述（仅用于日志输出）
     */
    static const char* TypeName(Type type) noexcept {
        switch (type) {
        case Type::kMajor:
            return "Major";
        case Type::kMinorNatural:
            return "Minor(Natural)";
        case Type::kMinorHarmonic:
            return "Minor(Harmonic)";
        case Type::kMinorMelodic:
            return "Minor(Melodic)";
        case Type::kChromatic:
            return "Chromatic";
        }
        return "Unknown";
    }

    static constexpr const char* kNoteNames[12] = {"C", "C#", "D", "D#", "E", "F",
                                                    "F#", "G", "G#", "A", "A#", "B"};
};
