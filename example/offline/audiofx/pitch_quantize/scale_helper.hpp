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
 *
 * 调性分类：
 *   核心基础 — Major (Ionian), Minor (Aeolian)
 *   教会调式 — Dorian, Phrygian, Lydian, Mixolydian, Locrian
 *   古典小调 — Harmonic Minor, Melodic Minor (Jazz Minor)
 *   五声特殊 — Major Pentatonic, Minor Pentatonic, Blues, Whole Tone
 *   半音阶   — Chromatic
 */
struct ScaleHelper {
    enum class Type : uint8_t {
        // 核心基础
        kMajor,         ///< 自然大调 (Ionian)
        kMinorNatural,  ///< 自然小调 (Aeolian)
        // 教会调式
        kDorian,
        kPhrygian,
        kLydian,
        kMixolydian,
        kLocrian,
        // 古典小调
        kMinorHarmonic, ///< 和声小调
        kMinorMelodic,  ///< 旋律小调（爵士小调，上行）
        // 五声与特殊
        kMajorPentatonic,
        kMinorPentatonic,
        kBlues,
        kWholeTone,
        // 半音阶
        kChromatic,     ///< 半音阶（全部 12 个音）
    };

    static constexpr std::array<int, 7> kMajorIntervals{0, 2, 4, 5, 7, 9, 11};
    static constexpr std::array<int, 7> kMinorNaturalIntervals{0, 2, 3, 5, 7, 8, 10};
    static constexpr std::array<int, 7> kDorianIntervals{0, 2, 3, 5, 7, 9, 10};
    static constexpr std::array<int, 7> kPhrygianIntervals{0, 1, 3, 5, 7, 8, 10};
    static constexpr std::array<int, 7> kLydianIntervals{0, 2, 4, 6, 7, 9, 11};
    static constexpr std::array<int, 7> kMixolydianIntervals{0, 2, 4, 5, 7, 9, 10};
    static constexpr std::array<int, 7> kLocrianIntervals{0, 1, 3, 5, 6, 8, 10};
    static constexpr std::array<int, 7> kMinorHarmonicIntervals{0, 2, 3, 5, 7, 8, 11};
    static constexpr std::array<int, 7> kMinorMelodicIntervals{0, 2, 3, 5, 7, 9, 11};
    static constexpr std::array<int, 5> kMajorPentatonicIntervals{0, 2, 4, 7, 9};
    static constexpr std::array<int, 5> kMinorPentatonicIntervals{0, 3, 5, 7, 10};
    static constexpr std::array<int, 6> kBluesIntervals{0, 3, 5, 6, 7, 10};
    static constexpr std::array<int, 6> kWholeToneIntervals{0, 2, 4, 6, 8, 10};
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
        case Type::kDorian:
            intervals = kDorianIntervals;
            break;
        case Type::kPhrygian:
            intervals = kPhrygianIntervals;
            break;
        case Type::kLydian:
            intervals = kLydianIntervals;
            break;
        case Type::kMixolydian:
            intervals = kMixolydianIntervals;
            break;
        case Type::kLocrian:
            intervals = kLocrianIntervals;
            break;
        case Type::kMinorHarmonic:
            intervals = kMinorHarmonicIntervals;
            break;
        case Type::kMinorMelodic:
            intervals = kMinorMelodicIntervals;
            break;
        case Type::kMajorPentatonic:
            intervals = kMajorPentatonicIntervals;
            break;
        case Type::kMinorPentatonic:
            intervals = kMinorPentatonicIntervals;
            break;
        case Type::kBlues:
            intervals = kBluesIntervals;
            break;
        case Type::kWholeTone:
            intervals = kWholeToneIntervals;
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
     * @brief 调性描述
     */
    static const char* TypeName(Type type) noexcept {
        switch (type) {
        case Type::kMajor:
            return "Major";
        case Type::kMinorNatural:
            return "Minor";
        case Type::kDorian:
            return "Dorian";
        case Type::kPhrygian:
            return "Phrygian";
        case Type::kLydian:
            return "Lydian";
        case Type::kMixolydian:
            return "Mixolydian";
        case Type::kLocrian:
            return "Locrian";
        case Type::kMinorHarmonic:
            return "Harmonic Minor";
        case Type::kMinorMelodic:
            return "Melodic Minor";
        case Type::kMajorPentatonic:
            return "Major Pentatonic";
        case Type::kMinorPentatonic:
            return "Minor Pentatonic";
        case Type::kBlues:
            return "Blues";
        case Type::kWholeTone:
            return "Whole Tone";
        case Type::kChromatic:
            return "Chromatic";
        }
        return "Unknown";
    }

    static constexpr size_t kNumTypes = 14;

    static constexpr const char* kNoteNames[12] = {"C", "C#", "D", "D#", "E", "F",
                                                    "F#", "G", "G#", "A", "A#", "B"};
};
