#pragma once

#include <qwqdsp/colormap/jet.hpp>
#include "raylib.h"

// ----------------------------------------
// jet 色表 (matplotlib jet, 256 级)
// 适配器: 将 qwqdsp_colormap::Jet 转为 raylib::Color
// ----------------------------------------

namespace detail {
    constexpr std::array<Color, 256> MakeJetTable() {
        std::array<Color, 256> table{};
        for (int i = 0; i < 256; ++i) {
            table[i] = Color{
                qwqdsp_colormap::Jet::kTable[i][0],
                qwqdsp_colormap::Jet::kTable[i][1],
                qwqdsp_colormap::Jet::kTable[i][2],
                255
            };
        }
        return table;
    }
}

struct JetColormap {
    static constexpr auto kTable = detail::MakeJetTable();
};
