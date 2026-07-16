#pragma once

#include <qwqdsp/colormap/viridis.hpp>
#include "raylib.h"

// ----------------------------------------
// viridis 色表 (matplotlib viridis, 256 级)
// 适配器: 将 qwqdsp_colormap::Viridis 转为 raylib::Color
// ----------------------------------------

namespace detail {
    constexpr std::array<Color, 256> MakeViridisTable() {
        std::array<Color, 256> table{};
        for (int i = 0; i < 256; ++i) {
            table[i] = Color{
                qwqdsp_colormap::Viridis::kTable[i][0],
                qwqdsp_colormap::Viridis::kTable[i][1],
                qwqdsp_colormap::Viridis::kTable[i][2],
                255
            };
        }
        return table;
    }
}

struct ViridisColormap {
    static constexpr auto kTable = detail::MakeViridisTable();
};
