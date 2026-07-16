#pragma once

#include <array>
#include "raylib.h"

// ------------------------------------------------------------
// 模板适配器: 将 qwqdsp_colormap 中的 RGB 表转为 raylib::Color
// 用法: using MagmaColormap = ColormapAdapter<qwqdsp_colormap::Magma>;
// ------------------------------------------------------------
template <typename Src>
struct ColormapAdapter {
    static constexpr std::array<Color, 256> MakeTable() {
        std::array<Color, 256> table{};
        for (int i = 0; i < 256; ++i) {
            table[i] = Color{
                Src::kTable[i][0],
                Src::kTable[i][1],
                Src::kTable[i][2],
                255
            };
        }
        return table;
    }
    static constexpr auto kTable = MakeTable();
};
