#pragma once

#include <qwqdsp/colormap/magma.hpp>
#include "raylib.h"

// ----------------------------------------
// magma 色表 (matplotlib magma, 256 级)
// 适配器: 将 qwqdsp_colormap::Magma 转为 raylib::Color
// ----------------------------------------

namespace detail {
    constexpr std::array<Color, 256> MakeMagmaTable() {
        std::array<Color, 256> table{};
        for (int i = 0; i < 256; ++i) {
            table[i] = Color{
                qwqdsp_colormap::Magma::kTable[i][0],
                qwqdsp_colormap::Magma::kTable[i][1],
                qwqdsp_colormap::Magma::kTable[i][2],
                255
            };
        }
        return table;
    }
}

struct MagmaColormap {
    static constexpr auto kTable = detail::MakeMagmaTable();
};
