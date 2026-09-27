#pragma once
#include <matplot/matplot.h>

#include <filesystem>
#include <string>

namespace qwqdsp_support {
/**
 * @brief 关掉 matplot++ 的中间重画，并把它兜底的输出终端从交互窗口改到临时 PNG
 *
 * 慢与"闪窗"的根源在 `figure_type::touch()`：它只要图被改动就立刻 `draw()` 一次整幅图，
 * 而 `plot()` / `semilogy()` / `title()` / `xlim()` / `legend()` / `grid()` ... 每个都触发
 * `touch()`。于是画 N 条曲线就有 N+ 次重画：先出现"只有曲线、没有图例"的帧，再空白，
 * 再完整帧——默认终端在 Windows 上是 `wxt`（其余平台 `qt`），每帧都要建/绘一次窗口。
 * 实测 `rls_lp_spectrum`（2 条曲线）28 s、`burg_lp_spectrum`（5 条曲线）64 s。
 *
 * ① `gcf(true)` 打开静默模式，`touch()` 不再 draw：只有显式 `matplot::save()`/`show()` 才画。
 * ② 再把后端输出指向临时 PNG（`set terminal pngcairo` + `set output "..."`）：万一真去画了
 *    （比如调了 `show()`），也是写文件而不是弹窗。
 *
 * @note 必须在**第一次 plot 之前**调用，否则第一帧已经画过了
 */
static inline void MakeMatplotHeadless() {
    auto const fig = matplot::gcf(true); // quiet mode：不再每次改图都重画
    auto const path = (std::filesystem::temp_directory_path() / "qwqdsp_matplot_tmp.png").string();
    fig->backend()->output(path); // 兜底：真画也画到文件，不弹窗
}
} // namespace qwqdsp_support
