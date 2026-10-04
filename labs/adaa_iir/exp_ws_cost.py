"""ADAA-IIR waveshaper 的**代价结构**：为什么它不能像 ADAA-LUT 那样预计算表。

要点
----
* FIR 版（ADAA-LUT）：原函数 F₀/F₁ 与输入无关，可预先数值积分成表 → 每次只需插值读取；
* IIR 版：积分核 ``e^{β(1-u)}`` 里的 β 是滤波器极点，但**积分区间由输入逐样本差分 Δ 决定**，
  因此不能预计算「与 Δ 无关的原函数表」。若把 f 拟合成 k 段分段线性，则每采样的代价 ∝
  **该采样区间跨过的段数**（= |Δ| / 段宽），这正是 [19] §III-F 说的 ``l = i_max − i_min``。

本脚本量化：跨段数统计、耗时、以及对 f 做 K 段拟合时代价如何随 K 增长。

用法: python qwqdsp/labs/adaa_iir/exp_ws_cost.py
"""
from __future__ import annotations

import time

import numpy as np

import aaiir as A
import ws as W

FS = 44100.0
F0 = 1000.0
GAIN = 10.0
DURATION = 5.0
REPEAT = 3


def timed(fn):
    best = float("inf")
    for _ in range(REPEAT):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def crossings(x, f: W.Pwl):
    """每采样区间跨过的段数（不含两端的半无限段）。"""
    i0 = np.searchsorted(f.X, x[:-1], side="right")
    i1 = np.searchsorted(f.X, x[1:], side="right")
    return np.abs(i1 - i0) + 1


def exp_calls_per_sample(f: W.Pwl, pairs: int) -> float:
    """每采样需要的 exp 次数估计 = 极点数 × 跨段数（每段一次）。"""
    return 0.0 * pairs  # 占位，main 里按实测跨段数计算


def main() -> int:
    n = int(FS * DURATION)
    x = GAIN * np.sin(2.0 * np.pi * F0 * np.arange(n) / FS)
    clip = W.hard_clip()

    print(f"输入：{GAIN:.0f}·sin({F0:.0f} Hz)，{DURATION:.0f} s @ {FS:.0f} Hz")
    print("\n== 每采样跨段数（决定 IIR 路线的实时代价） ==")
    print(f"{'非线性':<26}{'段数':>8}{'平均跨段数':>12}{'最大':>8}")
    variants = [("hard_clip（2 折点）", clip)]
    for k in (8, 32, 128):
        variants.append((f"hard_clip 的 {k} 段拟合", W.fit_from_callable(
            lambda v: np.clip(v, -1.0, 1.0), -10.0, 10.0, k)))
    stats = {}
    for name, f in variants:
        c = crossings(x, f)
        stats[name] = c.mean()
        print(f"{name:<26}{len(f.m) + 2:>8}{c.mean():>12.2f}{c.max():>8d}")
    print("（硬削波本身只有 2 个折点，跨段数 ≈ 1–3；拟合得越细，跨段数按 K 线性增长）")

    print("\n== 耗时（Python 批量实现，3 次取最快；看倍数不要看绝对值） ==")
    print(f"{'方法':<12}{'耗时 [s]':>10}{'µs/采样':>10}{'exp 次数/采样':>14}")
    rows = {
        "trivial": lambda: W.trivial(x, clip),
        "AA-FIR-1": lambda: W.aa_fir1(x, clip),
        "AA-FIR-2": lambda: W.aa_fir2(x, clip),
        "OVS-8": lambda: W.oversampled(x, 8, clip),
        "AA-IIR-1": lambda: W.aa_iir_pwl(x, W.aa_iir_1(), clip),
        "AA-IIR-2": lambda: W.aa_iir_pwl(x, W.aa_iir_2(), clip),
    }
    nc = stats["hard_clip（2 折点）"]
    n_runs = {
        "trivial": 0.0, "AA-FIR-1": 0.0, "AA-FIR-2": 0.0, "OVS-8": 0.0,
        "AA-IIR-1": 1 * nc, "AA-IIR-2": 5 * nc,
    }
    for name, fn in rows.items():
        t = timed(fn)
        print(f"{name:<12}{t:>10.2f}{t / n * 1e6:>10.1f}{n_runs[name]:>14.1f}")

    print("\n== 对 tanh 做 K 段拟合：代价随 K 增长 ==")
    print(f"{'K 段':>6}{'跨段数':>10}{'AA-IIR-1 耗时 [s]':>20}{'µs/采样':>10}")
    for k in (8, 32, 128):
        f = W.fit_from_callable(np.tanh, -10.0, 10.0, k)
        c = crossings(x, f).mean()
        t = timed(lambda f=f: W.aa_iir_pwl(x, W.aa_iir_1(), f))
        print(f"{k:>6}{c:>10.2f}{t:>20.2f}{t / n * 1e6:>10.1f}")

    print("\n结论：IIR 路线的每采样代价 ∝ 跨段数 ∝ |Δ|/段宽。要让代价与信号无关，只能")
    print("      (a) 用少量折点（硬削波/波折叠这类本来就是 2 折点）；")
    print("      (b) 用*数值求积*路线（固定点数，与折点数无关，见 ws.mean_integral_numeric）；")
    print("      (c) 或接受「密集波表用 FIR/ADAA-LUT、解析或粗分段用 IIR」的分工。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
