"""论文 §IV-D 末尾：计算代价对比与「代价随极点对数线性增长」的验证。

论文在 3 代 i5 + MATLAB（逐采样循环）上给出：10 s / 44.1 kHz 音频，
AA-FIR-1 ≈ 0.6 s、AA-IIR-1 ≈ 0.4 s、AA-IIR-2 ≈ 2 s（5 对极点）。
本复现是 numpy 批量实现，**绝对时间不可与论文对比**；可验证的是论文的
结构性结论：AA-IIR 的代价随（复共轭极点对数量）线性增长。

用法: python qwqdsp/labs/adaa_iir/exp_cost.py
"""
from __future__ import annotations

import time

import numpy as np

import aaiir as A

FS = 44100.0
DURATION = 10.0          # 与论文一致的时长
F0 = 1000.0
REPEAT = 3


def timed(fn):
    best = float("inf")
    for _ in range(REPEAT):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def main() -> int:
    tab = A.escalation_ii_w3_table()
    saw = A.saw_table()
    n = int(FS * DURATION)
    x = A.phase_ramp(n, F0, FS)

    print(f"目标：{DURATION:.0f} s / {FS:.0f} Hz 音频（Python 3 + numpy 批量实现，"
          f"{REPEAT} 次取最快；绝对时间不可与论文的 MATLAB 数字对比）")
    print(f"{'方法':<28}{'耗时 [s]':>10}{'µs/采样':>12}")
    cases = {
        "trivial（锯齿）": lambda: A.trivial(saw, n, F0, FS),
        "AA-FIR-1（锯齿）": lambda: A.aa_fir(saw, n, F0, FS, 1),
        "AA-FIR-2（锯齿）": lambda: A.aa_fir(saw, n, F0, FS, 2),
        "OVS-8（Escalation）": lambda: A.oversampled(tab, n, F0, FS, 8),
        "AA-IIR-1（Escalation）": lambda: A.aa_iir_osc(x, A.aa_iir_1(), tab),
        "AA-IIR-2（Escalation）": lambda: A.aa_iir_osc(x, A.aa_iir_2(), tab),
    }
    for name, fn in cases.items():
        t = timed(fn)
        print(f"{name:<28}{t:>10.2f}{t / n * 1e6:>12.3f}")

    print("\nAA-IIR 代价 vs 极点对数（同规格 Chebyshev II，仅改阶数）：")
    print(f"{'阶数':>6}{'极点对':>8}{'耗时 [s]':>12}{'µs/采样':>12}{'相对 1 对':>12}")
    base = None
    for order in (2, 4, 6, 8, 10):
        filt = A.aa_filter_cheby2(order, A.AAIIR2_RS, A.AAIIR2_WN)
        pairs = len(filt[0])
        t = timed(lambda f=filt: A.aa_iir_osc(x, f, tab))
        base = base or t
        print(f"{order:>6}{pairs:>8}{t:>12.2f}{t / n * 1e6:>12.3f}{t / base:>12.2f}")
    print("（线性增长的判据：耗时倍数 ≈ 极点对数倍数；批量实现里隐含常数开销会让小阶数偏大）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
