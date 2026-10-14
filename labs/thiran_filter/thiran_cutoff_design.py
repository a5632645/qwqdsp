# -*- coding: utf-8 -*-
"""
thiran_cutoff_design.py
=======================

数字 Thiran 全极点低通的「指定截止频率」设计：固定阶数 n，给定 -3dB 截止频率 wc，
求满足该截止频率的系数。

背景（文献结论，见 thiran_cutoff_design.md）
-------------------------------------------
文献里的 Thiran 设计只有一个参数（延迟 δ / β = 2δ）：
    H(z) = 1/A(z),  A(z) = Σ_k a_k z^-k,
    a_k = (-1)^k C(n,k) ∏_{i=0..n} (β+i)/(β+k+i),   a_0 = 1
δ 同时决定直流群延迟（= δ）与幅度响应（= 截止频率），不能各自指定 ——
Laakso 的综述明确写 "the magnitude response of the all-pole lowpass filter cannot be
controlled"。所以「指定截止频率」= 在 β 上解一个标量方程，系数公式本身不变。

推导出来的方程
--------------
-3dB 定义：|H(e^jwc)| = |H(1)|/√2  <=>  |A(e^jwc)|² = 2·A(1)²，即

    F(β, wc) = Σ_{k,l} a_k(β) a_l(β) cos((k-l)wc) - 2·(Σ_k a_k(β))²  = 0

把分母通分后得到 β 的 2n 次多项式方程（记 u = 1-cos wc = 2sin²(wc/2)）：

    n=1:  u·β² + 2u·β - 2 = 0
          β = -1 + sqrt(1 + 2/u)
    n=2:  u²·β⁴ + 8u²·β³ + u(6+19u)·β² + 12u(2+u)·β - 36 = 0   （四次，取正根）

一般 n 用数值求根：F 在 β ∈ (0,∞) 上恰好一个符号变化（实测 n=1…10、25 个 wc 都是 1 次），
唯一根；β→0⁺ 时 F→-1，β→∞ 时 F→+∞，所以"从 0 往外扩括号 + Brent"总能收敛。
"""

import math
from math import comb

import numpy as np
from scipy.optimize import brentq

SQRT2 = math.sqrt(2.0)


# ------------------------------------------------------------
# 系数与截止频率
# ------------------------------------------------------------

def thiran_allpole_coef(beta, n):
    """全极点 Thiran 的分母系数 a[0..n]（a_0 = 1），beta = 2*delta。"""
    a = np.ones(n + 1)
    for k in range(1, n + 1):
        prod = 1.0
        for i in range(n + 1):
            prod *= (beta + i) / (beta + k + i)
        a[k] = (-1) ** k * comb(n, k) * prod
    return a


def residual(beta, n, wc):
    """-3dB 条件残差 F(β, wc)：|A(e^jwc)|² - 2·A(1)²，随 β 单调（恰好一次过零）。"""
    a = thiran_allpole_coef(beta, n)
    mag2 = float(np.abs(a @ np.exp(-1j * wc * np.arange(n + 1))) ** 2)
    return mag2 - 2.0 * float(a.sum()) ** 2


def _bracket(n, wc):
    """向外扩出 F<0 / F>0 的括号。"""
    lo = 1e-12
    hi = 1e-6
    for _ in range(200):
        if residual(hi, n, wc) > 0.0:
            return lo, hi
        lo = hi
        hi *= 2.0
    raise RuntimeError(f"n={n}, wc={wc}: 未找到括号")


def solve_beta(n, wc):
    """解 F(β, wc) = 0，返回 β = 2δ。"""
    lo, hi = _bracket(n, wc)
    return brentq(lambda b: residual(b, n, wc), lo, hi, xtol=1e-15, rtol=1e-14)


def beta_closed_form_n1(wc):
    """n=1 的闭式解：β = -1 + sqrt(1 + 2/(1-cos wc))。"""
    u = 1.0 - math.cos(wc)
    return -1.0 + math.sqrt(1.0 + 2.0 / u)


def beta_closed_form_n2(wc, lo=1e-12, hi=1e6):
    """n=2 的正根：u²β⁴+8u²β³+u(6+19u)β²+12u(2+u)β-36 = 0，u = 1-cos wc。"""
    u = 1.0 - math.cos(wc)

    def poly(b):
        return u * u * b**4 + 8.0 * u * u * b**3 + u * (6.0 + 19.0 * u) * b**2 + 12.0 * u * (2.0 + u) * b - 36.0

    return brentq(poly, lo, hi, xtol=1e-15, rtol=1e-14)


def measure_cutoff(beta, n):
    """用二分法测 -3dB 频率（|A| 在 [0,π] 上单调增，故可二分），返回 0 表示未过 -3dB。"""
    a = thiran_allpole_coef(beta, n)
    dc = float(a.sum())
    g = lambda w: float(np.abs(a @ np.exp(-1j * w * np.arange(n + 1)))) / dc
    if g(math.pi) < SQRT2:
        return math.nan
    return brentq(lambda w: g(w) - SQRT2, 0.0, math.pi, xtol=1e-15, rtol=1e-14)


# ------------------------------------------------------------
# 自检
# ------------------------------------------------------------

def check_closed_forms():
    print("--- 闭式解 vs 数值解 ---")
    for wc_pi in (0.05, 0.2, 0.5, 0.8, 0.95):
        wc = wc_pi * math.pi
        b_num = solve_beta(1, wc)
        b_cf = beta_closed_form_n1(wc)
        print(f"  n=1 wc/pi={wc_pi}: 闭式 β={b_cf:.12f}  数值 β={b_num:.12f}  差 {abs(b_cf-b_num):.2e}")
    for wc_pi in (0.05, 0.2, 0.5, 0.8, 0.95):
        wc = wc_pi * math.pi
        b_num = solve_beta(2, wc)
        b_cf = beta_closed_form_n2(wc)
        print(f"  n=2 wc/pi={wc_pi}: 四次正根 β={b_cf:.12f}  数值 β={b_num:.12f}  差 {abs(b_cf-b_num):.2e}")


def check_uniqueness():
    print("--- 唯一性：F(β) 在 β∈[1e-4,1e6] 上的符号变化次数 ---")
    betas = np.logspace(-4, 6, 2001)
    for n in range(1, 11):
        counts = []
        for wc_pi in np.linspace(0.02, 0.98, 25):
            wc = wc_pi * math.pi
            v = np.array([residual(b, n, wc) for b in betas])
            v = v[np.abs(v) > 1e-12]
            counts.append(int(np.sum(np.diff(np.sign(v)) != 0)))
        print(f"  n={n:2d}: 符号变化次数 {min(counts)}…{max(counts)}（<=1 即唯一根）")


def check_design_sweep():
    print("--- 设计扫描：解 β 后实测 -3dB 频率误差 ---")
    worst = 0.0
    for n in (1, 2, 3, 4, 6, 8, 10):
        errs = []
        for wc_pi in np.linspace(0.03, 0.97, 19):
            wc = wc_pi * math.pi
            beta = solve_beta(n, wc)
            got = measure_cutoff(beta, n) / math.pi
            errs.append(abs(got - wc_pi))
        worst = max(worst, max(errs))
        print(f"  n={n:2d}: max|wc_measured - wc_target|/pi = {max(errs):.2e}")
    print(f"  总最大误差 = {worst:.2e}")


def table():
    print("--- 设计表：δ = β/2（同时就是直流群延迟）与直流增益 1/A(1) ---")
    header = "  n\\wc/pi " + " ".join(f"{w:>9.2f}" for w in (0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9))
    print(header)
    for n in (1, 2, 3, 4, 6, 8, 10):
        cells = []
        for wc_pi in (0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9):
            beta = solve_beta(n, wc_pi * math.pi)
            cells.append(f"{beta/2:9.4f}")
        print(f"  {n:>2}     " + " ".join(cells))

    print("--- 同一批设计的直流增益 1/A(1)（该族里截止频率与增益绑定，无法各自指定）---")
    print(header)
    for n in (1, 2, 3, 4, 6, 8, 10):
        cells = []
        for wc_pi in (0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9):
            beta = solve_beta(n, wc_pi * math.pi)
            cells.append(f"{1.0/thiran_allpole_coef(beta, n).sum():9.4f}")
        print(f"  {n:>2}     " + " ".join(cells))


def conditioning(beta, n):
    """浮点条件数：Σ|a_k|/|A(1)|。A(1) 是交替和，β 越大消位越狠，精度损失约 eps*cond。"""
    a = thiran_allpole_coef(beta, n)
    return float(np.abs(a).sum() / abs(a.sum()))


def asymptotics():
    print("--- 小截止频率渐近：K(n) = tan(wc/2)·δ（δ→∞ 时趋于常数）---")
    print("    只列数值可信的档位，括号里是条件数 Σ|a_k|/|A(1)|（精度损失 ≈ eps·cond）")
    for n in (1, 2, 3, 4, 6, 8, 10):
        cells = []
        for wc_pi in (0.05, 0.10, 0.20):
            wc = wc_pi * math.pi
            beta = solve_beta(n, wc)
            cells.append(f"{math.tan(wc/2)*beta/2:7.4f}({conditioning(beta, n):.0e})")
        print(f"  n={n:2d}: " + "  ".join(cells))
    print("  n=1 的闭式极限：wc→0 时 K→1/2（= tan(wc/2)*δ 的极限）")


def conditioning_limits():
    print("--- 数值可用范围：条件数随 wc 的变化（eps·cond 是可达精度）---")
    for n in (4, 8, 10):
        cells = []
        for wc_pi in (0.3, 0.1, 0.05, 0.02):
            beta = solve_beta(n, wc_pi * math.pi)
            cells.append(f"wc/pi={wc_pi:g}: cond={conditioning(beta, n):.1e}")
        print(f"  n={n:2d}: " + "  ".join(cells))


if __name__ == "__main__":
    check_closed_forms()
    print()
    check_uniqueness()
    print()
    check_design_sweep()
    print()
    asymptotics()
    print()
    conditioning_limits()
    print()
    table()
