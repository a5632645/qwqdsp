# -*- coding: utf-8 -*-
"""
zero_pole_lowpass.py
====================

零极点低通的两种直流最大平坦群延迟设计：

1. 线性相位零点 + Thiran 全极点（闭式）:
       B(z) = (1 + z^-1)^L
       A(z) = Thiran 全极点分母，参数 delta = tau0 - L/2
       H(z) = B(z) / A(z)
   分子群延迟恒为 L/2，分母的 2M 阶平坦性原样保留。

2. Zhang 线性方程组（平坦幅度 + 平坦延迟）:
   给定分子阶数 N、分母阶数 M、z=-1 零点数 L、直流延迟 tau0，
   解 N+M+1 个线性方程。N=M、L=0 时退化为 Thiran 全通。

脚本做三件事：固定 delta 扫描零点数、对比级联与 Zhang、数值自检。
图写到 output/：

    zero_pole_scan_L.png   固定 delta=1、M=4，扫 L（群延迟误差应完全重合）
    zero_pole_compare.png  同一分子下，级联极点 vs Zhang 极点（幅度/延迟取舍）
"""

import math
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from explore_scan import W, WX, ap_coef, group_delay, lp_coef, max_pole_radius

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei"]
matplotlib.rcParams["axes.unicode_minus"] = False

OUT_DIR = Path(__file__).parent / "output"


# ------------------------------------------------------------
# 设计器
# ------------------------------------------------------------

def binom_row(n):
    return np.array([math.comb(n, k) for k in range(n + 1)], dtype=float)


def cascade_coef(poles, zeros, tau0, normalize=True):
    """z=-1 的 zeros 个零点 + Thiran 全极点，直流群延迟 tau0。"""
    delta = tau0 - 0.5 * zeros
    a = lp_coef(poles, delta)
    b = binom_row(zeros)
    if normalize:
        b *= a.sum() / b.sum()
    return b, a, delta


def zhang_coef(N, M, tau0, L):
    """Zhang 线性方程组：分子阶数 N、分母阶数 M、z=-1 零点数 L。"""
    K = N + M + 1 - L
    rows = np.zeros((K + L, N + 1 + M))
    rhs = np.zeros(K + L)
    for i in range(K):
        for n in range(N + 1):
            rows[i, n] = (n - tau0) ** i
        for m in range(1, M + 1):
            rows[i, N + m] = -(m ** i)
        rhs[i] = 1.0 if i == 0 else 0.0
    for i in range(L):
        for n in range(N + 1):
            rows[K + i, n] = ((-1) ** n) * (n ** i)
        rhs[K + i] = 0.0
    x = np.linalg.solve(rows, rhs)
    return x[:N + 1], np.concatenate(([1.0], x[N + 1:]))


def response(b, a, w):
    w = np.atleast_1d(w)
    z = np.exp(-1j * w)
    num = np.sum(b[None, :] * z[:, None] ** np.arange(len(b))[None, :], axis=1)
    den = np.sum(a[None, :] * z[:, None] ** np.arange(len(a))[None, :], axis=1)
    return num / den


def delay_error(b, a, tau0, w):
    with np.errstate(divide="ignore", invalid="ignore"):
        return group_delay(b, a, np.atleast_1d(w)) - tau0


def flatness_slope(b, a, tau0, w1=0.05, w2=0.10):
    e1 = abs(delay_error(b, a, tau0, w1)[0])
    e2 = abs(delay_error(b, a, tau0, w2)[0])
    return math.log(e2 / e1) / math.log(w2 / w1)


# ------------------------------------------------------------
# 自检与表
# ------------------------------------------------------------

def check_cascade():
    print("--- 级联：固定 delta=1、M=4，扫描 z=-1 零点数 L ---")
    print("  L  tau0  max|p|   err@0.1pi    |H|@0.1pi dB  |H|@0.5pi dB")
    poles = 4
    delta = 1.0
    for L in range(5):
        tau0 = delta + 0.5 * L
        b, a, _ = cascade_coef(poles, L, tau0)
        err = abs(delay_error(b, a, tau0, 0.1 * np.pi)[0])
        m1 = 20.0 * math.log10(abs(response(b, a, 0.1 * np.pi)[0]))
        m2 = 20.0 * math.log10(abs(response(b, a, 0.5 * np.pi)[0]))
        print(f"  {L:>2}  {tau0:4.1f}  {max_pole_radius(a):6.4f}   {err:10.3e}    {m1:12.5f}    {m2:12.5f}")

    b, a, delta = cascade_coef(4, 2, 2.0)
    slope = flatness_slope(b, a, 2.0)
    print(f"  平坦阶数实测：M=4, L=2, tau0=2, delta=1 -> {slope:.3f}（期望 2M=8）")
    print(f"  B(-1) = {sum(((-1) ** k) * b[k] for k in range(len(b))):.3e}")
    print()


def check_zhang():
    print("--- Zhang：N=M、L=0 应精确退化为 Thiran 全通 ---")
    for n, delta in ((2, 0.5), (3, 0.4), (3, -0.7)):
        b, a = zhang_coef(n, n, n + delta, 0)
        a_ref = ap_coef(n, delta)
        b_ref = a_ref[::-1]
        print(f"  n={n} delta={delta:+.1f}: max|b-b_ref|={np.max(np.abs(b - b_ref)):.2e}"
              f"  max|a-a_ref|={np.max(np.abs(a - a_ref)):.2e}")
    print()

    print("--- Zhang：低通例子（N=8, M=3, L=6, K=6）---")
    for tau0 in (5.9, 5.2, 4.5, 3.8):
        b, a = zhang_coef(8, 3, tau0, 6)
        err = abs(delay_error(b, a, tau0, 0.1 * np.pi)[0])
        mag = 20.0 * math.log10(abs(response(b, a, 0.1 * np.pi)[0]))
        print(f"  tau0={tau0:3.1f}: max|p|={max_pole_radius(a):.4f}  err@0.1pi={err:.3e}"
              f"  |H|@0.1pi={mag:+.4f} dB")
    print()

    print("--- Zhang：N=7, M=3, L=3, K=8（K 为偶数，群延迟误差阶数 = K）---")
    N, M, L, tau0 = 7, 3, 3, 5.2
    b, a = zhang_coef(N, M, tau0, L)
    slope = flatness_slope(b, a, tau0, 0.05, 0.10)
    print(f"  max|p|={max_pole_radius(a):.4f}  实测幂次={slope:.2f}  预期 {N + M + 1 - L}")
    print(f"  |H|@0.1pi={20.0 * math.log10(abs(response(b, a, 0.1 * np.pi)[0])):+.5f} dB")
    print()

    print("--- Zhang：N=3, M=1, L=2, K=3（K 为奇数，群延迟误差阶数 = K-1）---")
    b, a = zhang_coef(3, 1, 2.0, 2)
    slope = flatness_slope(b, a, 2.0, 0.02, 0.04)
    print(f"  max|p|={max_pole_radius(a):.4f}  实测幂次={slope:.2f}  预期 2")
    print()


# ------------------------------------------------------------
# 图
# ------------------------------------------------------------

def plot_scan_L():
    fig, (ax_m, ax_g) = plt.subplots(1, 2, figsize=(11, 4.2))
    poles, delta = 4, 1.0
    for L in range(5):
        tau0 = delta + 0.5 * L
        b, a, _ = cascade_coef(poles, L, tau0)
        label = f"L={L}, τ0={tau0:g}"
        ax_m.plot(WX, 20.0 * np.log10(np.maximum(np.abs(response(b, a, W)), 1e-12)), label=label)
        ax_g.plot(WX[:-1], delay_error(b, a, tau0, W[:-1]), label=label)
    ax_m.set_title("固定 δ=1、M=4：零点越多，幅度越低")
    ax_m.set_xlabel("ω / π")
    ax_m.set_ylabel("|H| / dB")
    ax_m.set_ylim(-80, 3)
    ax_m.grid(alpha=0.3)
    ax_m.legend(fontsize=8)
    ax_g.set_title("固定 δ=1：群延迟误差曲线完全重合")
    ax_g.set_xlabel("ω / π")
    ax_g.set_ylabel("τ$_g$ − τ0 / 样本")
    ax_g.grid(alpha=0.3)
    ax_g.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "zero_pole_scan_L.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def plot_compare():
    fig, (ax_m, ax_g) = plt.subplots(1, 2, figsize=(11, 4.2))
    N, M, L, tau0 = 3, 3, 3, 3.0
    bz, az = zhang_coef(N, M, tau0, L)
    bc, ac, delta = cascade_coef(M, L, tau0)
    for b, a, name in ((bz, az, "Zhang"), (bc, ac, "级联 Thiran")):
        ax_m.plot(WX, 20.0 * np.log10(np.maximum(np.abs(response(b, a, W)), 1e-12)), label=name)
        ax_g.plot(WX[:-1], delay_error(b, a, tau0, W[:-1]), label=name)
    ax_m.set_title("N=M=3, L=3, τ0=3：幅度")
    ax_m.set_xlabel("ω / π")
    ax_m.set_ylabel("|H| / dB")
    ax_m.set_ylim(-80, 3)
    ax_m.grid(alpha=0.3)
    ax_m.legend(fontsize=8)
    ax_g.set_title("同一分子：Zhang 平坦幅度，级联平坦延迟")
    ax_g.set_xlabel("ω / π")
    ax_g.set_ylabel("τ$_g$ − τ0 / 样本")
    ax_g.set_yscale("symlog", linthresh=1e-6)
    ax_g.grid(alpha=0.3)
    ax_g.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "zero_pole_compare.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    OUT_DIR.mkdir(exist_ok=True)
    check_cascade()
    check_zhang()
    plot_scan_L()
    plot_compare()
    print(f"图已写入 {OUT_DIR}")
