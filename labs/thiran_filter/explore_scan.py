# -*- coding: utf-8 -*-
"""
explore_scan.py
===============

纯 Python 扫描 Thiran 两类滤波器（全通 AP / 全极点 LP）的群延迟与幅度表现。
公式直接照抄文献，与本库的 C++ 实现无关：

    AP: a_k = (-1)^k C(n,k) prod_i (d+i)/(d+k+i)，H(z) = z^-n A(1/z)/A(z)，直流群延迟 n+d
    LP: a_k = (-1)^k C(n,k) prod_i (2d+i)/(2d+k+i)，H(z) = 1/A(z)，直流群延迟 d

目标群延迟不做 [-0.5, 0.5] 的限制，扫描范围一直开到远大于阶数，
用极点模长判断稳定域（AP: D > n-1；LP: d > -0.5）。

四张图（输出在脚本同目录的 output/ 下）：
    ap_scan_delay.png  定阶数 n，扫描目标群延迟 D
    ap_scan_order.png  定目标群延迟 D，扫描阶数 n
    lp_scan_delta.png  定阶数 n，扫描 LP 的目标群延迟 d（含幅度响应）
    lp_scan_order.png  定目标群延迟 d，扫描阶数 n（含幅度响应）
"""

import math
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei"]
matplotlib.rcParams["axes.unicode_minus"] = False

OUT_DIR = Path(__file__).parent / "output"
W = np.linspace(0.0, np.pi, 2001)
WX = W / np.pi


# ------------------------------------------------------------
# 系数与频域量
# ------------------------------------------------------------

def _binom(n, k):
    return np.array([math.comb(n, int(kk)) for kk in k], dtype=float)


def ap_coef(n, d):
    """全通形式的分母系数 a[0..n]，直流群延迟为 n+d。"""
    k = np.arange(1, n + 1)
    i = np.arange(0, n + 1)
    prod = np.prod((d + i)[None, :] / (d + k[:, None] + i[None, :]), axis=1)
    return np.concatenate(([1.0], (-1.0) ** k * _binom(n, k) * prod))


def lp_coef(n, d):
    """全极点形式 1/A(z) 的分母系数 a[0..n]，直流群延迟为 d（乘积项用 2d）。"""
    return ap_coef(n, 2.0 * d)


def _tau_of(c, w):
    """多项式 sum_k c_k e^{-jkw} 对群延迟的贡献。"""
    k = np.arange(len(c))
    terms = c[None, :] * np.exp(-1j * np.outer(w, k))
    return np.real((k[None, :] * terms).sum(axis=1) / terms.sum(axis=1))


def group_delay(b, a, w):
    return _tau_of(b, w) - _tau_of(a, w)


def magnitude(b, a, w):
    k = np.arange(len(a))
    den = (a[None, :] * np.exp(-1j * np.outer(w, k))).sum(axis=1)
    if len(b) == 1:
        return np.full_like(w, abs(b[0]), dtype=float) / np.abs(den)
    num = (b[None, :] * np.exp(-1j * np.outer(w, k))).sum(axis=1)
    return np.abs(num / den)


def max_pole_radius(a):
    """A(z) = sum a_k z^-k 的极点模长最大值（A 退化为常数时无极点）。"""
    if np.allclose(a[1:], 0.0):
        return np.nan
    return float(np.max(np.abs(np.roots(a))))


# ------------------------------------------------------------
# 出图
# ------------------------------------------------------------

def ap_scan_delay(n=6):
    deltas = [1.0, 2.0, 3.0, 4.0, 5.25, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 13.0]
    stable = [d for d in deltas if d > n - 1]
    unstable = [d for d in deltas if d <= n - 1]

    fig = plt.figure(figsize=(12.0, 7.6))
    grid = fig.add_gridspec(2, 2, height_ratios=[2.4, 1.0], hspace=0.32, wspace=0.22)
    ax_g = fig.add_subplot(grid[0, 0])
    ax_e = fig.add_subplot(grid[0, 1])
    ax_p = fig.add_subplot(grid[1, :])

    cmap = plt.get_cmap("viridis")
    colors = {d: cmap(i / max(1, len(stable) - 1)) for i, d in enumerate(stable)}

    for d in unstable:
        a = ap_coef(n, d)
        tau = group_delay(a[::-1], a, W)
        ax_g.plot(WX, tau, color="0.72", ls="--", lw=1.0)
        ax_e.plot(WX, tau - d, color="0.72", ls="--", lw=1.0)

    for d in stable:
        a = ap_coef(n, d)
        tau = group_delay(a[::-1], a, W)
        lw = 2.6 if d == float(n) else 1.4
        ax_g.plot(WX, tau, color=colors[d], lw=lw, label=f"{d:g}")
        ax_e.plot(WX, tau - d, color=colors[d], lw=lw)

    ax_g.axhline(n, color="k", lw=0.8, ls=":")
    ax_g.set_title(f"全通：定阶数 n={n}，扫描目标群延迟 D（曲线颜色 = D）")
    ax_g.set_ylabel("群延迟 τ$_g$ / 样本")
    ax_g.set_xlabel("ω / π")
    ax_g.grid(alpha=0.3)
    ax_g.legend(ncol=2, fontsize=8, title="稳定 (D > n−1)，虚线为不稳定", title_fontsize=8)

    ax_e.axhline(0.0, color="k", lw=0.8, ls=":")
    ax_e.set_title("相对目标的偏差 τ$_g$ − D")
    ax_e.set_ylabel("偏差 / 样本")
    ax_e.set_xlabel("ω / π")
    ax_e.grid(alpha=0.3)

    dd = np.linspace(0.2, 14.0, 400)
    radii = np.array([max_pole_radius(ap_coef(n, x)) for x in dd])
    ax_p.plot(dd, radii, color="tab:blue", lw=1.6)
    ax_p.axhline(1.0, color="tab:red", lw=1.0, ls="--", label="|p| = 1")
    ax_p.axvline(n - 1, color="0.5", lw=1.0, ls=":", label=f"D = n−1 = {n-1}")
    ax_p.axvspan(n - 0.5, n + 0.5, color="tab:green", alpha=0.10, label="|D − n| ≤ 0.5")
    for d in deltas:
        r = max_pole_radius(ap_coef(n, d))
        ax_p.plot([d], [r], "o", ms=4, color=colors.get(d, "0.6"))
    ax_p.set_xlabel("目标群延迟 D")
    ax_p.set_ylabel("max|p|")
    ax_p.set_title("极点模长扫描（D ≤ n−1 不稳定，D = n−1 时公式奇异）")
    ax_p.grid(alpha=0.3)
    ax_p.legend(fontsize=8)

    fig.savefig(OUT_DIR / "ap_scan_delay.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return deltas


def ap_scan_order(target=6.5, orders=(1, 2, 3, 4, 5, 6, 7, 8, 9)):
    stable = [n for n in orders if target > n - 1]
    unstable = [n for n in orders if target <= n - 1]

    fig = plt.figure(figsize=(12.0, 7.6))
    grid = fig.add_gridspec(2, 2, height_ratios=[2.4, 1.0], hspace=0.32, wspace=0.22)
    ax_g = fig.add_subplot(grid[0, 0])
    ax_e = fig.add_subplot(grid[0, 1])
    ax_p = fig.add_subplot(grid[1, :])

    cmap = plt.get_cmap("plasma")
    colors = {n: cmap(i / max(1, len(stable) - 1)) for i, n in enumerate(stable)}

    for n in unstable:
        a = ap_coef(n, target - n)
        tau = group_delay(a[::-1], a, W)
        ax_g.plot(WX, tau, color="0.72", ls="--", lw=1.0)
        ax_e.plot(WX, tau - target, color="0.72", ls="--", lw=1.0)

    for n in stable:
        a = ap_coef(n, target - n)
        tau = group_delay(a[::-1], a, W)
        ax_g.plot(WX, tau, color=colors[n], lw=1.6, label=f"n={n}")
        ax_e.plot(WX, tau - target, color=colors[n], lw=1.6)

    ax_g.axhline(target, color="k", lw=1.0, ls="--")
    ax_g.set_title(f"全通：定目标群延迟 D={target}，扫描阶数 n（虚线 = 理想延迟）")
    ax_g.set_ylabel("群延迟 τ$_g$ / 样本")
    ax_g.set_xlabel("ω / π")
    ax_g.grid(alpha=0.3)
    ax_g.legend(ncol=2, fontsize=8, title="稳定 (D > n−1)，虚线为不稳定", title_fontsize=8)

    ax_e.axhline(0.0, color="k", lw=0.8, ls=":")
    ax_e.set_title("相对目标的偏差 τ$_g$ − D")
    ax_e.set_ylabel("偏差 / 样本")
    ax_e.set_xlabel("ω / π")
    ax_e.grid(alpha=0.3)

    radii = [max_pole_radius(ap_coef(n, target - n)) for n in orders]
    ax_p.plot(list(orders), radii, "-o", color="tab:blue", ms=5, label="max|p|")
    ax_p.axhline(1.0, color="tab:red", lw=1.0, ls="--", label="|p| = 1")
    ax_p.axvline(target + 1, color="0.5", lw=1.0, ls=":", label=f"D+1 = {target + 1:g}")
    ax_p.axvspan(target - 0.5, target + 0.5, color="tab:green", alpha=0.10, label="|D − n| ≤ 0.5")
    ax_p.set_xlabel("阶数 n")
    ax_p.set_ylabel("max|p|")
    ax_p.set_title("极点模长扫描（n ≥ D+1 时不稳定）")
    ax_p.grid(alpha=0.3)
    ax_p.legend(fontsize=8)

    fig.savefig(OUT_DIR / "ap_scan_order.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return orders


def lp_scan_delta(n=6, deltas=(0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0)):
    fig = plt.figure(figsize=(12.0, 7.6))
    grid = fig.add_gridspec(2, 2, height_ratios=[2.4, 1.0], hspace=0.32, wspace=0.22)
    ax_m = fig.add_subplot(grid[0, 0])
    ax_g = fig.add_subplot(grid[0, 1])
    ax_p = fig.add_subplot(grid[1, :])

    cmap = plt.get_cmap("viridis")
    colors = {d: cmap(i / max(1, len(deltas) - 1)) for i, d in enumerate(deltas)}

    for d in deltas:
        a = lp_coef(n, d)
        mag_db = 20.0 * np.log10(magnitude(np.array([1.0]), a, W))
        dc_db = mag_db[0]
        ax_m.plot(WX, mag_db - dc_db, color=colors[d], lw=1.6, label=f"d={d:g} ({dc_db:+.1f} dB)")
        ax_g.plot(WX, group_delay(np.array([1.0]), a, W), color=colors[d], lw=1.6)

    ax_m.axhline(0.0, color="k", lw=0.8, ls=":")
    ax_m.set_title(f"全极点：定阶数 n={n}，扫描目标群延迟 d（颜色 = d）")
    ax_m.set_ylabel("|H| / dB（直流归一化）")
    ax_m.set_xlabel("ω / π")
    ax_m.grid(alpha=0.3)
    ax_m.legend(ncol=2, fontsize=8, title="直流群延迟 = d，括号内为实际直流增益", title_fontsize=8)

    ax_g.axhline(0.0, color="k", lw=0.8, ls=":")
    ax_g.set_title("群延迟 τ$_g$（直流处 = d，带内变号）")
    ax_g.set_ylabel("群延迟 / 样本")
    ax_g.set_xlabel("ω / π")
    ax_g.grid(alpha=0.3)

    dd = np.linspace(-0.9, 6.0, 400)
    radii = np.array([max_pole_radius(lp_coef(n, x)) for x in dd])
    ax_p.plot(dd, radii, color="tab:blue", lw=1.6)
    ax_p.axhline(1.0, color="tab:red", lw=1.0, ls="--", label="|p| = 1")
    ax_p.axvline(-0.5, color="0.5", lw=1.0, ls=":", label="d = −0.5（奇异/边界）")
    ax_p.axvline(0.0, color="tab:green", lw=1.0, ls=":", label="d = 0（无极点，H≡1）")
    for d in deltas:
        ax_p.plot([d], [max_pole_radius(lp_coef(n, d))], "o", ms=4, color=colors[d])
    ax_p.set_xlabel("目标群延迟 d")
    ax_p.set_ylabel("max|p|")
    ax_p.set_title("极点模长扫描（d > −0.5 稳定）")
    ax_p.grid(alpha=0.3)
    ax_p.legend(fontsize=8)

    fig.savefig(OUT_DIR / "lp_scan_delta.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return deltas


def lp_scan_order(target=2.5, orders=(1, 2, 3, 4, 5, 6, 8, 10)):
    fig = plt.figure(figsize=(12.0, 7.6))
    grid = fig.add_gridspec(2, 2, height_ratios=[2.4, 1.0], hspace=0.32, wspace=0.22)
    ax_m = fig.add_subplot(grid[0, 0])
    ax_g = fig.add_subplot(grid[0, 1])
    ax_p = fig.add_subplot(grid[1, :])

    cmap = plt.get_cmap("plasma")
    colors = {n: cmap(i / max(1, len(orders) - 1)) for i, n in enumerate(orders)}

    for n in orders:
        a = lp_coef(n, target)
        mag_db = 20.0 * np.log10(magnitude(np.array([1.0]), a, W))
        dc_db = mag_db[0]
        ax_m.plot(WX, mag_db - dc_db, color=colors[n], lw=1.6, label=f"n={n} ({dc_db:+.1f} dB)")
        ax_g.plot(WX, group_delay(np.array([1.0]), a, W), color=colors[n], lw=1.6)

    ax_m.axhline(0.0, color="k", lw=0.8, ls=":")
    ax_m.set_title(f"全极点：定目标群延迟 d={target}，扫描阶数 n（颜色 = n）")
    ax_m.set_ylabel("|H| / dB（直流归一化）")
    ax_m.set_xlabel("ω / π")
    ax_m.grid(alpha=0.3)
    ax_m.legend(ncol=2, fontsize=8, title=f"直流群延迟 d = {target:g}，括号内为实际直流增益", title_fontsize=8)

    ax_g.axhline(0.0, color="k", lw=0.8, ls=":")
    ax_g.set_title("群延迟 τ$_g$（所有阶数在直流处都等于 d）")
    ax_g.set_ylabel("群延迟 / 样本")
    ax_g.set_xlabel("ω / π")
    ax_g.grid(alpha=0.3)

    ax_p.plot(list(orders), [max_pole_radius(lp_coef(n, target)) for n in orders], "-o", color="tab:blue", ms=5)
    ax_p.axhline(1.0, color="tab:red", lw=1.0, ls="--", label="|p| = 1")
    ax_p.set_xlabel("阶数 n")
    ax_p.set_ylabel("max|p|")
    ax_p.set_title("极点模长扫描（d > −0.5 时各阶数都稳定）")
    ax_p.grid(alpha=0.3)
    ax_p.legend(fontsize=8)

    fig.savefig(OUT_DIR / "lp_scan_order.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return orders


# ------------------------------------------------------------
# 数值自检
# ------------------------------------------------------------

def report():
    print("--- 全通：直流群延迟 vs 目标 ---")
    for n in (2, 6):
        for d in (0.5, 3.0, 6.0, 12.0):
            a = ap_coef(n, d)
            tau = group_delay(a[::-1], a, np.array([0.0]))[0]
            mean = np.trapezoid(group_delay(a[::-1], a, W), W) / np.pi
            print(f"  n={n} D={n+d:5.2f}  tau_g(0)={tau:10.6f}  带回均值={mean:8.5f} (=n)  max|p|={max_pole_radius(a):.4f}")

    print("--- 全极点：直流群延迟 vs 目标 ---")
    for n in (2, 6):
        for d in (0.1, 0.5, 2.0, 5.0):
            a = lp_coef(n, d)
            tau = group_delay(np.array([1.0]), a, np.array([0.0]))[0]
            print(f"  n={n} d={d:4.2f}  tau_g(0)={tau:10.6f}  |H(0)|dB={20 * np.log10(abs(1 / a.sum())):8.3f}"
                  f"  max|p|={max_pole_radius(a):.4f}")

    print("--- 全通 n=6：目标离阶数越远，带内偏差越大 ---")
    print("   D-n   |dev|@0.5pi  |dev|@0.75pi  |dev|@pi   max|dev|  mean|dev|  max|p|")
    for offset in (-0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 7.0):
        n = 6
        d = offset
        a = ap_coef(n, d)
        tau = group_delay(a[::-1], a, W)
        dev = tau - (n + d)
        print(f"  {offset:+6.2f}   {abs(dev[np.searchsorted(W, 0.5 * np.pi)]):11.4f}   "
              f"{abs(dev[np.searchsorted(W, 0.75 * np.pi)]):11.4f}   {abs(dev[-1]):9.4f}   "
              f"{np.max(np.abs(dev)):9.4f}   {np.trapezoid(np.abs(dev), W) / np.pi:9.4f}   {max_pole_radius(a):.4f}")


if __name__ == "__main__":
    OUT_DIR.mkdir(exist_ok=True)
    report()
    ap_scan_delay()
    ap_scan_order()
    lp_scan_delta()
    lp_scan_order()
    print(f"\n图已写入 {OUT_DIR}")
