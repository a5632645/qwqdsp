# -*- coding: utf-8 -*-
"""
allpole_bandpass.py
===================

全极点带通谐振器：在中心频率 ω0 附近最大平坦群延迟。

H(z) = 1 / A(z)，A 由 R=N/2 个共轭极点对组成：

    A(z) = prod_k (1 - 2 r_k cos(theta_k) z^-1 + r_k^2 z^-2)

未知量 r_k, theta_k 共 N 个。条件：

    tau_g(ω0) = tau0
    tau_g^(m)(ω0) = 0, m = 1..N-2
    d|H|/dω |_{ω0} = 0              （中心是幅度峰）

共 N 个方程，数值求解。若不加峰值条件，可用满 N-1 个延迟导数，
但幅度峰会偏离 ω0；脚本默认加峰值条件。

注意：这是全极点谐振器，不是真正的带通；H(±1) 不为零，阻带只有有限衰减。
"""

import math
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei"]
matplotlib.rcParams["axes.unicode_minus"] = False

OUT_DIR = Path(__file__).parent / "output"
W = np.linspace(0.0, np.pi, 2001)
WX = W / np.pi


# ------------------------------------------------------------
# 系数、频域量、Taylor 展开
# ------------------------------------------------------------

def pairs_to_coef(x):
    a = np.array([1.0])
    for k in range(0, len(x), 2):
        r, th = x[k], x[k + 1]
        a = np.convolve(a, np.array([1.0, -2.0 * r * math.cos(th), r * r]))
    return a


def response(a, w):
    z = np.exp(-1j * np.atleast_1d(w))
    k = np.arange(len(a))
    return 1.0 / np.sum(a[None, :] * z[:, None] ** k[None, :], axis=1)


def group_delay(a, w):
    z = np.exp(-1j * np.atleast_1d(w))
    k = np.arange(len(a))
    terms = a[None, :] * z[:, None] ** k[None, :]
    return -np.real((k[None, :] * terms).sum(axis=1) / terms.sum(axis=1))


def log_series(a, w0, order):
    """A(e^{j(w0+δ)})/A(e^{jw0}) 的 log 级数 L[0..order]。"""
    n = np.arange(len(a))
    e = np.exp(-1j * n * w0)
    A = [np.sum(a * (-1j * n) ** m * e) / math.factorial(m) for m in range(order + 1)]
    b = [1.0] + [A[m] / A[0] for m in range(1, order + 1)]
    L = [0.0] * (order + 1)
    for m in range(1, order + 1):
        s = sum(j * L[j] * b[m - j] for j in range(1, m))
        L[m] = b[m] - s / m
    return np.array(L)


def tau_derivs(a, w0, order):
    L = log_series(a, w0, order)
    return np.array([math.factorial(k + 1) * L[k + 1].imag for k in range(order)])


def max_pole(a):
    if np.allclose(a[1:], 0.0):
        return 0.0
    return float(np.max(np.abs(np.roots(a))))


# ------------------------------------------------------------
# 设计器
# ------------------------------------------------------------

def design(N, w0, tau0, peak=True, trials=6, seed=0):
    if N % 2:
        raise ValueError("N 需为偶数")
    R = N // 2
    r0 = tau0 / (tau0 + R)
    base = np.zeros(N)
    for k in range(R):
        base[2 * k] = min(0.95, r0 * (0.7 + 0.5 * (k + 1) / R))
        base[2 * k + 1] = w0 + (k - (R - 1) / 2.0) * 0.08

    def residual(x):
        a = pairs_to_coef(x)
        t = tau_derivs(a, w0, N)
        vals = [(t[0] - tau0) / max(1.0, tau0)]
        stop = N - 1 if peak else N
        vals += [t[k] / (max(1.0, tau0) * w0 ** k) for k in range(1, stop)]
        if peak:
            vals.append(log_series(a, w0, 1)[1].real)
        return np.array(vals)

    bounds = ([0.001, w0 * 0.3] * R, [0.999, w0 * 1.7] * R)
    rng = np.random.default_rng(seed)
    best = None
    for trial in range(trials):
        g = base.copy()
        if trial:
            g[0::2] = np.clip(g[0::2] * (0.75 + 0.05 * trial), 0.01, 0.99)
            g[1::2] = np.clip(g[1::2] + 0.04 * rng.standard_normal(R), 0.02, np.pi - 0.02)
        r = least_squares(residual, g, bounds=bounds, xtol=1e-14, ftol=1e-14,
                          gtol=1e-14, max_nfev=30000)
        a = pairs_to_coef(r.x)
        if max_pole(a) < 1.0 and (best is None or r.cost < best[0].cost):
            best = (r, a)
    if best is None:
        raise RuntimeError("未找到稳定解")
    r, a = best
    return r, a


# ------------------------------------------------------------
# 自检与表
# ------------------------------------------------------------

def check():
    print("--- 全极点谐振器：峰值约束，中心 0.4π ---")
    print("  N  tau0  max|p|   residual     |H(w0)|   |H(0)|   |H(pi)|  tau(w0)")
    cases = [(4, 2.0), (4, 3.0), (4, 4.0), (6, 2.0), (6, 3.0), (6, 4.0)]
    for N, tau0 in cases:
        r, a = design(N, 0.4 * np.pi, tau0, peak=True, trials=4)
        h0 = abs(response(a, 0.4 * np.pi)[0])
        tau = group_delay(a, 0.4 * np.pi)[0]
        print(f"  {N:>2} {tau0:5.1f}  {max_pole(a):6.4f}  {r.cost:10.2e}  {h0:8.4f}  "
              f"{abs(response(a, 0.0)[0]):7.4f}  {abs(response(a, np.pi)[0]):8.4f}  {tau:8.4f}")

    print("--- 平坦性：N=4/6/8，tau0=3，中心 0.4π ---")
    print("  N  d=0.01π 的 tau 误差   d=0.02π")
    for N in (4, 6, 8):
        r, a = design(N, 0.4 * np.pi, 3.0, peak=True, trials=5)
        e1 = group_delay(a, 0.4 * np.pi + 0.01 * np.pi)[0] - 3.0
        e2 = group_delay(a, 0.4 * np.pi + 0.02 * np.pi)[0] - 3.0
        print(f"  {N:>2}  {e1:+.3e}              {e2:+.3e}")


# ------------------------------------------------------------
# 图
# ------------------------------------------------------------

def plot_scan_tau():
    fig, (ax_m, ax_g) = plt.subplots(1, 2, figsize=(11, 4.2))
    N, w0 = 4, 0.4 * np.pi
    for tau0 in (2.0, 3.0, 4.0):
        r, a = design(N, w0, tau0, peak=True, trials=4)
        h = np.abs(response(a, W))
        ax_m.plot(WX, 20.0 * np.log10(h / h[np.argmin(np.abs(W - w0))]), label=f"τ0={tau0:g}")
        ax_g.plot(WX, group_delay(a, W), label=f"τ0={tau0:g}")
    ax_m.axvline(w0 / np.pi, color="tab:red", lw=0.8, ls="--")
    ax_m.set_title("N=4 全极点谐振器：幅度（中心归一化）")
    ax_m.set_xlabel("ω / π")
    ax_m.set_ylabel("|H| / dB")
    ax_m.set_ylim(-40, 3)
    ax_m.grid(alpha=0.3)
    ax_m.legend(fontsize=8)
    ax_g.axvline(w0 / np.pi, color="tab:red", lw=0.8, ls="--")
    ax_g.set_title("群延迟：中心处平顶")
    ax_g.set_xlabel("ω / π")
    ax_g.set_ylabel("τ$_g$ / 样本")
    ax_g.grid(alpha=0.3)
    ax_g.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "allpole_resonator_scan_tau.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def plot_scan_order():
    fig, (ax_m, ax_g) = plt.subplots(1, 2, figsize=(11, 4.2))
    w0, tau0 = 0.4 * np.pi, 3.0
    for N in (4, 6, 8):
        r, a = design(N, w0, tau0, peak=True, trials=5)
        h = np.abs(response(a, W))
        ax_m.plot(WX, 20.0 * np.log10(h / h[np.argmin(np.abs(W - w0))]), label=f"N={N}")
        ax_g.plot(WX, group_delay(a, W) - tau0, label=f"N={N}")
    ax_m.axvline(w0 / np.pi, color="tab:red", lw=0.8, ls="--")
    ax_m.set_title("τ0=3：阶数越高，滚降越陡")
    ax_m.set_xlabel("ω / π")
    ax_m.set_ylabel("|H| / dB")
    ax_m.set_ylim(-40, 3)
    ax_m.grid(alpha=0.3)
    ax_m.legend(fontsize=8)
    ax_g.axvline(w0 / np.pi, color="tab:red", lw=0.8, ls="--")
    ax_g.set_title("群延迟误差：中心处平坦阶数 N-2")
    ax_g.set_xlabel("ω / π")
    ax_g.set_ylabel("τ$_g$ − τ0 / 样本")
    ax_g.set_ylim(-3, 1)
    ax_g.grid(alpha=0.3)
    ax_g.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "allpole_resonator_scan_order.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    OUT_DIR.mkdir(exist_ok=True)
    check()
    plot_scan_tau()
    plot_scan_order()
    print(f"图已写入 {OUT_DIR}")
