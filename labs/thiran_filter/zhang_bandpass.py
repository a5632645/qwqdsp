# -*- coding: utf-8 -*-
"""
zhang_bandpass.py
=================

Zhang 2008 的带通线性方程组：零极点带通，在中心 ω0 附近平坦幅度 + 平坦群延迟。

H(z) = B(z)/A(z)，B 阶数 N，A 阶数 M，a0=1。要求

    H(e^{jω}) e^{j(τ0 ω + θ0)} = 1 + O((ω-ω0)^K)

并在 z=1 处放 L1 个零点、z=-1 处放 L2 个零点。
实部/虚部展开后共 2K+L1+L2 个方程，令

    2K + L1 + L2 = N + M + 1

即可解出 b0..bN 和 a1..aM。

注意：这里的指数采用标准 z^-1 约定（A(e^{jω}) = Σ a_m e^{-jmω}）。
文献里若把频率响应写成 Σ a_m e^{+jmω}，分母的指数符号要相应翻转。
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


def design(N, M, w0, tau0, K, L1, L2, theta0=0.0):
    if 2 * K + L1 + L2 != N + M + 1:
        raise ValueError("需满足 2K + L1 + L2 = N + M + 1")
    n = np.arange(N + 1)
    m = np.arange(1, M + 1)
    rows = []
    rhs = []
    for i in range(K):
        rn = (n - tau0) ** i * np.exp(1j * ((tau0 - n) * w0 + theta0))
        sm = m ** i * np.exp(-1j * m * w0)
        row = np.zeros(N + 1 + M)
        row[:N + 1] = rn.real
        row[N + 1:] = -sm.real
        rows.append(row)
        rhs.append(1.0 if i == 0 else 0.0)
        row = np.zeros(N + 1 + M)
        row[:N + 1] = rn.imag
        row[N + 1:] = -sm.imag
        rows.append(row)
        rhs.append(0.0)
    for i in range(L1):
        row = np.zeros(N + 1 + M)
        row[:N + 1] = n ** i
        rows.append(row)
        rhs.append(0.0)
    for i in range(L2):
        row = np.zeros(N + 1 + M)
        row[:N + 1] = ((-1) ** n) * (n ** i)
        rows.append(row)
        rhs.append(0.0)
    mat = np.array(rows)
    vec = np.array(rhs)
    x = np.linalg.solve(mat, vec)
    return x[:N + 1], np.concatenate(([1.0], x[N + 1:])), mat, vec, x


def response(b, a, w):
    z = np.exp(-1j * np.atleast_1d(w))
    num = np.sum(b[None, :] * z[:, None] ** np.arange(len(b))[None, :], axis=1)
    den = np.sum(a[None, :] * z[:, None] ** np.arange(len(a))[None, :], axis=1)
    return num / den


def group_delay(b, a, w):
    z = np.exp(-1j * np.atleast_1d(w))

    def tau(c):
        k = np.arange(len(c))
        terms = c[None, :] * z[:, None] ** k[None, :]
        return np.real((k[None, :] * terms).sum(axis=1) / terms.sum(axis=1))

    return tau(b) - tau(a)


def max_pole(a):
    if np.allclose(a[1:], 0.0):
        return 0.0
    return float(np.max(np.abs(np.roots(a))))


def check():
    print("--- Zhang 带通：论文例 N=6 M=4 K=3 L1=2 L2=3 ---")
    N, M, w0, tau0, K, L1, L2 = 6, 4, 0.2 * np.pi, 12.8, 3, 2, 3
    for theta0 in (0.0, 0.25 * np.pi, 0.5 * np.pi):
        b, a, mat, vec, x = design(N, M, w0, tau0, K, L1, L2, theta0)
        print(f"  theta0={theta0/np.pi:.2f}π: residual={np.max(np.abs(mat @ x - vec)):.2e}"
              f"  max|p|={max_pole(a):.4f}  |H(w0)|={abs(response(b, a, w0)[0]):.12f}"
              f"  tau(w0)={group_delay(b, a, w0)[0]:.9f}"
              f"  B(1)={abs(np.sum(b)):.2e}  B(-1)={abs(np.sum(b * (-1) ** np.arange(N + 1))):.2e}")
    print()

    print("--- Zhang 带通：简单例 N=4 M=4 K=3 L1=1 L2=2，中心 0.3π ---")
    N, M, w0, tau0, K, L1, L2 = 4, 4, 0.3 * np.pi, 5.0, 3, 1, 2
    for t0 in (2.0, 4.0, 6.0, 8.0):
        b, a, mat, vec, x = design(N, M, w0, t0, K, L1, L2)
        print(f"  tau0={t0:4.1f}: max|p|={max_pole(a):.4f}  |H(w0)|={abs(response(b, a, w0)[0]):.9f}"
              f"  B(1)={abs(np.sum(b)):.1e}  B(-1)={abs(np.sum(b * (-1) ** np.arange(N + 1))):.1e}")
    print()

    print("--- 平坦性：论文例，theta0=0，K=3 ---")
    b, a, mat, vec, x = design(6, 4, 0.2 * np.pi, 12.8, 3, 2, 3)
    for d in (0.001, 0.002, 0.005):
        err_tau = group_delay(b, a, 0.2 * np.pi + d)[0] - 12.8
        err_mag = abs(response(b, a, 0.2 * np.pi + d)[0]) - 1.0
        print(f"  δ={d:g}: tau error={err_tau:+.3e}  |H| error={err_mag:+.3e}")
    print()


def plot_example():
    N, M, w0, tau0, K, L1, L2 = 10, 6, 0.4 * np.pi, 8.0, 6, 2, 3
    b, a, _, _, _ = design(N, M, w0, tau0, K, L1, L2)
    fig, (ax_m, ax_g) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax_m.plot(WX, 20.0 * np.log10(np.maximum(np.abs(response(b, a, W)), 1e-12)))
    ax_m.axvline(w0 / np.pi, color="tab:red", lw=0.8, ls="--")
    ax_m.set_title("Zhang 带通 N=10 M=6：幅度（中心 0.4π）")
    ax_m.set_xlabel("ω / π")
    ax_m.set_ylabel("|H| / dB")
    ax_m.set_ylim(-100, 5)
    ax_m.grid(alpha=0.3)
    ax_g.plot(WX[1:-1], group_delay(b, a, W[1:-1]))
    ax_g.axvline(w0 / np.pi, color="tab:red", lw=0.8, ls="--")
    ax_g.axhline(tau0, color="0.6", lw=0.8, ls=":")
    ax_g.set_title("群延迟：中心 τ0=8，平坦阶数 K-1=5")
    ax_g.set_xlabel("ω / π")
    ax_g.set_ylabel("τ$_g$ / 样本")
    ax_g.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "zhang_bandpass_example.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def plot_scan_tau():
    N, M, w0, K, L1, L2 = 10, 6, 0.4 * np.pi, 6, 2, 3
    fig, (ax_m, ax_g) = plt.subplots(1, 2, figsize=(11, 4.2))
    for tau0 in (8.0, 10.0, 12.0):
        b, a, _, _, _ = design(N, M, w0, tau0, K, L1, L2)
        ax_m.plot(WX, 20.0 * np.log10(np.maximum(np.abs(response(b, a, W)), 1e-12)), label=f"τ0={tau0:g}")
        ax_g.plot(WX[1:-1], group_delay(b, a, W[1:-1]), label=f"τ0={tau0:g}")
    ax_m.axvline(w0 / np.pi, color="tab:red", lw=0.8, ls="--")
    ax_m.set_title("N=10 M=6：扫描中心延迟")
    ax_m.set_xlabel("ω / π")
    ax_m.set_ylabel("|H| / dB")
    ax_m.set_ylim(-100, 5)
    ax_m.grid(alpha=0.3)
    ax_m.legend(fontsize=8)
    ax_g.axvline(w0 / np.pi, color="tab:red", lw=0.8, ls="--")
    ax_g.set_title("群延迟")
    ax_g.set_xlabel("ω / π")
    ax_g.set_ylabel("τ$_g$ / 样本")
    ax_g.grid(alpha=0.3)
    ax_g.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "zhang_bandpass_scan_tau.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    OUT_DIR.mkdir(exist_ok=True)
    check()
    plot_example()
    plot_scan_tau()
    print(f"图已写入 {OUT_DIR}")
