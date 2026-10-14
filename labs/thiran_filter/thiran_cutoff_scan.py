# -*- coding: utf-8 -*-
"""
thiran_cutoff_scan.py
=====================

用 thiran_cutoff_design.py 的设计器（固定阶数、指定 -3dB 截止频率）做两张扫描图：

    thiran_cutoff_scan_wc.png   固定阶数 n=6，扫截止频率 wc
    thiran_cutoff_scan_n.png    固定截止频率 wc=0.2π，扫阶数 n

每张图两栏：幅度响应（按直流归一化，图例里给出实际直流增益）与群延迟（样本）。
幅度栏在每个目标 wc 处画一个点（实测 -3dB）与一条虚线，直观核对设计是否落在目标上。
"""

import math
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from thiran_cutoff_design import measure_cutoff, solve_beta, thiran_allpole_coef

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei"]
matplotlib.rcParams["axes.unicode_minus"] = False

OUT_DIR = Path(__file__).parent / "output"
W = np.linspace(0.0, np.pi, 4001)
WX = W / np.pi


def response(a, w):
    k = np.arange(len(a))
    return np.sum(a[None, :] * np.exp(-1j * np.outer(w, k)), axis=1)


def group_delay(a, w):
    """H = 1/A 的群延迟（样本）：tau = -Re[Σ k a_k e^-jkw / Σ a_k e^-jkw]"""
    k = np.arange(len(a))
    terms = a[None, :] * np.exp(-1j * np.outer(w, k))
    return -np.real((k[None, :] * terms).sum(axis=1) / terms.sum(axis=1))


def draw(n, beta, wc, label, color, ax_m, ax_g, symlog=False):
    a = thiran_allpole_coef(beta, n)
    mag_db = -20.0 * np.log10(np.abs(response(a, W)))  # H = 1/A
    dc_db = mag_db[0]
    ax_m.plot(WX, mag_db - dc_db, color=color, lw=1.6, label=f"{label}（直流增益 {dc_db:+.1f} dB）")
    ax_m.plot([wc / np.pi], [-3.0103], "o", color=color, ms=5, zorder=5)
    ax_g.plot(WX, group_delay(a, W), color=color, lw=1.6, label=label)
    if symlog:
        ax_g.set_yscale("symlog", linthresh=0.1)
    return dc_db


def scan_fixed_order(n=6, cutoffs=(0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.5, 0.7, 0.9)):
    fig, (ax_m, ax_g) = plt.subplots(2, 1, figsize=(9.5, 8.0), sharex=True)
    colors = plt.get_cmap("viridis")(np.linspace(0.03, 0.92, len(cutoffs)))
    print(f"--- 固定 n={n}，扫截止频率 ---")
    for wc_pi, color in zip(cutoffs, colors):
        wc = wc_pi * math.pi
        beta = solve_beta(n, wc)
        got = measure_cutoff(beta, n)
        draw(n, beta, wc, f"$\\omega_c/\\pi$={wc_pi:g}", color, ax_m, ax_g)
        print(f"  wc/pi={wc_pi:g}: delta={beta/2:8.4f}  实测 -3dB={got/math.pi:.12f}（目标 {wc_pi}）"
              f"  误差 {abs(got/math.pi - wc_pi):.2e}")

    ax_m.axhline(0.0, color="0.6", lw=0.8, ls=":")
    ax_m.axhline(-3.0103, color="tab:red", lw=0.8, ls="--")
    ax_m.set_ylim(-60, 2)
    ax_m.set_title(f"数字 Thiran 低通（指定截止频率，n={n} 固定）：幅度（直流归一化）")
    ax_m.set_ylabel("|H| / dB")
    ax_m.grid(alpha=0.3)
    ax_m.legend(fontsize=8, ncol=2)
    ax_g.axhline(0.0, color="0.6", lw=0.8, ls=":")
    ax_g.set_ylim(-3.0, 16.5)
    ax_g.set_title("群延迟（直流处 = δ = β/2，由截止频率决定）")
    ax_g.set_xlabel("ω / π")
    ax_g.set_ylabel("τ$_g$ / 样本")
    ax_g.grid(alpha=0.3)
    ax_g.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "thiran_cutoff_scan_wc.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def scan_fixed_cutoff(wc_pi=0.2, orders=(1, 2, 3, 4, 5, 6, 8, 10, 12)):
    fig, (ax_m, ax_g) = plt.subplots(2, 1, figsize=(9.5, 8.0), sharex=True)
    colors = plt.get_cmap("plasma")(np.linspace(0.03, 0.92, len(orders)))
    wc = wc_pi * math.pi
    print(f"--- 固定 wc={wc_pi}π，扫阶数 ---")
    ax_m.axvline(wc_pi, color="tab:red", lw=0.8, ls="--")
    for n, color in zip(orders, colors):
        beta = solve_beta(n, wc)
        got = measure_cutoff(beta, n)
        draw(n, beta, wc, f"n={n}", color, ax_m, ax_g)
        print(f"  n={n:2d}: delta={beta/2:8.4f}  实测 -3dB/pi={got/math.pi:.12f}（目标 {wc_pi}）"
              f"  误差 {abs(got/math.pi - wc_pi):.2e}")

    ax_m.axhline(0.0, color="0.6", lw=0.8, ls=":")
    ax_m.axhline(-3.0103, color="tab:red", lw=0.8, ls="--")
    ax_m.set_ylim(-80, 3)
    ax_m.set_title(f"数字 Thiran 低通（指定截止频率 $\\omega_c/\\pi$={wc_pi:g}，扫阶数）：幅度（直流归一化）")
    ax_m.set_ylabel("|H| / dB")
    ax_m.grid(alpha=0.3)
    ax_m.legend(fontsize=8, ncol=2)
    ax_g.set_title("群延迟（直流处 = δ = β/2；阶数越高峰越高）")
    ax_g.set_xlabel("ω / π")
    ax_g.set_ylabel("τ$_g$ / 样本")
    ax_g.grid(alpha=0.3)
    ax_g.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "thiran_cutoff_scan_n.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    OUT_DIR.mkdir(exist_ok=True)
    scan_fixed_order()
    scan_fixed_cutoff()
    print(f"\n图已写入 {OUT_DIR}")
