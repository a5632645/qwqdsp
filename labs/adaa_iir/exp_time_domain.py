"""论文 §IV 时域图复现：Fig. 1（锯齿波三方法叠图）、Fig. 5（波形）、Fig. 6（瞬态）。

用法: python qwqdsp/labs/adaa_iir/exp_time_domain.py
输出: qwqdsp/labs/adaa_iir/output/fig1_saw_time.png
      qwqdsp/labs/adaa_iir/output/fig5_wavetable.png
      qwqdsp/labs/adaa_iir/output/fig6_transient.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np

import aaiir as A

FS = 44100.0
OUT = Path(__file__).parent / "output"


def fig1() -> None:
    """Fig. 1：5 kHz 锯齿波，trivial（虚线）/ DPW-N（实线）/ AA-FIR（圆点）。"""
    f0, n = 5000.0, 20
    tab = A.saw_table()
    x = A.phase_ramp(n, f0, FS)
    trivial = A.trivial(tab, n, f0, FS)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
    for ax, order in zip(axes, (2, 3)):
        dpw = A.dpw(tab, n, f0, FS, order)
        fir = A.aa_fir(tab, n, f0, FS, order - 1)
        ax.plot(np.arange(n), trivial, "k--", lw=1, label="trivial")
        ax.plot(np.arange(n), dpw, "-", color="tab:red", lw=1.5, label=f"DPW-{order}")
        ax.plot(np.arange(n), fir, "o", mfc="none", color="tab:blue", ms=6,
                label=f"AA-FIR-{order - 1}")
        ax.set_title(f"阶数 {order}（DPW-{order} ≡ AA-FIR-{order - 1}）")
        ax.set_xlabel("采样 n")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        print(f"Fig.1 阶数 {order}: max|DPW-AA-FIR| = "
              f"{np.max(np.abs(dpw - fir)):.3e}, max|AA-FIR-trivial| = "
              f"{np.max(np.abs(fir - trivial)):.3e}")
    axes[0].set_ylabel("幅值")
    fig.suptitle(f"Fig. 1 复现：5 kHz 锯齿波（fs = {FS:.0f} Hz）")
    fig.tight_layout()
    fig.savefig(OUT / "fig1_saw_time.png", dpi=150)
    plt.close(fig)


def fig5() -> None:
    """Fig. 5：Escalation II 第 3 张表（2048 点）与它的分段线性表示。"""
    tab = A.escalation_ii_w3_table()
    fig, ax = plt.subplots(figsize=(6.5, 3))
    ax.plot(tab.wt, lw=1, label="Escalation II #3（2048 点）")
    xs = np.linspace(0.0, tab.T, 4000, endpoint=False)
    ax.plot(xs * len(tab.wt), A.pwl_value(xs, tab), "--", lw=1,
            label=f"分段线性表示（k = {tab.k} 段）")
    ax.set_xlabel("表内采样")
    ax.set_ylabel("幅值")
    ax.set_title("Fig. 5 复现：Escalation II 波形（Kangaroo 缺失，未复现）")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_wavetable.png", dpi=150)
    plt.close(fig)


def fig6() -> None:
    """Fig. 6：Escalation 波形在跳变处的瞬态，四种方法的抗混叠滤波效果。"""
    f0, n = 1000.0, 220
    tab = A.escalation_ii_w3_table()
    x = A.phase_ramp(n, f0, FS)
    sigs = {
        "trivial": A.trivial(tab, n, f0, FS),
        "OVS-2": A.oversampled(tab, n, f0, FS, 2),
        "AA-FIR-2": A.aa_fir(tab, n, f0, FS, 2),
        "AA-IIR-1": A.aa_iir_osc(x, A.aa_iir_1(), tab),
        "AA-IIR-2": A.aa_iir_osc(x, A.aa_iir_2(), tab),
    }
    lo, hi = 40, 120        # 跳变出现在第一个周期末尾附近
    fig, ax = plt.subplots(figsize=(9, 3.4))
    style = {"trivial": ("k", "-", 1.0), "OVS-2": ("tab:green", "-", 1.0),
             "AA-FIR-2": ("tab:orange", "-", 1.2), "AA-IIR-1": ("0.45", "-", 1.6),
             "AA-IIR-2": ("tab:blue", "-", 1.6)}
    for name, y in sigs.items():
        c, ls, lw = style[name]
        ax.plot(np.arange(lo, hi), y[lo:hi], ls, color=c, lw=lw, label=name)
    ax.set_xlabel("采样 n")
    ax.set_ylabel("幅值")
    ax.set_title("Fig. 6 复现：Escalation 波形跳变处的瞬态（f0 = 1 kHz）")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig6_transient.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    fig1()
    fig5()
    fig6()
    print(f"已写出 {OUT} 下 fig1/fig5/fig6 三张图")
