"""论文 §IV-A/§IV-B：锯齿波 SNR 复现（对应 Fig. 2 与侧表）。

对 88 键钢琴基频逐个生成 2^17 采样，按 ``aaiir.snr_db`` 的固定口径算 SNR
（定义见 README「SNR 定义」）。

**有效性门限**：定义里「谐波 ±4 bin 之外的功率 = 噪声」在低音区会失效——当
``fs/f0`` 接近整数时，折叠分量恰好落在谐波栅格上，被误判成信号（实测 SNR 虚高，
甚至 inf）。用**解析**值 ``aaiir.analytic_saw_snr``（理想锯齿波采样后的折叠分量
功率比）当锚点：trivial 的实测与解析差 > 2 dB 的音视为病态，**逐音剔除**，
剩余音集对全部方法统一使用。剔除的音多为最低的几个八度。

用法: python qwqdsp/labs/adaa_iir/exp_saw_snr.py
输出: qwqdsp/labs/adaa_iir/output/fig2_saw_snr.png
      qwqdsp/labs/adaa_iir/output/saw_snr.csv
      stdout: 平均 SNR 表（本复现 vs 论文）+ 门限说明
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np

import aaiir as A

FS = 44100.0
N = 2 ** 17            # 每音长度 [采样]
WARMUP = 4096          # 丢弃 AA 递推启动瞬态
GUARD = 4              # 谐波 ±guard 个 bin 计入信号
GATE_DB = 2.0          # 病态音判定阈值（trivial 实测 vs 解析）

# 论文 Fig. 2 侧表的平均 SNR [dB]（用于并列比较）
PAPER = {
    "trivial": 8.7,
    "AA-FIR-1/DPW-2": 18.8,
    "AA-FIR-2/DPW-3": 24.5,
    "OVS-2": 17.3,
    "OVS-8": 30.8,
    "AA-IIR-1": 18.8,
    "AA-IIR-2": 64.6,
}


def methods(tab, f0, n):
    """返回 {方法名: 波形}。"""
    x = A.phase_ramp(n, f0, FS)
    return {
        "trivial": A.trivial(tab, n, f0, FS),
        "AA-FIR-1/DPW-2": A.dpw(tab, n, f0, FS, 2),
        "AA-FIR-2/DPW-3": A.dpw(tab, n, f0, FS, 3),
        "OVS-2": A.oversampled(tab, n, f0, FS, 2),
        "OVS-8": A.oversampled(tab, n, f0, FS, 8),
        "AA-IIR-1": A.aa_iir_osc(x, A.aa_iir_1(), tab),
        "AA-IIR-2": A.aa_iir_osc(x, A.aa_iir_2(), tab),
    }


def main() -> int:
    tab = A.saw_table()
    notes = A.piano_frequencies()
    names = list(PAPER.keys())
    analytic = np.array([A.analytic_saw_snr(f0, FS) for f0 in notes])
    snr = {name: [] for name in names}

    for i, f0 in enumerate(notes):
        for name, y in methods(tab, f0, N).items():
            snr[name].append(A.snr_db(y, FS, f0, guard_bins=GUARD, warmup=WARMUP))
        if (i + 1) % 11 == 0:
            print(f"  {i + 1}/88 键完成")

    for name in names:
        snr[name] = np.array(snr[name])

    dev = snr["trivial"] - analytic
    keep = np.isfinite(dev) & (np.abs(dev) <= GATE_DB)
    print(f"\n有效性门限：{int((~keep).sum())}/88 个音被剔除"
          f"（trivial 实测 vs 解析偏离 > {GATE_DB:.0f} dB，折叠分量落在谐波栅格上）：")
    print("  " + ", ".join(f"{f:.1f} Hz" for f in notes[~keep]))

    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    with open(out / "saw_snr.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["f0_hz", "analytic_trivial", "kept"] + names)
        for i, f0 in enumerate(notes):
            w.writerow([f"{f0:.4f}", f"{analytic[i]:.3f}", int(keep[i])]
                       + [f"{snr[n][i]:.3f}" for n in names])

    print(f"\n平均 SNR [dB]（论文 Fig. 2 侧表 vs 本复现，{int(keep.sum())} 个音）")
    print(f"{'方法':<18}{'论文':>8}{'本复现':>10}{'差':>8}")
    for name in names:
        m = float(np.mean(snr[name][keep]))
        print(f"{name:<18}{PAPER[name]:>8.1f}{m:>10.1f}{m - PAPER[name]:>+8.1f}")
    print("论文未给 SNR 计算细节；本表口径已用解析锚点验证（见 README）。")

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.semilogx(notes, analytic, "k:", lw=1.5, label="trivial 解析值（锚点）")
    style = {
        "trivial": dict(color="k", ls="", marker=".", ms=6),
        "AA-FIR-1/DPW-2": dict(color="tab:red", ls="-", lw=1),
        "AA-FIR-2/DPW-3": dict(color="tab:blue", ls="-", lw=1),
        "OVS-2": dict(color="tab:red", ls="--", lw=1),
        "OVS-8": dict(color="tab:blue", ls="--", lw=1),
        "AA-IIR-1": dict(color="0.45", ls="-", lw=2.5),
        "AA-IIR-2": dict(color="tab:blue", ls="-", lw=2.5),
    }
    for name in names:
        ax.semilogx(notes[keep], snr[name][keep], label=name, **style[name])
        if (~keep).any():
            ax.semilogx(notes[~keep], snr[name][~keep], **style[name], alpha=0.25, label=None)
    ax.set_xlabel("基频 f0 [Hz]（88 键钢琴；淡色点 = 病态音，未计入平均）")
    ax.set_ylabel("SNR [dB]")
    ax.set_title("锯齿波 SNR（复现论文 Fig. 2 的排列方式）")
    ax.set_ylim(bottom=0)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    png = out / "fig2_saw_snr.png"
    fig.savefig(png, dpi=150)
    print(f"\n已写出 {png}\n已写出 {out / 'saw_snr.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
