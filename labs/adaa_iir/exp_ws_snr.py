"""论文 [19] §4.1 复现：硬削波的混叠抑制对比（对应其 Fig. 3）。

输入：增益 10 的正弦（每个钢琴键频率），过硬削波 ``clip(x, ±1)``；SNR = 谐波功率 / 其余分量
（排除 DC），统计 88 键。方法与 [19] Table 1 一致：

  trivial / OVS-2 / OVS-8 / AA-FIR-1（矩形核）/ AA-FIR-2（三角核）/ AA-IIR-1 / AA-IIR-2

**有效性门限**：与振荡器那篇同样的做法 —— 低音区折叠分量会正好落在谐波栅格上（`fs/f0` 接近
整数），谱上无法区分混叠与谐波，实测 SNR 虚高。用**解析锚点** `ws.analytic_clip_snr`
（削波正弦的解析傅里叶系数 → 折叠功率比）判掉这类音，剩余音集对全部方法统一使用。

用法: python qwqdsp/labs/adaa_iir/exp_ws_snr.py
输出: qwqdsp/labs/adaa_iir/output/ws_snr_hardclip.png（+ .csv）
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
import ws as W

FS = 44100.0
N = 2 ** 16
GAIN = 10.0
WARMUP = 4096
GUARD = 4
GATE_DB = 2.0

# [19] Fig. 3 的定性结论（论文未给数值表，只能按顺序核对）
ORDER_CLAIM = "trivial < OVS-2 < AA-FIR-1 ≲ AA-IIR-1 < AA-FIR-2 < OVS-8 ≲ AA-IIR-2"


def run_methods(x, f):
    """返回 (方法名 -> 输出)。"""
    return {
        "trivial": W.trivial(x, f),
        "OVS-2": W.oversampled(x, 2, f),
        "OVS-8": W.oversampled(x, 8, f),
        "AA-FIR-1": W.aa_fir1(x, f),
        "AA-FIR-2": W.aa_fir2(x, f),
        "AA-IIR-1": W.aa_iir_pwl(x, W.aa_iir_1(), f),
        "AA-IIR-2": W.aa_iir_pwl(x, W.aa_iir_2(), f),
    }


def main() -> int:
    f = W.hard_clip(1.0)
    notes = A.piano_frequencies()
    names = ["trivial", "OVS-2", "OVS-8", "AA-FIR-1", "AA-FIR-2", "AA-IIR-1", "AA-IIR-2"]
    snr = {k: [] for k in names}
    for i, f0 in enumerate(notes):
        x = GAIN * np.sin(2.0 * np.pi * f0 * np.arange(N) / FS)
        for k, y in run_methods(x, f).items():
            snr[k].append(A.snr_db(y, FS, f0, guard_bins=GUARD, warmup=WARMUP))
        if (i + 1) % 22 == 0:
            print(f"  {i + 1}/88 键完成")

    anchor = np.array([W.analytic_clip_snr(f0, FS, GAIN, 1.0) for f0 in notes])
    snr = {k: np.array(v) for k, v in snr.items()}
    dev = snr["trivial"] - anchor
    keep = np.isfinite(dev) & (np.abs(dev) <= GATE_DB)
    print(f"\n有效性门限：剔除 {int((~keep).sum())}/88 个音（trivial 实测 vs 解析偏离 > "
          f"{GATE_DB:.0f} dB）：")
    print("  " + ", ".join(f"{v:.1f}" for v in notes[~keep]))

    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    with open(out / "ws_snr_hardclip.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["f0_hz", "analytic_trivial", "kept"] + names)
        for i, f0 in enumerate(notes):
            wr.writerow([f"{f0:.4f}", f"{anchor[i]:.3f}", int(keep[i])]
                        + [f"{snr[k][i]:.3f}" for k in names])

    print(f"\n平均 SNR [dB]（{int(keep.sum())} 个音；[19] 未给数值表，只核对排序）")
    print(f"{'方法':<12}{'平均':>8}{'中位':>8}{'最小':>8}")
    for k in names:
        v = snr[k][keep]
        print(f"{k:<12}{v.mean():>8.1f}{np.median(v):>8.1f}{v.min():>8.1f}")
    print(f"\n论文定性排序（[19] §4.1）：{ORDER_CLAIM}")

    fig, ax = plt.subplots(figsize=(9, 5))
    style = {
        "trivial": dict(color="k", ls=":", lw=1.2),
        "OVS-2": dict(color="tab:green", ls="--", lw=1),
        "OVS-8": dict(color="tab:green", ls="-.", lw=1),
        "AA-FIR-1": dict(color="tab:red", ls="-", lw=1),
        "AA-FIR-2": dict(color="tab:orange", ls="-", lw=1),
        "AA-IIR-1": dict(color="0.45", ls="-", lw=2.2),
        "AA-IIR-2": dict(color="tab:blue", ls="-", lw=2.2),
    }
    ax.semilogx(notes, anchor, "k", lw=0.8, alpha=0.5, label="trivial 解析值（锚点）")
    for k in names:
        ax.semilogx(notes[keep], snr[k][keep], label=k, **style[k])
    ax.set_xlabel("基频 f0 [Hz]（88 键钢琴，正弦增益 10 过硬削波）")
    ax.set_ylabel("SNR [dB]")
    ax.set_title("复现 [19] Fig. 3：硬削波混叠抑制对比")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    png = out / "ws_snr_hardclip.png"
    fig.savefig(png, dpi=150)
    print(f"\n已写出 {png}\n已写出 {out / 'ws_snr_hardclip.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
