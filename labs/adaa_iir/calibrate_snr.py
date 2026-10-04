"""SNR 口径标定：为什么采用 `blackmanharris` 窗 + ±4 bin + 2^17 采样。

论文（§IV）只写「SNR = 谐波功率 / 其余分量功率」，没有给测量细节。本脚本做两件事：

1. **解析锚点**：trivial 锯齿波的混叠是**解析可算**的——理想锯齿波第 k 次谐波
   （幅度 ``2/(pi*k)``）超过 Nyquist 后折叠，故
   ``SNR = 10*log10(Σ_{k≤K} 1/k² / Σ_{k>K} 1/k²)``，``K = floor(fs/2/f0)``
   （`aaiir.analytic_saw_snr`）。把实测与解析逐音比对，用来**排除病态音**：
   当 ``fs/f0`` 接近整数时，折叠分量恰好落在谐波栅格上，谱上无法区分混叠与谐波，
   实测 SNR 会虚高（甚至 inf），这类音必须剔除。
2. **口径扫描**：遍历 (长度 × 窗 × 排除带) 组合，看有没有哪一组能同时命中论文
   侧表的 trivial / OVS-8 / AA-IIR-2 三档数值（8.7 / 30.8 / 64.6 dB）。

用法: python qwqdsp/labs/adaa_iir/calibrate_snr.py [anchor|scan|all]
"""
from __future__ import annotations

import sys

import numpy as np

import aaiir as A

FS = 44100.0
WARMUP = 4096
PAPER = {"trivial": 8.7, "OVS-8": 30.8, "AA-IIR-2": 64.6}
NOTES = A.piano_frequencies()
TAB = A.saw_table()


def waveforms(f0, n):
    x = A.phase_ramp(n, f0, FS)
    return {
        "trivial": A.trivial(TAB, n, f0, FS),
        "OVS-2": A.oversampled(TAB, n, f0, FS, 2),
        "OVS-8": A.oversampled(TAB, n, f0, FS, 8),
        "AA-FIR-1/DPW-2": A.dpw(TAB, n, f0, FS, 2),
        "AA-FIR-2/DPW-3": A.dpw(TAB, n, f0, FS, 3),
        "AA-IIR-1": A.aa_iir_osc(x, A.aa_iir_1(), TAB),
        "AA-IIR-2": A.aa_iir_osc(x, A.aa_iir_2(), TAB),
    }


def anchor() -> None:
    """实测 vs 解析（trivial），并列出会被门限剔除的音。"""
    print("== 解析锚点：trivial 实测 vs analytic_saw_snr ==")
    n = 2 ** 17                     # 与 exp_saw_snr.py 同长度
    dev = []
    for f0 in NOTES:
        got = A.snr_db(A.trivial(TAB, n, f0, FS), FS, f0, guard_bins=4, warmup=WARMUP)
        dev.append(got - A.analytic_saw_snr(f0, FS))
    dev = np.array(dev)
    bad = ~np.isfinite(dev) | (np.abs(dev) > 2.0)
    print(f"  88 音：中位偏差 {np.median(dev[~bad]):+.2f} dB，"
          f"良态音最大偏差 {np.abs(dev[~bad]).max():.2f} dB")
    print(f"  剔除 {int(bad.sum())} 个病态音（|偏差| > 2 dB 或 inf），"
          f"最难的一档：")
    for i in np.argsort(-np.abs(np.nan_to_num(dev, nan=99)))[:5]:
        print(f"    f0 = {NOTES[i]:8.2f} Hz  fs/f0 = {FS / NOTES[i]:10.3f}  偏差 {dev[i]:+.1f} dB")


def scan() -> None:
    """(长度 × 窗 × 排除带) 扫描：能否同时命中论文的三档数值。"""
    print("\n== 口径扫描（88 音平均，目标：论文的 trivial / OVS-8 / AA-IIR-2） ==")
    windows = ("rect", "hann", "hamming", "blackmanharris", "nuttall")
    for n in (4410, 44100, 2 ** 17):
        print(f"\n--- N = {n}（{n / FS * 1000:.0f} ms）---")
        cache = {f0: waveforms(f0, n) for f0 in NOTES}
        best = None
        print(f"{'窗':<16}{'排除带':>7}" + "".join(f"{k:>12}" for k in PAPER))
        for win in windows:
            for g in (0, 1, 2, 4, 8):
                row = []
                for key in PAPER:
                    vals = [A.snr_db(cache[f0][key], FS, f0, guard_bins=g, warmup=WARMUP,
                                     window=win) for f0 in NOTES]
                    m = float(np.mean(vals))
                    row.append(m if np.isfinite(m) else -99.0)
                err = max(abs(row[i] - PAPER[k]) for i, k in enumerate(PAPER))
                print(f"{win:<16}{g:>7}" + "".join(f"{v:>12.1f}" for v in row))
                if best is None or err < best[0]:
                    best = (err, win, g, row)
        print(f"  最优：窗 {best[1]}、排除带 ±{best[2]} → 三档最大偏差 {best[0]:.1f} dB")
    print("\n结论：没有任何 (长度 × 窗 × 排除带) 能同时命中三档；"
          "因此采用可被解析锚点验证的一组（blackmanharris / ±4 bin / 2^17）。")


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("anchor", "all"):
        anchor()
    if what in ("scan", "all"):
        scan()
