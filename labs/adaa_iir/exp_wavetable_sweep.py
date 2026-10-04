"""论文 §IV-C/§IV-D：Escalation II 波形的扫频测试（对应 Fig. 3 的 5 张谱图）。

做法：对一串指数分布的基频 f0（20 Hz → 10 kHz），每种方法各生成一小段定频信号，
取中间一段加 Hann 窗做 FFT，按 dB 堆成二维图（横轴 f0，纵轴频率），下限 -80 dB。
论文的 Kangaroo（商业 wavetable）无法获取，故只复现 Escalation II（见 README）。

用法: python qwqdsp/labs/adaa_iir/exp_wavetable_sweep.py
输出: qwqdsp/labs/adaa_iir/output/fig3_escalation_sweep.png
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
N_F0 = 256
F0_MIN, F0_MAX = 20.0, 10000.0
N_CHUNK = 4410          # 每档 0.1 s
SKIP = 1000             # 跳过 AA 递推启动瞬态
N_FFT = 2048
N_PAD = 8192
FLOOR_DB = -80.0

METHODS = ("trivial", "OVS-2", "OVS-8", "AA-IIR-1", "AA-IIR-2")


def spectrum(y, f0):
    """取缓中段做 Hann 加窗 FFT，返回幅度谱（已按窗增益归一）。"""
    seg = y[SKIP:SKIP + N_FFT]
    seg = seg - seg.mean()
    w = np.hanning(len(seg))
    spec = np.abs(np.fft.rfft(seg * w, n=N_PAD)) / (w.sum() / 2.0)
    return spec


def spurious_db(spec, f0, guard_hz=60.0):
    """最坏非谐波分量相对基波的 dB（谐波 ±guard_hz 之外的最大谱线）。

    低音区谐波间隔小于 2*guard_hz 时谐波带互相重叠、覆盖全谱，此时返回 nan。
    """
    freqs = np.fft.rfftfreq(N_PAD, 1.0 / FS)
    mask = freqs < guard_hz                       # 去掉 DC 附近
    k = 1
    while k * f0 < FS / 2:
        mask |= np.abs(freqs - k * f0) <= guard_hz
        k += 1
    rest = spec[~mask]
    if rest.size == 0:
        return np.nan
    base = spec[freqs <= f0 * (1 + 0.05)].max()    # 基波峰值（含邻近泄漏）
    return 20.0 * np.log10(rest.max() / base)


def sweep(method: str, tab: A.PwlTable, f0s):
    """扫频二维谱：shape = (len(f0s), N_PAD//2+1)。"""
    out = []
    for f0 in f0s:
        x = A.phase_ramp(N_CHUNK, f0, FS)
        if method == "trivial":
            y = A.trivial(tab, N_CHUNK, f0, FS)
        elif method.startswith("OVS"):
            y = A.oversampled(tab, N_CHUNK, f0, FS, int(method[4:]))
        elif method == "AA-IIR-1":
            y = A.aa_iir_osc(x, A.aa_iir_1(), tab)
        elif method == "AA-IIR-2":
            y = A.aa_iir_osc(x, A.aa_iir_2(), tab)
        elif method == "AA-FIR-2":
            y = A.aa_fir(tab, N_CHUNK, f0, FS, 2)
        else:
            raise ValueError(method)
        out.append(spectrum(y, f0))
    return np.array(out)


def main() -> int:
    tab = A.escalation_ii_w3_table()
    f0s = np.exp(np.linspace(np.log(F0_MIN), np.log(F0_MAX), N_F0))
    freqs = np.fft.rfftfreq(N_PAD, 1.0 / FS) / 1000.0     # [kHz]

    specs = {m: sweep(m, tab, f0s) for m in METHODS}
    peak = max(float(s.max()) for s in specs.values())    # 全体共用同一参考（便于跨面板比较）

    fig, axes = plt.subplots(len(METHODS), 1, figsize=(7.5, 11), sharex=True, sharey=True)
    for ax, m in zip(axes, METHODS):
        db = 20.0 * np.log10(np.maximum(specs[m] / peak, 10 ** (FLOOR_DB / 20)))
        db = np.maximum(db, FLOOR_DB)
        mesh = ax.pcolormesh(f0s, freqs, db.T, cmap="magma", vmin=FLOOR_DB, vmax=0.0,
                             shading="auto")
        ax.set_title(m, fontsize=9, loc="left")
        ax.set_ylim(0, 20)
    axes[-1].set_xlabel("基频 f0 [Hz]（指数扫频）")
    for ax in axes:
        ax.set_ylabel("频率 [kHz]", fontsize=8)
    fig.suptitle("Fig. 3 复现：Escalation II 扫频谱图（下限 -80 dB）")
    fig.colorbar(mesh, ax=axes, label="幅度 [dB]", pad=0.02)
    fig.savefig(Path(__file__).parent / "output" / "fig3_escalation_sweep.png", dpi=150)
    print("每档 f0 的「最坏非谐波分量 / 基波」[dB]（低音区谐波带重叠，已跳过）：")
    print(f"{'方法':<12}{'有效档数':>10}{'中位数':>10}{'最差':>10}")
    for m in METHODS:
        s = np.array([spurious_db(specs[m][i], f0) for i, f0 in enumerate(f0s)])
        ok = np.isfinite(s)
        print(f"{m:<12}{int(ok.sum()):>10}{np.median(s[ok]):>10.1f}{s[ok].max():>10.1f}")
    print("已写出 output/fig3_escalation_sweep.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
