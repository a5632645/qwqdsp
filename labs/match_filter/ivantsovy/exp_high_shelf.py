"""High-shelf 的 §4 一般形式离散化与幅度对照。

四个面板：Butterworth N=8 的 +12 / −12 dB，以及 Chebyshev II N=8、ripple=1 dB
的 +12 / −12 dB。fs = 48 kHz、f0 = 1 kHz、sigma = pi*sqrt(2/3) + 预畸变。

Chebyshev II 使用 ``even_modify=False``（scipy 没有偶数阶修正原型），因此
高频端在 ``boost - sign(boost)*ripple`` dB（面板里约 ±11 dB），棚带内的极值
仍在有限频率处到达 ±12 dB。

用法: python qwqdsp/labs/match_filter/ivantsovy/exp_high_shelf.py
产物: output/high_shelf.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import ivantsovy as I
import shelf as S

FS = 48000.0
F0 = 1000.0
SIGMA = np.pi * np.sqrt(2.0 / 3.0)
ORDER = 8
RIPPLE = 1.0
OUT = Path(__file__).parent / "output"

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

CASES = (
    ("butter", "Butterworth +12 dB", "Butterworth N=8 high shelf，+12 dB，f0 = 1 kHz", 12.0),
    ("butter", "Butterworth −12 dB", "Butterworth N=8 high shelf，−12 dB，f0 = 1 kHz", -12.0),
    ("cheb2", "Chebyshev II +12 dB",
     "Chebyshev II N=8，ripple=1 dB，+12 dB，f0 = 1 kHz（高频端 ≈11 dB）", 12.0),
    ("cheb2", "Chebyshev II −12 dB",
     "Chebyshev II N=8，ripple=1 dB，−12 dB，f0 = 1 kHz（高频端 ≈−11 dB）", -12.0),
)


def to_db(magnitude: np.ndarray) -> np.ndarray:
    return np.maximum(20.0 * np.log10(np.maximum(magnitude, 1e-300)), -60.0)


def build(kind: str, boost: float):
    if kind == "butter":
        zeros, poles = S.butterworth_shelf(ORDER // 2, boost)
    else:
        zeros, poles = S.chebyshev2_shelf(ORDER // 2, boost, RIPPLE)
    return zeros, poles


def main() -> int:
    f_plot = np.geomspace(20.0, FS / 2.0, 16000)
    f_measure = np.linspace(20.0, FS / 2.0, 100001)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.0))
    for ax, (kind, short, title, boost) in zip(axes.flat, CASES):
        zeros, poles = build(kind, boost)
        design = I.design(zeros, poles, FS, F0, SIGMA, gain=1.0, prewarp=True)
        analog_db = to_db(I.analog_magnitude(zeros, poles, f_plot, F0))
        digital_db = to_db(I.digital_magnitude(design, f_plot, FS))
        analog = I.analog_magnitude(zeros, poles, f_measure, F0)
        digital = I.digital_magnitude(design, f_measure, FS)
        db_err = 20.0 * np.log10(np.maximum(digital, 1e-300)) - 20.0 * np.log10(
            np.maximum(analog, 1e-300))
        max_db = float(np.max(np.abs(db_err)))

        ax.semilogx(f_plot, analog_db, color="C0", lw=1.6, label="模拟原型")
        ax.semilogx(f_plot, digital_db, color="C1", lw=1.4, ls="--",
                    label=f"Ivantsov 离散（σ={SIGMA:.3f}，预畸变）")
        ax.axvline(F0, color="0.5", lw=0.8, ls=":", zorder=0)
        ax.axhline(0.0, color="0.8", lw=0.6, zorder=0)
        ax.set_xlim(20.0, FS / 2.0)
        ax.set_ylim(-24.0, 15.0)
        ax.set_xlabel("频率 [Hz]")
        ax.set_ylabel("幅度 [dB]")
        ax.set_title(title, fontsize=9.6)
        ax.grid(True, which="major", alpha=0.35)
        ax.grid(True, which="minor", alpha=0.15)
        ax.legend(loc="lower left", fontsize=8.3)
        ax.text(0.97, 0.05, f"max |ΔdB| = {max_db:.3f} dB",
                transform=ax.transAxes, fontsize=8.5, ha="right", va="bottom",
                bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.85,
                          edgecolor="0.75"))
        print(f"{short:22} max |ΔdB| = {max_db:.4f} dB")

    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    path = OUT / "high_shelf.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"图：{path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
