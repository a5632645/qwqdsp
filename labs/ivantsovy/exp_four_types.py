"""LP / HP / BP / BS 四种经典类型的 §4 一般形式离散化与幅度对照。

两组原型：

- Butterworth N=4；
- Chebyshev II N=8、rs=80 dB（带通 / 带阻为 8→16 阶）。

fs = 48 kHz、f_ref = 5 kHz、sigma = pi*sqrt(2/3) + 预畸变。每个家族一张 2x2
图，分别画模拟原型与离散滤波器。

注意：scipy 的 zpk 增益经 ``normalized_gain`` 换算后才能作为“原点零点 +
常数项 1 因子”形式的增益；Chebyshev II 低通的 zpk 增益是 1e-4，直接使用会
整体偏移 80 dB。

用法: python qwqdsp/labs/ivantsovy/exp_four_types.py
产物: output/four_types.png、output/four_types_cheby2.png
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

import ivantsovy as I

FS = 48000.0
F_REF = 5000.0
SIGMA = np.pi * np.sqrt(2.0 / 3.0)
WC = 2.0 * np.pi * F_REF
OUT = Path(__file__).parent / "output"

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False


@dataclass(frozen=True)
class Family:
    slug: str
    name: str
    order: int
    rs: float | None
    png: str


FAMILIES = (
    Family("butterworth", "Butterworth", 4, None, "four_types.png"),
    Family("cheby2", "Chebyshev II", 8, 80.0, "four_types_cheby2.png"),
)

SPECS = (
    ("低通 LP", "lowpass", WC),
    ("高通 HP", "highpass", WC),
    ("带通 BP", "bandpass", [2.0 * np.pi * 4000.0, 2.0 * np.pi * 6000.0]),
    ("带阻 BS", "bandstop", [2.0 * np.pi * 4000.0, 2.0 * np.pi * 6000.0]),
)


def to_db(magnitude: np.ndarray) -> np.ndarray:
    return np.maximum(20.0 * np.log10(np.maximum(magnitude, 1e-300)), -240.0)


def subtitle(family: Family, btype: str) -> str:
    rs = f"，rs={family.rs:g} dB" if family.rs is not None else ""
    if btype in ("lowpass", "highpass"):
        edge = "fc" if family.rs is None else "阻带边沿"
        return f"{family.order} 阶 {family.name}{rs}，{edge} = 5 kHz"
    return f"{family.order} 阶原型 → {2 * family.order} 阶{rs}，4–6 kHz"


def scaled_prototype(family: Family, btype: str, wn):
    if family.rs is None:
        z, p, k = signal.butter(family.order, wn, btype=btype, analog=True,
                                output="zpk")
    else:
        z, p, k = signal.cheby2(family.order, family.rs, wn, btype=btype,
                                analog=True, output="zpk")
    z0 = np.asarray(z, dtype=complex)
    p0 = np.asarray(p, dtype=complex)
    z = z0 / WC
    p = p0 / WC
    k_s = k * WC ** (len(z0) - len(p0))
    gain = I.normalized_gain(z, p, k_s)
    return z, p, gain


def main() -> int:
    f_plot = np.geomspace(20.0, FS / 2.0, 16000)
    f_measure = np.linspace(0.0, FS / 2.0, 100001)

    for family in FAMILIES:
        fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.0))
        print(f"\n[{family.name}]")
        for ax, (title, btype, wn) in zip(axes.flat, SPECS):
            zeros, poles, gain = scaled_prototype(family, btype, wn)
            design = I.design(zeros, poles, FS, F_REF, SIGMA,
                              prewarp=True, gain=gain)
            analog_db = to_db(I.analog_magnitude(zeros, poles, f_plot, F_REF, gain))
            digital_db = to_db(I.digital_magnitude(design, f_plot, FS))
            analog = I.analog_magnitude(zeros, poles, f_measure, F_REF, gain)
            digital = I.digital_magnitude(design, f_measure, FS)
            max_linear = float(np.max(np.abs(digital - analog)))
            db_err = 20.0 * np.log10(np.maximum(digital, 1e-300)) - 20.0 * np.log10(
                np.maximum(analog, 1e-300))
            mask = analog >= 10.0 ** (-20.0 / 20.0)
            max_db = float(np.max(np.abs(db_err[mask])))

            ax.semilogx(f_plot, analog_db, color="C0", lw=1.6, label="模拟原型")
            ax.semilogx(f_plot, digital_db, color="C1", lw=1.4, ls="--",
                        label=f"Ivantsov 离散（σ={SIGMA:.3f}，预畸变）")
            ax.axvline(F_REF, color="0.5", lw=0.8, ls=":", zorder=0)
            ax.set_xlim(20.0, FS / 2.0)
            ax.set_ylim(-200.0, 5.0)
            ax.set_xlabel("频率 [Hz]")
            ax.set_ylabel("幅度 [dB]")
            ax.set_title(f"{title}  |  {subtitle(family, btype)}", fontsize=10.2)
            ax.grid(True, which="major", alpha=0.35)
            ax.grid(True, which="minor", alpha=0.15)
            ax.legend(loc="lower left", fontsize=8.3)
            ax.text(0.97, 0.05,
                    f"max |Hd−Ha| = {max_linear:.1e}\n"
                    f"模拟 ≥ −20 dB：≤ {max_db:.3f} dB",
                    transform=ax.transAxes, fontsize=8.2, ha="right", va="bottom",
                    bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.85,
                              edgecolor="0.75"))
            print(f"  {title:6} max |Hd-Ha| = {max_linear:.3e}  "
                  f"max dB (analog >= -20 dB) = {max_db:.3f}")

        fig.tight_layout()
        OUT.mkdir(exist_ok=True)
        path = OUT / family.png
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"  图：{path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
