"""双点匹配的一 / 二阶测试图：5 个测试截止频率，整条响应。

测试尺度统一为 ``fc = 20 / 117.6 / 691.4 / 4065 / 23900 Hz``（对数 5 点、含端点，
fs 一律 48 kHz）；每张图画整条幅度响应（模拟实线、双点解虚线）与两者的偏差。
``fc = 30 kHz`` 及以上在 48 kHz 下无解（截止频率达不到 Nyquist），只做数值断言。

用法: python qwqdsp/labs/match_filter/ivantsovy/exp_two_point.py
产物: output/two_point_scales.png、output/two_point_lowpass.png、
      output/two_point_shelf.png、output/two_point_shelf2.png、
      output/two_point_second_order.png、output/two_point_bands.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import ivantsovy as I
import two_point as TP

CASES = tuple((float(v), 48000.0) for v in np.geomspace(20.0, 23900.0, 5))
"""测试截止频率：20 Hz … 23.9 kHz，对数 5 点（含端点）。"""

POLES = np.array([-1.0])
HP_ZEROS = np.array([0.0])
LP_ZEROS = np.array([])
OUT = Path(__file__).parent / "output"
SIGMA_OPT = np.pi * np.sqrt(2.0 / 3.0)
ZETA = 1.0 / np.sqrt(2.0)
COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False


def to_db(magnitude: np.ndarray) -> np.ndarray:
    return np.maximum(20.0 * np.log10(np.maximum(magnitude, 1e-300)), -240.0)


def case_label(fc: float) -> str:
    return f"fc={fc / 1000:.3g} kHz" if fc >= 1000 else f"fc={fc:.4g} Hz"


def grid(fs: float) -> np.ndarray:
    return np.geomspace(20.0, fs / 2.0, 12000)


def first_order(zeros, build, fc: float, fs: float):
    f = grid(fs)
    design = build(fc, fs)
    return (f, I.analog_magnitude(zeros, POLES, f, fc),
            I.digital_magnitude(design, f, fs))


def shelf(kind: str, boost_db: float, fc: float, fs: float):
    f = grid(fs)
    build = TP.high_shelf if kind == "high" else TP.low_shelf
    zeros, poles, gain = TP.shelf_prototype(kind, boost_db)
    return (f, I.analog_magnitude(zeros, poles, f, fc, gain),
            I.digital_magnitude(build(fc, fs, boost_db), f, fs))


SHELF2_QS = (np.sqrt(2.0) / 2.0, 6.0)


def shelf2(kind: str, boost_db: float, zeta: float, fc: float, fs: float):
    f = grid(fs)
    build = TP.high_shelf2 if kind == "high" else TP.low_shelf2
    zeros, poles, gain = TP.shelf_prototype2(kind, boost_db, zeta)
    return (f, I.analog_magnitude(zeros, poles, f, fc, gain),
            I.digital_magnitude(build(fc, fs, boost_db, zeta), f, fs))


def second_order(kind: str, fc: float, fs: float):
    f = grid(fs)
    x = f / fc
    b = 2.0 * ZETA ** 2 - 1.0
    den = np.sqrt(x ** 4 + 2.0 * x ** 2 * b + 1.0)
    analog = x ** 2 / den if kind == "hp" else 1.0 / den
    build = TP.highpass2 if kind == "hp" else TP.lowpass2
    return f, analog, I.digital_magnitude(build(fc, fs, ZETA), f, fs)


def band(kind: str, fc: float, fs: float):
    f = grid(fs)
    s = 1j * f / fc
    den = np.abs(s ** 2 + 2.0 * ZETA * s + 1.0)
    analog = np.abs(2.0 * ZETA * s) / den if kind == "bp" else np.abs(s ** 2 + 1.0) / den
    build = TP.bandpass if kind == "bp" else TP.bandstop
    return f, analog, I.digital_magnitude(build(fc, fs, ZETA), f, fs)


def panel_magnitude(ax, curves, title: str, floor: float | None = None) -> None:
    values: list[float] = []
    for (fc, fs), color in zip(CASES, COLORS):
        f, analog, digital = curves(fc, fs)
        adb, ddb = to_db(analog), to_db(digital)
        if floor is not None:
            adb, ddb = np.maximum(adb, floor), np.maximum(ddb, floor)
        ax.semilogx(f, adb, color=color, lw=1.1)
        ax.semilogx(f, ddb, color=color, lw=1.5, ls="--", label=case_label(fc))
        ax.axvline(fc, color=color, lw=0.6, ls=":", alpha=0.6)
        values += list(adb[::20]) + list(ddb[::20])
    ax.set_xlim(20.0, 30000.0)
    ax.set_ylim(min(values) - 3.0, max(values) + 3.0)
    ax.set_xlabel("频率 [Hz]")
    ax.set_ylabel("幅度 [dB]")
    ax.set_title(f"{title}：实线模拟原型，虚线双点匹配")
    ax.grid(True, which="major", alpha=0.3)
    ax.grid(True, which="minor", alpha=0.12)
    ax.legend(loc="lower left" if title.startswith(("一阶高通", "二阶高通", "二阶带通"))
              else "lower right", fontsize=8)


def panel_deviation(ax, curves, title: str, note: str,
                    mask: float | None = None) -> None:
    values: list[float] = []
    for (fc, fs), color in zip(CASES, COLORS):
        f, analog, digital = curves(fc, fs)
        dev = to_db(digital) - to_db(analog)
        if mask is not None:
            # 陷波两侧 |Ha| -> 0，dB 比值无意义；只在 |Ha| >= mask 处绘制。
            dev = np.where(analog >= mask, dev, np.nan)
        ax.semilogx(f, dev, color=color, lw=1.4, label=case_label(fc))
        values += [v for v in dev[::20] if np.isfinite(v)]
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_xlim(20.0, 30000.0)
    ax.set_ylim(min(values) - 1.0, max(values) + 1.0)
    ax.set_xlabel("频率 [Hz]")
    ax.set_ylabel("Hd − Ha [dB]")
    ax.set_title(f"{title}：偏差（匹配点处为 0）")
    ax.grid(True, which="major", alpha=0.3)
    ax.grid(True, which="minor", alpha=0.12)
    ax.legend(loc="best", fontsize=8)
    ax.text(0.02, 0.97, note, transform=ax.transAxes, fontsize=8.2, va="top",
            ha="left",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      alpha=0.85, edgecolor="0.75"))


def panel_sigma(ax) -> None:
    u = np.geomspace(2.0 + 1e-9, 100000.0, 4000)
    sigma = np.array([TP.omega_sigma(48000.0 / v, 48000.0)[1] for v in u])
    ax.plot(1.0 / u, sigma, color="C0", lw=1.8, label="双点匹配的 σ(fc/fs)")
    ax.axhspan(2.0, SIGMA_OPT, color="C2", alpha=0.15)
    ax.axhline(2.0, color="C2", lw=0.9, ls="--")
    ax.text(0.02, SIGMA_OPT + 0.05, "论文幅度最优区间 [2, π√(2/3)]",
            color="C2", fontsize=8.5, va="bottom")
    for fc, fs in CASES:
        _, sigma = TP.omega_sigma(fc, fs)
        ax.plot([fc / fs], [sigma], "o", color="C0", ms=4.5, zorder=3)
    ax.set_xlim(0.0, 0.5)
    ax.set_ylim(0.0, 2.8)
    ax.set_xlabel("fc / fs")
    ax.set_ylabel("σ")
    ax.set_title("σ 被双点条件唯一确定，恒 ≤ 2")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=8.5)


def save(fig, name: str) -> Path:
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    path = OUT / name
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main() -> int:
    print("测试截止频率 [Hz]：" + ", ".join(f"{fc:.1f}" for fc, _ in CASES))

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.7))
    panel_magnitude(axes[0], lambda fc, fs: first_order(HP_ZEROS, TP.highpass, fc, fs),
                    "一阶高通")
    panel_deviation(axes[1], lambda fc, fs: first_order(HP_ZEROS, TP.highpass, fc, fs),
                    "一阶高通",
                    "偏差以低频最大（= 20log10(ω/u)）；\n"
                    "fc 越接近 Nyquist 塌得越深")
    panel_sigma(axes[2])
    print(f"图：{save(fig, 'two_point_scales.png')}")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    panel_magnitude(axes[0], lambda fc, fs: first_order(LP_ZEROS, TP.lowpass, fc, fs),
                    "一阶低通")
    panel_deviation(axes[1], lambda fc, fs: first_order(LP_ZEROS, TP.lowpass, fc, fs),
                    "一阶低通",
                    "DC 精确（三点精确）；残差集中在 Nyquist 附近")
    print(f"图：{save(fig, 'two_point_lowpass.png')}")

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))
    panel_magnitude(axes[0][0], lambda fc, fs: shelf("high", 12.0, fc, fs),
                    "一阶高 shelf +12 dB")
    panel_deviation(axes[0][1], lambda fc, fs: shelf("high", 12.0, fc, fs),
                    "一阶高 shelf",
                    "+12 dB；−12 dB 是其严格镜像\n（|H+b|·|H−b| = 1）")
    panel_magnitude(axes[1][0], lambda fc, fs: shelf("low", 12.0, fc, fs),
                    "一阶低 shelf +12 dB")
    panel_deviation(axes[1][1], lambda fc, fs: shelf("low", 12.0, fc, fs),
                    "一阶低 shelf",
                    "+12 dB；DC / 中点 / Nyquist 三点精确，\n残差集中在拐点附近")
    print(f"图：{save(fig, 'two_point_shelf.png')}")

    fig, axes = plt.subplots(4, 2, figsize=(11.5, 15.0))
    for row, (kind, label) in enumerate((("high", "二阶高 shelf"), ("low", "二阶低 shelf"))):
        for sub, Q in enumerate(SHELF2_QS):
            zeta = 1.0 / (2.0 * Q)
            panel_magnitude(axes[row * 2 + sub][0],
                            lambda fc, fs, k=kind, z=zeta: shelf2(k, 12.0, z, fc, fs),
                            f"{label} +12 dB（Q={Q:.3g}）")
            panel_deviation(axes[row * 2 + sub][1],
                            lambda fc, fs, k=kind, z=zeta: shelf2(k, 12.0, z, fc, fs),
                            f"{label}（Q={Q:.3g}）",
                            "中点 / Nyquist / DC 三点精确；\nQ 越大，偏离模拟越明显")
    print(f"图：{save(fig, 'two_point_shelf2.png')}")

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))
    panel_magnitude(axes[0][0], lambda fc, fs: second_order("hp", fc, fs), "二阶高通")
    panel_deviation(axes[0][1], lambda fc, fs: second_order("hp", fc, fs), "二阶高通",
                    "低频偏移 = 40log10(ω/u) dB\n（一阶的两倍，每个原点零点一份 ω/u）")
    panel_magnitude(axes[1][0], lambda fc, fs: second_order("lp", fc, fs), "二阶低通")
    panel_deviation(axes[1][1], lambda fc, fs: second_order("lp", fc, fs), "二阶低通",
                    "DC 精确；残差集中在 Nyquist 附近")
    print(f"图：{save(fig, 'two_point_second_order.png')}")

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))
    panel_magnitude(axes[0][0], lambda fc, fs: band("bp", fc, fs), "二阶带通")
    panel_deviation(axes[0][1], lambda fc, fs: band("bp", fc, fs), "二阶带通",
                    "峰值精确在 fc；低频侧同高通\n多衰减 20log10(ω/u) dB")
    panel_magnitude(axes[1][0], lambda fc, fs: band("bs", fc, fs), "二阶带阻",
                    floor=-80.0)
    panel_deviation(axes[1][1], lambda fc, fs: band("bs", fc, fs), "二阶带阻",
                    "陷波精确在 fc（|H(fc)| ≈ 1e-16），零点在单位圆上；\n"
                    "只画 |Ha| ≥ 0.05 处（陷波两侧比值无意义）", mask=0.05)
    print(f"图：{save(fig, 'two_point_bands.png')}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
