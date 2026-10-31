"""高阶 Butterworth / Chebyshev II 的 §4 一般形式离散化与幅度对照。

两组配置：

1. ``sigma = pi*sqrt(2/3)``、不做预畸变：论文的幅度最优区间上端点，整体
   幅度误差最小。
2. ``sigma = 2``、做预畸变：``Hd(f_ref)`` 精确等于 ``Ha(f_ref)``，但 Nyquist
   的精确匹配点从 ``fs/sigma`` 移到 ``f_ref * omega_d / sigma``（本配置为
   23.67 kHz，而不是 24 kHz）。

每个案例、每种配置导出一张图，同一坐标系画模拟原型与离散滤波器。

用法: python qwqdsp/labs/ivantsovy/exp_magnitude.py
产物: output/*.png
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
ORDER = 8
SIGMA_OPT = np.pi * np.sqrt(2.0 / 3.0)
OUT = Path(__file__).parent / "output"

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False


@dataclass(frozen=True)
class Case:
    slug: str
    title: str
    zeros: np.ndarray
    poles: np.ndarray
    cutoff_label: str


@dataclass(frozen=True)
class Config:
    sigma: float
    prewarp: bool
    suffix: str
    label: str


CONFIGS = (
    Config(SIGMA_OPT, False, "", "σ=π√(2/3)，无预畸变"),
    Config(2.0, True, "_prewarp_sigma2", "σ=2 + 预畸变"),
)

COMPARISON = (
    Config(SIGMA_OPT, False, "", "2.565 / no prewarp"),
    Config(2.0, False, "", "2.000 / no prewarp"),
    Config(2.0, True, "", "2.000 / prewarp"),
    Config(SIGMA_OPT, True, "", "2.565 / prewarp"),
)


def make_cases() -> tuple[Case, ...]:
    z, p, _ = signal.buttap(ORDER)
    butter = Case("butterworth_N8", f"Butterworth  N={ORDER}, fc=5 kHz",
                  I.complete_conjugates(z), I.complete_conjugates(p), "fc")
    z, p, _ = signal.cheb2ap(ORDER, 80.0)
    cheby2 = Case("chebyshev2_N8", f"Chebyshev II  N={ORDER}, rs=80 dB, 阻带边沿 5 kHz",
                  I.complete_conjugates(z), I.complete_conjugates(p), "阻带边沿")
    return butter, cheby2


def to_db(magnitude: np.ndarray) -> np.ndarray:
    return np.maximum(20.0 * np.log10(np.maximum(magnitude, 1e-300)), -240.0)


def measure(case: Case, config: Config) -> tuple[I.Design, dict]:
    design = I.design(case.zeros, case.poles, FS, F_REF, config.sigma,
                      prewarp=config.prewarp)
    f = np.linspace(0.0, FS / 2.0, 200001)
    analog = I.analog_magnitude(case.zeros, case.poles, f, F_REF)
    digital = I.digital_magnitude(design, f, FS)
    dev = digital - analog
    db_dev = to_db(digital) - to_db(analog)

    omega = (I.prewarped_omega(F_REF, FS, config.sigma) if config.prewarp
             else FS / F_REF)
    f_match = F_REF * omega / config.sigma
    hd_ref = float(I.digital_magnitude(design, np.array([F_REF]), FS)[0])
    ha_ref = float(I.analog_magnitude(case.zeros, case.poles, np.array([F_REF]), F_REF)[0])
    ha_match = float(I.analog_magnitude(case.zeros, case.poles, np.array([f_match]), F_REF)[0])

    stats: dict[str, object] = {
        "max_linear": float(np.max(np.abs(dev))),
        "max_linear_hz": float(f[np.argmax(np.abs(dev))]),
        "dc_error": float(abs(digital[0] - analog[0])),
        "ref_db": float(20.0 * np.log10(hd_ref / ha_ref)),
        "hd_ref": hd_ref,
        "ha_ref": ha_ref,
        "omega": float(omega),
        "f_match": float(f_match),
        "ha_match": ha_match,
        "ha_nyquist": float(analog[-1]),
        "hd_nyquist": float(digital[-1]),
    }
    for threshold in (20, 40, 60, 80):
        mask = analog >= 10.0 ** (-threshold / 20.0)
        stats[f"db_{threshold}"] = float(np.max(np.abs(db_dev[mask])))

    zeros, poles = I.digital_zeros_poles(design)
    stats["pole_radius"] = float(np.max(np.abs(poles)))
    stats["zeros"] = zeros
    return design, stats


def describe_zeros(case: Case, zeros: np.ndarray) -> None:
    if not len(case.zeros):
        print(f"  数字零点：{len(zeros)} 个都来自 s -> infinity 的阶数亏欠，"
              f"位于 z = {zeros[0].real:+.4f}（Nyquist 轴，单位圆内）")
        return

    on_circle = zeros[np.abs(np.abs(zeros) - 1.0) <= 1e-9]
    real_axis = zeros[np.abs(zeros.imag) <= 1e-6]
    freq = np.sort(np.abs(np.angle(on_circle)) / (2.0 * np.pi) * FS)[::2]
    analog_hz = np.sort(np.abs(case.zeros) * F_REF)[::2]
    print("  模拟零点 [Hz]：", ", ".join(f"{v:.0f}" for v in analog_hz))
    print("  数字单位圆零点 [Hz]：", ", ".join(f"{v:.0f}" for v in freq))
    print("  数字实轴零点：", ", ".join(f"z={z.real:+.4f}" for z in real_axis),
          f"（对应模拟零点 {analog_hz[-1]:.0f} Hz，高于 Nyquist）")


def report(case: Case, config: Config, stats: dict) -> None:
    print(f"\n[{case.title} | {config.label}]  "
          f"sigma={config.sigma:.4f}, fs={FS:.0f} Hz, f_ref={F_REF:.0f} Hz")
    print(f"  参考频率：Hd/Ha = {stats['ref_db']:+.4f} dB "
          f"(Hd={stats['hd_ref']:.6e}, Ha={stats['ha_ref']:.6e})")
    print(f"  Nyquist：Hd={stats['hd_nyquist']:.4e}；"
          f"A(f_ref*omega/sigma={stats['f_match']:.0f} Hz)={stats['ha_match']:.4e}；"
          f"A(Nyquist)={stats['ha_nyquist']:.4e}")
    print(f"  全频段 max |Hd - Ha| = {stats['max_linear']:.3e} "
          f"@ {stats['max_linear_hz']:.1f} Hz")
    print("  max |Delta dB|："
          f"模拟 >= -20 dB: {stats['db_20']:.3f} dB，"
          f">= -40 dB: {stats['db_40']:.3f} dB，"
          f">= -60 dB: {stats['db_60']:.3f} dB，"
          f">= -80 dB: {stats['db_80']:.3f} dB")
    print(f"  DC 误差 {stats['dc_error']:.3e}；最大数字极点模 {stats['pole_radius']:.6f}")


def plot(case: Case, config: Config, design: I.Design, stats: dict) -> Path:
    f = np.geomspace(20.0, FS / 2.0, 16000)
    analog_db = to_db(I.analog_magnitude(case.zeros, case.poles, f, F_REF))
    digital_db = to_db(I.digital_magnitude(design, f, FS))

    fig, ax = plt.subplots(figsize=(9.5, 5.4))
    ax.semilogx(f, analog_db, color="C0", lw=1.7, label="模拟原型")
    ax.semilogx(f, digital_db, color="C1", lw=1.5, ls="--",
                label=f"Ivantsov 离散（{config.label}）")
    ax.axvline(F_REF, color="0.5", lw=0.9, ls=":", zorder=0)
    ax.text(F_REF * 1.04, -196.0, f"{case.cutoff_label} {F_REF / 1000:.0f} kHz",
            color="0.35", fontsize=9, va="bottom")

    text = (f"全频段 max |Hd−Ha|（线性）= {stats['max_linear']:.1e}\n"
            f"模拟 ≥ −20 dB：≤ {stats['db_20']:.3f} dB\n"
            f"模拟 ≥ −60 dB：≤ {stats['db_60']:.3f} dB\n"
            f"f_ref 处 Δ = {stats['ref_db']:+.3f} dB；"
            f"Nyquist ≈ A({stats['f_match'] / 1000:.2f} kHz)")
    ax.text(0.02, 0.04, text, transform=ax.transAxes, fontsize=8.8,
            va="bottom", ha="left",
            bbox=dict(boxstyle="round,pad=0.45", facecolor="white", alpha=0.85,
                      edgecolor="0.75"))

    ax.set_xlim(20.0, FS / 2.0)
    ax.set_ylim(-200.0, 5.0)
    ax.set_xlabel("频率 [Hz]")
    ax.set_ylabel("幅度 [dB]")
    ax.set_title(case.title)
    ax.grid(True, which="major", alpha=0.35)
    ax.grid(True, which="minor", alpha=0.15)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()

    OUT.mkdir(exist_ok=True)
    path = OUT / f"{case.slug}{config.suffix}.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def comparison(cases: tuple[Case, ...]) -> None:
    rows: list[tuple[str, dict, dict]] = []
    for config in COMPARISON:
        stats = {}
        for case in cases:
            _, st = measure(case, config)
            stats[case.slug] = st
        rows.append((config.label, stats["butterworth_N8"], stats["chebyshev2_N8"]))

    print("\n配置对照（maxlin = 全频段 max|Hd-Ha| 线性；dB60 = 模拟 >= -60 dB 内 max|ΔdB|；"
          "dRef = f_ref 处 Hd/Ha [dB]）")
    print(f"{'config':<24}{'butter maxlin':>15}{'butter dB60':>13}"
          f"{'cheb2 maxlin':>15}{'cheb2 dB60':>12}{'butter dRef':>12}{'cheb2 dRef':>12}")
    for label, butter, cheby2 in rows:
        print(f"{label:<24}{butter['max_linear']:>15.3e}{butter['db_60']:>13.3f}"
              f"{cheby2['max_linear']:>15.3e}{cheby2['db_60']:>12.3f}"
              f"{butter['ref_db']:>+12.4f}{cheby2['ref_db']:>+12.4f}")


def main() -> int:
    cases = make_cases()
    for case in cases:
        for config in CONFIGS:
            design, stats = measure(case, config)
            report(case, config, stats)
            describe_zeros(case, stats["zeros"])
            print(f"  图：{plot(case, config, design, stats)}")
    comparison(cases)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
