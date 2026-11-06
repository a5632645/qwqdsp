"""用三点功率匹配和迭代优化，把一阶模拟低通拟合为二阶数字 biquad。

模拟目标为 H_a(jw) = 1 / (1 + jw/wc)。分子、分母的功率响应都可写成
phi(w) = (1-s, s, 4s(1-s)) 的线性组合，s = sin^2(w/2)。
在 DC、wc 和 Nyquist 处匹配后，分子功率系数由分母 a1/a2 唯一确定，
剩下的两个自由度用稳定分母参数化，再最小化全带最大幅度误差（dB）。

对与其余 match_filter 实验相同的 5 个 fc 逐个做差分进化，输出
p/q、两个极点（分母根）、全带 max|ΔdB| 以及 5 条响应/误差曲线。
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import differential_evolution

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

FS = 48000.0
FC_LIST = (20.0, 117.6, 691.4, 4065.0, 23900.0)
OPTIMIZATION_GRID_SIZE = 4097
REPORT_GRID_SIZE = 262145
REFLECTION_LIMIT = 1.0 - 1e-9
OUTPUT_DIR = Path(__file__).parent / "output"


def get_phi(w):
    s = np.sin(np.asarray(w, dtype=float) / 2.0) ** 2
    return np.stack((1.0 - s, s, 4.0 * (1.0 - s) * s), axis=-1)


def get_power_coefficients(a1: float, a2: float) -> np.ndarray:
    return np.array([(1.0 + a1 + a2) ** 2,
                     (1.0 - a1 + a2) ** 2,
                     -4.0 * a2])


def fit_numerator(a1: float, a2: float, wc: float):
    """由 DC、wc、Nyquist 三个目标功率解出可实现的二阶分子。"""
    denominator = get_power_coefficients(a1, a2)
    points = np.array([0.0, wc, np.pi])
    point_phi = get_phi(points)
    target_power = 1.0 / (1.0 + (points / wc) ** 2)

    B0 = denominator[0] * target_power[0]
    B1 = denominator[1] * target_power[2]
    cutoff_power = target_power[1] * float(np.dot(point_phi[1], denominator))
    B2 = (cutoff_power - B0 * point_phi[1, 0] - B1 * point_phi[1, 1]) / point_phi[1, 2]

    if B0 <= 0.0 or B1 <= 0.0:
        return None
    root_B0 = np.sqrt(B0)
    root_B1 = np.sqrt(B1)
    W = (root_B0 + root_B1) / 2.0
    radicand = W * W + B2
    if radicand <= 0.0:
        return None

    b0 = (W + np.sqrt(radicand)) / 2.0
    b1 = (root_B0 - root_B1) / 2.0
    b2 = -B2 / (4.0 * b0)
    return np.array([b0, b1, b2]), np.array([B0, B1, B2])


def denominator_from_reflection(k1: float, k2: float) -> tuple[float, float]:
    """用 |k1|, |k2| < 1 参数化稳定的二阶分母。"""
    return k1 * (1.0 + k2), k2


def maximum_error_db(reflection: np.ndarray, w: np.ndarray, wc: float) -> float:
    a1, a2 = denominator_from_reflection(*reflection)
    fitted = fit_numerator(a1, a2, wc)
    if fitted is None:
        return 1e6

    _, numerator = fitted
    denominator = get_power_coefficients(a1, a2)
    phi = get_phi(w)
    digital_power = phi @ numerator / (phi @ denominator)
    analog_power = 1.0 / (1.0 + (w / wc) ** 2)
    if np.any(digital_power <= 0.0) or not np.all(np.isfinite(digital_power)):
        return 1e6
    error_db = 10.0 * np.log10(digital_power / analog_power)
    return float(np.max(np.abs(error_db)))


def design(fc: float, fs: float = FS):
    """对单个 fc 做差分进化、三点匹配与密集核验。"""
    wc = 2.0 * np.pi * fc / fs
    optimization_grid = np.unique(np.concatenate((
        np.linspace(0.0, np.pi, OPTIMIZATION_GRID_SIZE),
        np.array([wc]),
    )))
    result = differential_evolution(
        maximum_error_db,
        bounds=((-REFLECTION_LIMIT, REFLECTION_LIMIT),) * 2,
        args=(optimization_grid, wc),
        seed=42,
        tol=1e-10,
        atol=1e-10,
        maxiter=500,
        polish=True,
        workers=1,
    )

    k1, k2 = result.x
    a1, a2 = denominator_from_reflection(k1, k2)
    fitted = fit_numerator(a1, a2, wc)
    if fitted is None:
        raise RuntimeError(f"fc={fc:g}: 优化结果没有可实现的实数二阶分子")
    numerator, numerator_power = fitted
    denominator_power = get_power_coefficients(a1, a2)

    report_grid = np.unique(np.concatenate((
        np.linspace(0.0, np.pi, REPORT_GRID_SIZE),
        np.array([wc]),
    )))
    report_phi = get_phi(report_grid)
    report_digital_power = report_phi @ numerator_power / (report_phi @ denominator_power)
    report_analog_power = 1.0 / (1.0 + (report_grid / wc) ** 2)
    report_error_db = 10.0 * np.log10(report_digital_power / report_analog_power)
    report_error = float(np.max(np.abs(report_error_db)))

    match_points = np.array([0.0, wc, np.pi])
    match_phi = get_phi(match_points)
    digital_power = match_phi @ numerator_power / (match_phi @ denominator_power)
    analog_power = 1.0 / (1.0 + (match_points / wc) ** 2)
    poles = np.roots([1.0, a1, a2])
    if np.max(np.abs(poles)) >= 1.0:
        raise RuntimeError(f"fc={fc:g}: 优化分母不稳定")
    if not np.allclose(digital_power, analog_power, rtol=1e-10, atol=1e-14):
        raise RuntimeError(f"fc={fc:g}: 三点匹配未通过")

    return {
        "fc": fc,
        "wc": wc,
        "a1": float(a1),
        "a2": float(a2),
        "p": 1.0 + a1 + a2,
        "q": 1.0 - a1 + a2,
        "b": numerator,
        "beta": numerator_power,
        "poles": poles,
        "max_db": report_error,
        "wc_grid": report_grid,
        "error_db": report_error_db,
        "iterations": int(result.nit),
        "success": bool(result.success),
    }


def save_error_figure(rows, fs: float) -> Path:
    """保存 5 个 fc 的幅频响应和全带误差双面板图。"""
    OUTPUT_DIR.mkdir(exist_ok=True)
    fig, (response_ax, error_ax) = plt.subplots(1, 2, figsize=(12.5, 4.8))
    for index, row in enumerate(rows):
        color = f"C{index}"
        w = row["wc_grid"]
        f = w / np.pi * fs / 2.0
        phi = get_phi(w)
        denominator_power = get_power_coefficients(row["a1"], row["a2"])
        digital_power = (phi @ row["beta"]) / (phi @ denominator_power)
        analog_power = 1.0 / (1.0 + (w / row["wc"]) ** 2)
        response_ax.semilogx(f[1:], 10.0 * np.log10(analog_power[1:]),
                             color=color, lw=3.4, ls="--", alpha=0.30)
        response_ax.semilogx(f[1:], 10.0 * np.log10(digital_power[1:]),
                             color=color, lw=1.5, label=f"{row['fc']:g} Hz")
        error_ax.plot(w / np.pi, row["error_db"], color=color, lw=1.4,
                      label=f"{row['fc']:g} Hz")

    response_ax.set_xlim(10.0, fs / 2.0)
    response_ax.set_ylim(-80.0, 3.0)
    response_ax.set_xlabel("频率 / Hz")
    response_ax.set_ylabel("幅度 / dB")
    response_ax.set_title("幅频响应（实线：三点匹配 biquad；虚线：模拟原型）")
    response_ax.grid(True, which="both", alpha=0.25)
    response_ax.legend(loc="lower left", fontsize=8)

    error_limit = max(0.05, 1.05 * max(np.max(np.abs(row["error_db"])) for row in rows))
    error_ax.axhline(0.0, color="0.5", lw=0.8, ls=":")
    error_ax.set_xlim(0.0, 1.0)
    error_ax.set_ylim(-error_limit, error_limit)
    error_ax.set_xlabel("f / (fs/2)")
    error_ax.set_ylabel("ΔdB（数字 − 模拟）")
    error_ax.set_title("三点匹配后的全带误差")
    error_ax.grid(True, alpha=0.25)
    error_ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()

    path = OUTPUT_DIR / "lowpass_1st_error.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main() -> int:
    rows = []
    for fc in FC_LIST:
        try:
            rows.append(design(fc, FS))
        except RuntimeError as error:
            print(f"{fc:8.1f}  失败：{error}")
    if not rows:
        return 1

    print(f"fs = {FS:g} Hz，三点功率匹配 + 差分进化（反射系数），全带 max|ΔdB|")
    print(f"{'fc':>8}{'a1':>13}{'a2':>13}{'z1':>10}{'z2':>10}"
          f"{'p':>12}{'q':>12}{'maxdB':>9}{'iter':>5}")
    print("z1/z2 是分母 1+a1 z^-1+a2 z^-2 的两个根（极点）："
          "a1=-(z1+z2), a2=z1*z2")
    for row in rows:
        poles = np.sort(row["poles"])
        print(f"{row['fc']:8.1f}{row['a1']:13.9f}{row['a2']:13.9f}"
              f"{poles[0]:10.6f}{poles[1]:10.6f}"
              f"{row['p']:12.9f}{row['q']:12.9f}{row['max_db']:9.4f}"
              f"{row['iterations']:5d}")

    figure_path = save_error_figure(rows, FS)
    print(f"\n图：{figure_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
