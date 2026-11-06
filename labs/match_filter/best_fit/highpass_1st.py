"""用三点功率匹配和差分进化，把一阶模拟高通拟合为二阶数字 biquad。

模拟目标为 |Ha(jw)|^2 = (w/wc)^2 / (1 + (w/wc)^2)。在 DC、wc 和 Nyquist
处匹配后，分子功率系数由分母 a1/a2 唯一确定，剩下的两个自由度用稳定分母
参数化，再最小化全带最大幅度误差（dB）。DC 处两边功率都为零，误差取 w->0
的解析极限，避免把 0/0 误差误记成零。
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
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


def analog_power(w, wc: float) -> np.ndarray:
    ratio = np.asarray(w, dtype=float) / wc
    return ratio * ratio / (1.0 + ratio * ratio)


def fit_numerator(a1: float, a2: float, wc: float):
    """由 DC、wc、Nyquist 三个目标功率解出可实现的二阶分子。"""
    denominator = get_power_coefficients(a1, a2)
    points = np.array([0.0, wc, np.pi])
    point_phi = get_phi(points)
    target_power = analog_power(points, wc)

    B0 = denominator[0] * target_power[0]
    B1 = denominator[1] * target_power[2]
    cutoff_power = target_power[1] * float(np.dot(point_phi[1], denominator))
    B2 = (cutoff_power - B0 * point_phi[1, 0] - B1 * point_phi[1, 1]) / point_phi[1, 2]

    if B0 < 0.0 or B1 <= 0.0:
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


def error_db(w: np.ndarray, wc: float, numerator: np.ndarray,
             denominator: np.ndarray) -> np.ndarray | None:
    """返回相对幅度误差；DC 点使用两边高通斜率的比值极限。"""
    phi = get_phi(w)
    denominator_value = phi @ denominator
    if np.any(denominator_value <= 0.0) or not np.all(np.isfinite(denominator_value)):
        return None

    digital_power = phi @ numerator
    nonzero = w > 0.0
    if np.any(digital_power[nonzero] <= 0.0):
        return None
    analog = analog_power(w[nonzero], wc)
    ratio = np.empty_like(w)
    ratio[nonzero] = digital_power[nonzero] / (denominator_value[nonzero] * analog)

    digital_slope_power = -numerator[0] + numerator[1] + 4.0 * numerator[2]
    low_frequency_ratio = (digital_slope_power * wc * wc /
                           (4.0 * denominator[0]))
    if not np.isfinite(low_frequency_ratio) or low_frequency_ratio <= 0.0:
        return None
    ratio[~nonzero] = low_frequency_ratio
    if np.any(ratio <= 0.0) or not np.all(np.isfinite(ratio)):
        return None
    return 10.0 * np.log10(ratio)


def maximum_error_db(reflection: np.ndarray, w: np.ndarray, wc: float) -> float:
    a1, a2 = denominator_from_reflection(*reflection)
    fitted = fit_numerator(a1, a2, wc)
    if fitted is None:
        return 1e6

    _, numerator = fitted
    error = error_db(w, wc, numerator, get_power_coefficients(a1, a2))
    if error is None:
        return 1e6
    return float(np.max(np.abs(error)))


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
        maxiter=1000,
        polish=True,
        workers=1,
    )

    a1, a2 = denominator_from_reflection(*result.x)
    fitted = fit_numerator(a1, a2, wc)
    if fitted is None:
        raise RuntimeError(f"fc={fc:g}: 优化结果没有可实现的实数二阶分子")
    numerator, numerator_power = fitted
    denominator_power = get_power_coefficients(a1, a2)

    report_grid = np.unique(np.concatenate((
        np.linspace(0.0, np.pi, REPORT_GRID_SIZE),
        np.array([wc]),
    )))
    report_error_db = error_db(report_grid, wc, numerator_power, denominator_power)
    if report_error_db is None:
        raise RuntimeError(f"fc={fc:g}: 密集频率网格上出现不可行功率响应")
    poles = np.roots([1.0, a1, a2])

    match_points = np.array([0.0, wc, np.pi])
    point_phi = get_phi(match_points)
    digital_power = (point_phi @ numerator_power) / (point_phi @ denominator_power)
    target_power = analog_power(match_points, wc)
    z = np.exp(-1j * match_points)
    response = ((numerator[0] + numerator[1] * z + numerator[2] * z * z) /
                (1.0 + a1 * z + a2 * z * z))
    if np.max(np.abs(poles)) >= 1.0:
        raise RuntimeError(f"fc={fc:g}: 优化分母不稳定")
    if not np.allclose(digital_power, target_power, rtol=1e-10, atol=1e-14):
        raise RuntimeError(f"fc={fc:g}: 三点匹配未通过")
    if not np.allclose(np.abs(response) ** 2, target_power, rtol=1e-10, atol=1e-14):
        raise RuntimeError(f"fc={fc:g}: 直接传递函数的三点校验未通过")

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
        "max_db": float(np.max(np.abs(report_error_db))),
        "grid": report_grid,
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
        w = row["grid"]
        frequency = w / np.pi * fs / 2.0
        phi = get_phi(w)
        denominator = get_power_coefficients(row["a1"], row["a2"])
        digital_power = (phi @ row["beta"]) / (phi @ denominator)
        analog = analog_power(w, row["wc"])
        response_ax.semilogx(frequency[1:], 10.0 * np.log10(analog[1:]),
                             color=color, lw=3.4, ls="--", alpha=0.30)
        response_ax.semilogx(frequency[1:], 10.0 * np.log10(digital_power[1:]),
                             color=color, lw=1.5, label=f"{row['fc']:g} Hz")
        error_ax.plot(w / np.pi, row["error_db"], color=color, lw=1.4,
                      label=f"{row['fc']:g} Hz ({row['max_db']:.3g} dB)")

    response_ax.set_xlim(10.0, fs / 2.0)
    response_ax.set_ylim(-80.0, 3.0)
    response_ax.set_xlabel("频率 / Hz")
    response_ax.set_ylabel("幅度 / dB")
    response_ax.set_title("幅频响应（实线：三点匹配 biquad；虚线：模拟原型）")
    response_ax.grid(True, which="both", alpha=0.25)
    response_ax.legend(loc="lower right", fontsize=8)

    max_error = max(row["max_db"] for row in rows)
    error_ax.axhline(0.0, color="0.5", lw=0.8, ls=":")
    error_ax.set_xlim(0.0, 1.0)
    error_limit = max(0.3, 1.1 * max_error)
    error_ax.set_ylim(-error_limit, error_limit)
    error_ax.set_yscale("symlog", linthresh=1e-7)
    tick_values = np.array([-0.3, -0.1, -0.01, -0.001, -1e-4, -1e-5, -1e-6,
                            -1e-7, 0.0, 1e-7, 1e-6, 1e-5, 1e-4, 0.001,
                            0.01, 0.1, 0.3])
    error_ax.set_yticks(tick_values[np.abs(tick_values) <= error_limit])
    error_ax.yaxis.set_major_formatter(
        FuncFormatter(lambda value, _: "0" if value == 0.0 else
                      f"{value:.0e}" if abs(value) < 0.01 else f"{value:g}"))
    error_ax.set_xlabel("f / (fs/2)")
    error_ax.set_ylabel("ΔdB（数字 − 模拟；symlog）")
    error_ax.set_title("三点匹配后的全带误差（颜色同左图）")
    error_ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()

    path = OUTPUT_DIR / "highpass_1st_error.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main() -> int:
    rows = []
    for fc in FC_LIST:
        rows.append(design(fc, FS))

    print(f"fs = {FS:g} Hz，三点功率匹配 + 差分进化（反射系数），全带 max|ΔdB|")
    print("高通在 DC 的误差按 w→0 的比值极限计算。")
    print(f"{'fc':>8}{'a1':>13}{'a2':>13}{'z1':>10}{'z2':>10}"
          f"{'p':>12}{'q':>12}{'maxdB':>12}{'iter':>6}{'状态':>10}")
    for row in rows:
        poles = np.sort(row["poles"])
        print(f"{row['fc']:8.1f}{row['a1']:13.9f}{row['a2']:13.9f}"
              f"{poles[0]:10.6f}{poles[1]:10.6f}"
              f"{row['p']:12.9f}{row['q']:12.9f}{row['max_db']:12.4g}"
              f"{row['iterations']:6d}{'收敛' if row['success'] else '达到上限':>10}")

    figure_path = save_error_figure(rows, FS)
    print(f"\n图：{figure_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
