"""二阶模拟高通匹配：保留 DC 二阶零点和 fc，放弃 Nyquist 条件。

模拟目标功率为 r^4 / ((1-r^2)^2 + (r/Q)^2)，r = w/wc。数字分子固定为
k(1-z^-1)^2；其增益由 fc 处匹配确定，分母仍有 p/q 两个稳定自由度，按全带
max|Delta dB| 优化。差分进化以冲击响应不变分母作为显式基线候选。
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

PI = np.pi
FS = 48000.0
FC_LIST = (20.0, 117.6, 691.4, 4065.0, 23900.0)
Q_LIST = (0.5, 1.0 / np.sqrt(2.0), 1.0, 2.0)
OPTIMIZATION_GRID_SIZE = 8193
REPORT_GRID_SIZE = 262145
OUTPUT_DIR = Path(__file__).parent / "output"


def phi_from_s(s):
    s = np.asarray(s, dtype=float)
    return np.stack((1.0 - s, s, 4.0 * s * (1.0 - s)), axis=-1)


def denominator_power(p: float, q: float) -> np.ndarray:
    return np.array([p * p, q * q, 4.0 - 2.0 * (p + q)])


def analog_power(w, wc: float, Q: float) -> np.ndarray:
    r = np.asarray(w, dtype=float) / wc
    return r ** 4 / ((1.0 - r * r) ** 2 + (r / Q) ** 2)


def impulse_invariant_denominator(wc: float, Q: float) -> tuple[float, float]:
    zeta = 1.0 / (2.0 * Q)
    decay = np.exp(-zeta * wc)
    if zeta <= 1.0:
        a1 = -2.0 * decay * np.cos(np.sqrt(1.0 - zeta * zeta) * wc)
    else:
        a1 = -2.0 * decay * np.cosh(np.sqrt(zeta * zeta - 1.0) * wc)
    return float(a1), float(decay * decay)


def _encode_denominator(p: float, q: float) -> np.ndarray:
    p_ratio = np.clip(p / 4.0, 1e-14, 1.0 - 1e-14)
    q_ratio = np.clip(q / (4.0 - p), 1e-14, 1.0 - 1e-14)
    return np.log([p_ratio / (1.0 - p_ratio), q_ratio / (1.0 - q_ratio)])


def _decode_denominator(x: np.ndarray) -> tuple[float, float]:
    sigmoid = 1.0 / (1.0 + np.exp(-np.clip(x, -36.0, 36.0)))
    p = 4.0 * sigmoid[0]
    q = (4.0 - p) * sigmoid[1]
    return float(p), float(q)


def design_from_pq(fc: float, fs: float, Q: float, p: float, q: float):
    """给定稳定分母，以 fc 定标双零点分子；不匹配 Nyquist。"""
    if not (p > 0.0 and q > 0.0 and p + q < 4.0):
        return None
    wc = 2.0 * PI * fc / fs
    alpha = denominator_power(p, q)
    sc = np.sin(wc / 2.0) ** 2
    phi_fc = np.array([1.0 - sc, sc, 4.0 * (1.0 - sc) * sc])
    denominator_fc = float(phi_fc @ alpha)
    if denominator_fc <= 0.0:
        return None

    # |k(1-z^-1)^2|^2 = k^2 * 16 sin^4(w/2).
    k_squared = Q * Q * denominator_fc / (16.0 * sc * sc)
    if k_squared <= 0.0 or not np.isfinite(k_squared):
        return None
    beta = np.array([0.0, 16.0 * k_squared, -4.0 * k_squared])
    k = float(np.sqrt(k_squared))
    b = np.array([k, -2.0 * k, k])
    a1 = (p - q) / 2.0
    a2 = (p + q) / 2.0 - 1.0
    poles = np.roots([1.0, a1, a2])
    if np.max(np.abs(poles)) >= 1.0:
        return None
    return {"p": p, "q": q, "a1": float(a1), "a2": float(a2),
            "alpha": alpha, "beta": beta, "b": b, "k": k,
            "poles": poles, "wc": wc, "fc": fc, "Q": Q}


def error_db(design, s: np.ndarray):
    """计算全带误差，DC 用模拟/数字四阶功率零点的解析比值极限。"""
    w = 2.0 * np.arcsin(np.sqrt(s))
    ph = phi_from_s(s)
    digital_power = (ph @ design["beta"]) / (ph @ design["alpha"])
    target = analog_power(w, design["wc"], design["Q"])
    nonzero = w > 0.0
    if np.any(digital_power[nonzero] <= 0.0):
        return None

    ratio = np.empty_like(s)
    ratio[nonzero] = digital_power[nonzero] / target[nonzero]
    ratio[~nonzero] = design["k"] ** 2 * design["wc"] ** 4 / design["p"] ** 2
    if np.any(ratio <= 0.0) or not np.all(np.isfinite(ratio)):
        return None
    return 10.0 * np.log10(ratio)


def _impulse_invariant_seed(fc: float, fs: float, Q: float) -> np.ndarray:
    wc = 2.0 * PI * fc / fs
    a1, a2 = impulse_invariant_denominator(wc, Q)
    p = 1.0 + a1 + a2
    q = 1.0 - a1 + a2
    return _encode_denominator(p, q)


def design(fc: float, fs: float, Q: float):
    s_grid = np.linspace(0.0, 1.0, OPTIMIZATION_GRID_SIZE)
    baseline_x = _impulse_invariant_seed(fc, fs, Q)

    def objective(x):
        p, q = _decode_denominator(x)
        candidate = design_from_pq(fc, fs, Q, p, q)
        if candidate is None:
            return 1e6
        error = error_db(candidate, s_grid)
        return 1e6 if error is None else float(np.max(np.abs(error)))

    population_size = 40
    rng = np.random.default_rng(42)
    initial_population = rng.uniform(-16.0, 16.0,
                                     size=(population_size, 2))
    initial_population[0] = baseline_x
    result = differential_evolution(
        objective,
        bounds=((-24.0, 24.0), (-24.0, 24.0)),
        init=initial_population,
        seed=42,
        tol=1e-9,
        atol=1e-11,
        maxiter=1200,
        polish=True,
        workers=1,
    )

    p, q = _decode_denominator(result.x)
    candidate = design_from_pq(fc, fs, Q, p, q)
    baseline_p, baseline_q = _decode_denominator(baseline_x)
    baseline = design_from_pq(fc, fs, Q, baseline_p, baseline_q)
    if candidate is None or baseline is None:
        raise RuntimeError(f"Q={Q:g}, fc={fc:g}: 优化结果或冲击不变基线不可行")
    baseline_error = error_db(baseline, s_grid)
    candidate_error = error_db(candidate, s_grid)
    baseline_max = float(np.max(np.abs(baseline_error)))
    candidate_max = float(np.max(np.abs(candidate_error)))
    if baseline_max < candidate_max:
        candidate = baseline
        candidate_error = baseline_error
        candidate_max = baseline_max

    report_s = np.linspace(0.0, 1.0, REPORT_GRID_SIZE)
    report_error = error_db(candidate, report_s)
    baseline_report_error = error_db(baseline, report_s)
    match_points = np.array([0.0, candidate["wc"]])
    match_phi = phi_from_s(np.sin(match_points / 2.0) ** 2)
    digital_match = ((match_phi @ candidate["beta"]) /
                     (match_phi @ candidate["alpha"]))
    target_match = analog_power(match_points, candidate["wc"], Q)
    if not np.allclose(digital_match, target_match, rtol=2e-10, atol=1e-13):
        raise RuntimeError(f"Q={Q:g}, fc={fc:g}: DC/fc 匹配校验失败")
    z_fc = np.exp(-1j * candidate["wc"])
    response_fc = ((candidate["b"][0] + candidate["b"][1] * z_fc +
                    candidate["b"][2] * z_fc * z_fc) /
                   (1.0 + candidate["a1"] * z_fc + candidate["a2"] * z_fc * z_fc))
    if not np.isclose(abs(response_fc) ** 2, Q * Q, rtol=2e-10, atol=1e-13):
        raise RuntimeError(f"Q={Q:g}, fc={fc:g}: biquad fc 增益校验失败")

    candidate.update({"s": report_s, "w": 2.0 * np.arcsin(np.sqrt(report_s)),
                      "error_db": report_error,
                      "baseline_max_db": float(np.max(np.abs(baseline_report_error))),
                      "max_db": float(np.max(np.abs(report_error))),
                      "success": bool(result.success),
                      "iterations": int(result.nit),
                      "matched_power": digital_match,
                      "target_power": target_match})
    return candidate


def save_figure(rows, fs: float) -> Path:
    fig = plt.figure(figsize=(12.5, 8.8))
    outer = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.22)
    frequency = np.geomspace(10.0, fs / 2.0, 4001)
    w_plot = 2.0 * PI * frequency / fs
    ph_plot = phi_from_s(np.sin(w_plot / 2.0) ** 2)

    for index, Q in enumerate(Q_LIST):
        cell = outer[index // 2, index % 2].subgridspec(
            2, 1, height_ratios=(2.1, 1.1), hspace=0.25)
        response_ax = fig.add_subplot(cell[0])
        error_ax = fig.add_subplot(cell[1])
        for row_index, row in enumerate(r for r in rows if r["Q"] == Q):
            color = f"C{row_index}"
            analog = analog_power(w_plot, row["wc"], Q)
            digital = ((ph_plot @ row["beta"]) /
                       (ph_plot @ row["alpha"]))
            response_ax.semilogx(frequency, 10.0 * np.log10(np.maximum(analog, 1e-300)),
                                 color=color, lw=3.2, ls="--", alpha=0.3)
            response_ax.semilogx(frequency, 10.0 * np.log10(np.maximum(digital, 1e-300)),
                                 color=color, lw=1.4, label=f"{row['fc']:g} Hz")
            error_ax.plot(row["w"] / PI, row["error_db"], color=color, lw=1.2,
                          label=f"{row['fc']:g} Hz ({row['max_db']:.3g} dB)")
        response_ax.set_title(f"Q = {Q:g}")
        response_ax.set_ylabel("幅度 / dB")
        response_ax.set_xlim(10.0, fs / 2.0)
        response_ax.set_ylim(-100.0, 15.0)
        response_ax.grid(True, which="both", alpha=0.25)
        response_ax.legend(loc="upper left", fontsize=7)
        response_ax.set_xlabel("频率 / Hz")
        error_ax.set_xlim(0.0, 1.0)
        error_ax.set_ylabel("ΔdB（数字 − 模拟）")
        error_ax.set_xlabel("f / (fs/2)")
        error_ax.grid(True, alpha=0.25)
        error_ax.legend(fontsize=6, ncol=2)

    fig.suptitle("二阶模拟高通：DC/fc 匹配，Nyquist 留作优化")
    path = OUTPUT_DIR / "highpass_2nd_response.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main() -> int:
    rows = []
    print(f"fs = {FS:g} Hz；匹配 DC/fc，放弃 Nyquist；冲击不变极点作为 DE 基线")
    print(f"{'Q':>8}{'fc':>8}{'p':>13}{'q':>13}{'z1':>11}{'z2':>11}"
          f"{'DE maxdB':>12}{'II maxdB':>12}{'iter':>7}{'status':>12}")
    for Q in Q_LIST:
        for fc in FC_LIST:
            row = design(fc, FS, Q)
            rows.append(row)
            poles = np.sort(row["poles"])
            status = "converged" if row["success"] else "limit"
            print(f"{Q:8.4f}{fc:8.1f}{row['p']:13.7g}{row['q']:13.7g}"
                  f"{poles[0].real:11.6f}{poles[1].real:11.6f}"
                  f"{row['max_db']:12.5g}{row['baseline_max_db']:12.5g}"
                  f"{row['iterations']:7d}{status:>12}")
    path = save_figure(rows, FS)
    print(f"\n图：{path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
