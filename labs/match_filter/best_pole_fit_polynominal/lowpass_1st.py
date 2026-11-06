"""用归一化截止频率多项式近似一阶低通的 minimax 极点轨迹。

离线参照值来自 best_pole_fit 的逐 fc 最佳分母。对两个实极点分别做离散
minimax 多项式拟合，再由极点恢复二阶分母并用相同响应采样网格核验误差。
"""
from __future__ import annotations

import sys
import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
BEST_POLE_FIT_DIR = HERE.parent / "best_pole_fit"
sys.path.insert(0, str(BEST_POLE_FIT_DIR))

_REFERENCE_SPEC = importlib.util.spec_from_file_location(
    "best_pole_fit_lowpass_1st", BEST_POLE_FIT_DIR / "lowpass_1st.py")
if _REFERENCE_SPEC is None or _REFERENCE_SPEC.loader is None:
    raise ImportError("无法加载 best_pole_fit 的 lowpass_1st 参考实现")
reference = importlib.util.module_from_spec(_REFERENCE_SPEC)
_REFERENCE_SPEC.loader.exec_module(reference)

FS_REFERENCE = 48000.0
FC_REFERENCE = np.array([20.0, 117.6, 691.4, 4065.0, 23900.0])
X_MIN = float(FC_REFERENCE[0] / (FS_REFERENCE / 2.0))
X_MAX = float(FC_REFERENCE[-1] / (FS_REFERENCE / 2.0))

# Change these values to regenerate a piecewise polynomial with a different runtime cost.
PIECE_COUNT = 32
POLYNOMIAL_DEGREE = 7
FIT_CUTOFF_GRID_SIZE = 9
VALIDATION_CUTOFF_GRID_SIZE = 129
FIT_RESPONSE_GRID_SIZE = 4001
RESPONSE_GRID_SIZE = 400001
OUTPUT_DIR = HERE / "output"
REFERENCE_CACHE = OUTPUT_DIR / "lowpass_1st_reference.npz"


def design_reference(x: float):
    """在归一化 Nyquist 单位下计算逐点 minimax 设计。"""
    p, q = reference.best_denominator(float(x), 2.0)
    match = reference.build_match(float(x), 2.0, p, q)
    if match is None:
        raise RuntimeError(f"x={x:.9g}: minimax 分母不可行")
    poles = np.sort(np.real_if_close(match.poles()).real)
    return np.array([p, q]), poles, match


def fit_feasible_polynomials(x: np.ndarray, params: np.ndarray, degree: int,
                             anchor_x: np.ndarray, anchor_params: np.ndarray):
    """在每个频率段内拟合逐点 minimax 的 p/q 参数。"""
    root_edges = np.linspace(np.sqrt(X_MIN), np.sqrt(X_MAX), PIECE_COUNT + 1)
    edges = root_edges * root_edges
    edges[0] = X_MIN
    edges[-1] = X_MAX
    coefficients = np.empty((2, PIECE_COUNT, degree + 1))
    coefficients.fill(0.0)
    for segment in range(PIECE_COUNT):
        left, right = edges[segment:segment + 2]
        local_indices = np.flatnonzero((x >= left) & (x <= right))
        local_degree = min(degree, len(local_indices) - 1)
        local_x = x[local_indices]
        local_t = (2.0 * local_x - left - right) / (right - left)
        for parameter in range(2):
            local_coefficients = np.polynomial.chebyshev.chebfit(
                local_t, params[local_indices, parameter], local_degree)
            desired_endpoints = np.interp([left, right], anchor_x,
                                          anchor_params[:, parameter])
            actual_endpoints = np.polynomial.chebyshev.chebval(
                np.array([-1.0, 1.0]), local_coefficients)
            delta_left, delta_right = desired_endpoints - actual_endpoints
            local_coefficients[0] += 0.5 * (delta_left + delta_right)
            if local_degree >= 1:
                local_coefficients[1] += 0.5 * (delta_right - delta_left)
            coefficients[parameter, segment, :len(local_coefficients)] = local_coefficients
    return edges, coefficients


def eval_piecewise_polynomial(edges, coefficients, parameter: int, x):
    x = np.asarray(x, dtype=float)
    flat_x = x.ravel()
    segment = np.clip(np.searchsorted(edges, flat_x, side="right") - 1,
                      0, len(edges) - 2)
    left = edges[segment]
    right = edges[segment + 1]
    local_t = (2.0 * flat_x - left - right) / (right - left)
    result = np.array([
        np.polynomial.chebyshev.chebval(t, coefficients[parameter, index])
        for index, t in zip(segment, local_t)
    ])
    return result.reshape(x.shape)


def response_error(match, x: float, n: int):
    return reference.error_db(match, x, 2.0, n=n)[0]


def build_reference_data(x_grid: np.ndarray):
    denominator_params = np.empty((len(x_grid), 2))
    poles = np.empty((len(x_grid), 2))
    errors = np.empty(len(x_grid))
    matches = []
    for index, x in enumerate(x_grid):
        denominator_params[index], poles[index], match = design_reference(float(x))
        errors[index] = np.max(np.abs(response_error(match, float(x), FIT_RESPONSE_GRID_SIZE)))
        matches.append(match)
    return denominator_params, poles, errors, matches


def fit_and_validate():
    # Place reference designs inside every segment so local polynomials do not extrapolate.
    root_edges = np.linspace(np.sqrt(X_MIN), np.sqrt(X_MAX), PIECE_COUNT + 1)
    local_theta = np.linspace(0.0, np.pi, FIT_CUTOFF_GRID_SIZE)
    fit_nodes = []
    for left_root, right_root in zip(root_edges[:-1], root_edges[1:]):
        left, right = left_root * left_root, right_root * right_root
        fit_nodes.extend((left + right) / 2.0 + (right - left) / 2.0 *
                         np.cos(local_theta))
    fit_x = np.unique(np.asarray(fit_nodes))
    validation_roots = np.linspace(np.sqrt(X_MIN), np.sqrt(X_MAX),
                                   VALIDATION_CUTOFF_GRID_SIZE)
    validation_x = validation_roots * validation_roots
    validation_x[0] = X_MIN
    validation_x[-1] = X_MAX
    segment_edges = root_edges * root_edges
    segment_edges[0] = X_MIN
    segment_edges[-1] = X_MAX
    cached = None
    if REFERENCE_CACHE.exists():
        with np.load(REFERENCE_CACHE) as archive:
            if (np.array_equal(archive["fit_x"], fit_x) and
                    np.array_equal(archive["validation_x"], validation_x) and
                    int(archive["response_grid_size"]) == FIT_RESPONSE_GRID_SIZE):
                cached = {name: archive[name] for name in archive.files}

    if cached is None:
        fit_params, fit_poles, fit_optimal_error, _ = build_reference_data(fit_x)
        validation_params, validation_poles, _, _ = build_reference_data(validation_x)
        OUTPUT_DIR.mkdir(exist_ok=True)
        np.savez(REFERENCE_CACHE, fit_x=fit_x, fit_params=fit_params,
                 fit_poles=fit_poles, fit_optimal_error=fit_optimal_error,
                 validation_x=validation_x, validation_params=validation_params,
                 validation_poles=validation_poles,
                 response_grid_size=np.array(FIT_RESPONSE_GRID_SIZE))
    else:
        fit_params = cached["fit_params"]
        fit_poles = cached["fit_poles"]
        fit_optimal_error = cached["fit_optimal_error"]
        validation_params = cached["validation_params"]
        validation_poles = cached["validation_poles"]

    edges, coefficients = fit_feasible_polynomials(
        fit_x, fit_params, POLYNOMIAL_DEGREE, validation_x, validation_params)
    polynomial_params = np.column_stack([
        eval_piecewise_polynomial(edges, coefficients, i, validation_x) for i in range(2)
    ])
    polynomial_error = np.empty(len(validation_x))
    for index, x in enumerate(validation_x):
        p, q = polynomial_params[index]
        match = reference.build_match(float(x), 2.0, float(p), float(q))
        if match is None:
            polynomial_error[index] = np.inf
            continue
        polynomial_error[index] = np.max(np.abs(response_error(match, float(x),
                                                               RESPONSE_GRID_SIZE)))
    polynomial_poles = np.array([
        np.sort(np.real_if_close(np.roots([1.0, (p - q) / 2.0,
                                           (p + q) / 2.0 - 1.0])).real)
        for p, q in polynomial_params
    ])

    validation_optimal_error = np.empty(len(validation_x))
    for index, (x, (p, q)) in enumerate(zip(validation_x, validation_params)):
        match = reference.build_match(float(x), 2.0, float(p), float(q))
        validation_optimal_error[index] = np.max(np.abs(
            response_error(match, float(x), RESPONSE_GRID_SIZE)))

    boundary_errors = []
    boundary_param_jumps = []
    for segment, x in enumerate(segment_edges[1:-1], start=1):
        left_t = 1.0
        right_t = -1.0
        left_params = np.array([
            np.polynomial.chebyshev.chebval(left_t, coefficients[i, segment - 1])
            for i in range(2)
        ])
        right_params = np.array([
            np.polynomial.chebyshev.chebval(right_t, coefficients[i, segment])
            for i in range(2)
        ])
        boundary_param_jumps.append(float(np.max(np.abs(left_params - right_params))))
        for params_at_edge in (left_params, right_params):
            match = reference.build_match(float(x), 2.0,
                                          float(params_at_edge[0]),
                                          float(params_at_edge[1]))
            if match is None:
                boundary_errors.append(np.inf)
            else:
                boundary_errors.append(float(np.max(np.abs(
                    response_error(match, float(x), RESPONSE_GRID_SIZE)))))

    return {
        "fit_x": fit_x,
        "fit_poles": fit_poles,
        "fit_optimal_error": fit_optimal_error,
        "validation_x": validation_x,
        "validation_params": validation_params,
        "validation_poles": validation_poles,
        "validation_optimal_error": validation_optimal_error,
        "polynomial_poles": polynomial_poles,
        "polynomial_error": polynomial_error,
        "coefficients": coefficients,
        "edges": edges,
        "uniform_validation_count": len(validation_x),
        "boundary_max_error": float(np.max(boundary_errors)),
        "boundary_max_param_jump": float(np.max(boundary_param_jumps)),
    }


def save_figure(data) -> Path:
    x = data["validation_x"]
    poles = data["validation_poles"]
    polynomial_poles = data["polynomial_poles"]
    optimal = data["validation_optimal_error"]
    achieved = data["polynomial_error"]

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.0), sharex=True)
    axes[0].plot(x, poles[:, 0], "o", ms=2.5, color="C0", label="minimax pole 1")
    axes[0].plot(x, polynomial_poles[:, 0], color="C0", lw=1.5,
                 label=f"{PIECE_COUNT} x degree-{POLYNOMIAL_DEGREE} polynomial")
    axes[0].plot(x, poles[:, 1], "o", ms=2.5, color="C1", label="minimax pole 2")
    axes[0].plot(x, polynomial_poles[:, 1], color="C1", lw=1.5,
                 label=f"{PIECE_COUNT} x degree-{POLYNOMIAL_DEGREE} polynomial")
    axes[0].set_ylabel("pole location")
    axes[0].set_ylim(-1.02, 1.02)
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(ncol=2, fontsize=8)

    axes[1].plot(x, optimal, color="0.35", lw=1.3, label="pointwise minimax")
    axes[1].plot(x, achieved, color="C3", lw=1.3, label="polynomial poles")
    axes[1].set_xlabel("normalized cutoff x = 2 fc / fs")
    axes[1].set_ylabel("max |ΔdB|")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.suptitle("First-order lowpass: polynomial pole approximation")
    fig.tight_layout()
    OUTPUT_DIR.mkdir(exist_ok=True)
    path = OUTPUT_DIR / "lowpass_1st_polynomial.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def write_cpp_header(data) -> Path:
    edges = data["edges"]
    coefficients = data["coefficients"]

    def format_array(values, indent):
        rows = []
        for value in values:
            rows.append(" " * indent + "{" + ", ".join(f"{item:.17g}" for item in value) + "}")
        return ",\n".join(rows)

    header = f'''#pragma once
#include "qwqdsp/filter/biquad_coeff.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>

namespace qwqdsp_filter::best_fit_detail {{
inline constexpr double kLowpass1stXMin = {X_MIN:.17g};
inline constexpr double kLowpass1stXMax = {X_MAX:.17g};
inline constexpr size_t kLowpass1stPieceCount = {PIECE_COUNT};
inline constexpr size_t kLowpass1stPolynomialDegree = {POLYNOMIAL_DEGREE};
inline constexpr std::array<double, {len(edges)}> kLowpass1stEdges{{{{
    {", ".join(f"{value:.17g}" for value in edges)}
}}}};
inline constexpr std::array<std::array<double, {POLYNOMIAL_DEGREE + 1}>, {PIECE_COUNT}> kLowpass1stPCoeffs{{{{
{format_array(coefficients[0], 4)}
}}}};
inline constexpr std::array<std::array<double, {POLYNOMIAL_DEGREE + 1}>, {PIECE_COUNT}> kLowpass1stQCoeffs{{{{
{format_array(coefficients[1], 4)}
}}}};

inline double EvalChebyshev(std::array<double, {POLYNOMIAL_DEGREE + 1}> const& c,
                            double x) noexcept {{
    double b1 = 0.0;
    double b2 = 0.0;
    for (size_t i = {POLYNOMIAL_DEGREE}; i > 0; --i) {{
        double const b0 = 2.0 * x * b1 - b2 + c[i];
        b2 = b1;
        b1 = b0;
    }}
    return x * b1 - b2 + c[0];
}}

inline BiquadCoeff LowpassOnepolePolynomial(float wc) noexcept {{
    double const x = std::clamp(static_cast<double>(wc) / 3.14159265358979323846,
                                kLowpass1stXMin, kLowpass1stXMax);
    auto const upper = std::upper_bound(kLowpass1stEdges.begin(),
                                        kLowpass1stEdges.end(), x);
    size_t const segment = std::min(
        static_cast<size_t>(std::distance(kLowpass1stEdges.begin(), upper) - 1),
        kLowpass1stPieceCount - 1);
    double const left = kLowpass1stEdges[segment];
    double const right = kLowpass1stEdges[segment + 1];
    double const t = (2.0 * x - left - right) / (right - left);
    double const p = EvalChebyshev(kLowpass1stPCoeffs[segment], t);
    double const q = EvalChebyshev(kLowpass1stQCoeffs[segment], t);

    double const a1 = (p - q) / 2.0;
    double const a2 = (p + q) / 2.0 - 1.0;
    double const sc = std::sin(static_cast<double>(wc) / 2.0);
    double const s = sc * sc;
    double const phi0 = 1.0 - s;
    double const phi1 = s;
    double const phi2 = 4.0 * phi0 * phi1;
    double const A0 = p * p;
    double const A1 = q * q;
    double const A2 = 4.0 - 2.0 * (p + q);
    double const T = 1.0 / (1.0 + (3.14159265358979323846 / wc) *
                                  (3.14159265358979323846 / wc));
    double const B0 = A0;
    double const B1 = T * A1;
    double const R = 0.5 * (phi0 * A0 + phi1 * A1 + phi2 * A2);
    double const B2 = (R - B0 * phi0 - B1 * phi1) / phi2;
    double const W = (std::sqrt(B0) + std::sqrt(B1)) / 2.0;
    double const b0 = (W + std::sqrt(std::max(W * W + B2, 0.0))) / 2.0;
    double const b1 = (std::sqrt(B0) - std::sqrt(B1)) / 2.0;
    double const b2 = -B2 / (4.0 * b0);
    return DoubleBiquadCoeff{{b0, b1, b2, a1, a2}}.ToFloat();
}}
}} // namespace qwqdsp_filter::best_fit_detail
'''
    path = HERE / "lowpass_1st_polynomial.hpp"
    path.write_text(header, encoding="utf-8")
    return path


def main() -> int:
    data = fit_and_validate()
    stable = bool(np.all(np.max(np.abs(data["polynomial_poles"]), axis=1) < 1.0))
    if (not stable or not np.all(np.isfinite(data["polynomial_error"])) or
            not np.isfinite(data["boundary_max_error"]) or
            np.max(data["polynomial_error"]) > 0.14):
        raise RuntimeError("多项式分母在验证频率上不稳定或分子不可实现")
    coefficients = data["coefficients"]
    excess = data["polynomial_error"] - data["validation_optimal_error"]
    worst_index = int(np.argmax(excess))

    print(f"normalized cutoff domain: [{X_MIN:.12g}, {X_MAX:.12g}]")
    print(f"fit degree={POLYNOMIAL_DEGREE}, fit points={FIT_CUTOFF_GRID_SIZE}, "
          f"validation points={VALIDATION_CUTOFF_GRID_SIZE}, "
          f"response grid={RESPONSE_GRID_SIZE}")
    print(f"piece count={PIECE_COUNT}, local degree={POLYNOMIAL_DEGREE}")
    print(f"max pointwise minimax error: "
          f"{np.max(data['validation_optimal_error']):.9g} dB")
    print(f"max polynomial response error: "
          f"{np.max(data['polynomial_error']):.9g} dB")
    print(f"max response-error excess: {np.max(excess):.9g} dB at "
          f"x={data['validation_x'][worst_index]:.9g}")
    print(f"max segment-boundary response error: "
          f"{data['boundary_max_error']:.9g} dB")
    print(f"max segment-boundary p/q jump: "
          f"{data['boundary_max_param_jump']:.9g}")
    print(f"stable polynomial poles: {stable}")
    print(f"figure: {save_figure(data)}")
    print(f"C++ header: {write_cpp_header(data)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
