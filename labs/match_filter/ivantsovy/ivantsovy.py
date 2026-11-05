"""Ivantsov 理想双线性 / 双二次变换（论文 §4 一般形式）。

论文：Yuriy Ivantsov, *On the Ideal Bilinear and Biquadratic Digital Filter*,
2026, https://ivantsovy.com/research/paper1.pdf

做法：把连续域传递函数写成

    H(S) = gain * S^m * prod(1 - S/z_i) / prod(1 - S/p_i),

即分子先分离出 ``m`` 个位于原点的零点，剩余有限零点 / 极点配对成常数项为 1
的实系数一次 / 二次因子，再逐因子映射到 z 域（论文 (4.2)–(4.6)）：

    a S + 1          -> (1 + c z^-1) / (1 + c)
    a S^2 + b S + 1  -> (1 + c1 z^-1 + c2 z^-2) / (1 + c1 + c2)
    S                -> (1 - z^-1) * omega / (2 pi)        (论文 (4.2) 第二行)
    S^m              -> 上述因子重复 m 次

其中

    c     = phi(omega a)          (a >= 0)，或 1 / phi(omega a)   (a < 0)
    c1,c2 = {phi1, phi2}(x, y)    (b >= 0)，或 {phi1/phi2, 1/phi2} (b < 0)
    x = omega sqrt(a),  y = b / (2 sqrt(a)),  omega = fs / f_ref
    phi(x)     = (pi - sqrt(x^2 + sigma^2)) / (pi + sqrt(x^2 + sigma^2))
    phi1,phi2  = 论文 (3.7)

分子的阶数亏欠（s -> infinity 的零点）按论文 (4.7) 退化为 a = 0 的一次因子：
每个亏欠零点映射为 (1 + phi(0) z^-1) / (1 + phi(0))。S^m 因子负责 LP/HP/BP/BS
中原点处的零点；s = 0 处的极点（积分器）没有对应规则，会被拒绝。

约定：连续域变量按 S = s / (2 pi f_ref) 归一化，f_ref 为原型参考频率（Hz），
模拟幅度在频率 f 处用 S = j f / f_ref 求值。系数 sigma 是自由参数；论文给出的
幅度最优区间是 [2, pi*sqrt(2/3)]。本 lab 默认取上端点，并对照 sigma = 2 +
预畸变 (:func:`prewarped_omega`)：后者让参考频率精确匹配，但 Nyquist 的精确
匹配点从 fs/sigma 移到 f_ref * omega / sigma。

数值上始终按因子保存，不展开成高阶 z 多项式——高阶展开系数病态，而因子串联
既稳定又天然给出二阶节。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

TOL = 1e-9
"""判定实根 / 共轭对 / 原点零点的相对阈值。"""


@dataclass(frozen=True)
class Factor:
    """一个映射后的 z 域因子 ``gain * poly(z^-1) / norm``。

    ``norm`` 是 ``poly`` 在 z = 1 处的值，用于把一次 / 二次因子的直流增益
    归一化为 1。原点零点因子 ``(1 - z^-1) * omega/(2*pi)`` 的 ``norm`` 为 1，
    ``gain`` 为 ``omega/(2*pi)``。
    """

    poly: np.ndarray
    norm: float = 1.0
    gain: float = 1.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "poly", np.asarray(self.poly, dtype=float))


@dataclass(frozen=True)
class Design:
    """映射后的 z 域传递函数：``gain * prod(num) / prod(den)``。"""

    num: tuple[Factor, ...]
    den: tuple[Factor, ...]
    gain: float = 1.0


def phi(x: float, sigma: float) -> float:
    """论文 (2.6)：``(pi - sqrt(x^2 + sigma^2)) / (pi + sqrt(x^2 + sigma^2))``。"""
    r = np.hypot(x, sigma)
    return float((np.pi - r) / (np.pi + r))


def phi12(x: float, y: float, sigma: float) -> tuple[float, float]:
    """论文 (3.7) 的 ``phi1(x, y)`` 与 ``phi2(x, y)``。

    两者共用分母 ``d``；``1 + phi1 + phi2 = 4 pi^2 / d``（论文结论部分）。

    ``|y| <= 1`` 时 ``nu + kappa`` 在解析上非负，且等于
    ``4 x^2 y^2 (1-y^2) / (nu - kappa)``。``b = 0``（纯虚轴零点，例如带阻）
    且 ``x > sigma`` 时它恰好为 0，而直接算法是 ``nu`` 量级的相减：浮点舍入
    会给出约 ``1e-16 * nu`` 的负值，此时只能截断。``|y| > 1`` 不在本映射的
    适用范围（实根应按 (4.7) 拆成一次因子），仍按原样报错。
    """
    x2 = x * x
    y2 = y * y
    s2 = sigma * sigma
    t = 2.0 * y2 - 1.0
    nu = np.sqrt(x2 * x2 + 2.0 * x2 * s2 * t + s2 * s2)
    kappa = x2 * t + s2
    if kappa >= 0.0:
        radicand = nu + kappa
    else:
        # nu + kappa = (nu^2 - kappa^2) / (nu - kappa)，分母是两个正数之和。
        tail = x2 * x2 * (1.0 - t * t)
        radicand = tail / (nu - kappa) if tail > 0.0 else 0.0
    if radicand < 0.0:
        if abs(y) > 1.0:
            raise ValueError("nu + kappa < 0：需要 |y| <= 1 的二次因子")
        radicand = 0.0
    root = np.sqrt(radicand)
    d = np.pi * np.pi + np.pi * np.sqrt(2.0) * root + nu
    p1 = (2.0 * np.pi * np.pi - 2.0 * nu) / d
    p2 = (np.pi * np.pi - np.pi * np.sqrt(2.0) * root + nu) / d
    return float(p1), float(p2)


def complete_conjugates(roots) -> np.ndarray:
    """补全共轭对；实根（含原点零点）不受影响。"""
    roots = [complex(r) for r in np.asarray(roots, dtype=complex).ravel()]
    if not any(r.imag < -TOL for r in roots):
        roots += [r.conjugate() for r in roots if r.imag > TOL]
    return np.asarray(roots, dtype=complex)


def split_origin(roots) -> tuple[int, np.ndarray]:
    """分离 s = 0 处的根；返回 ``(重数, 非零根)``。"""
    roots = complete_conjugates(roots)
    if roots.size == 0:
        return 0, roots
    scale = max(1.0, float(np.max(np.abs(roots))))
    finite = roots[np.abs(roots) > TOL * scale]
    return int(roots.size - finite.size), finite


def factorize(roots) -> list[tuple[str, float, float]]:
    """把非零根组成常数项为 1 的一次 / 二次因子。

    返回 ``("linear", a, 0)`` 或 ``("quadratic", a, b)``，对应
    ``1 + a S`` 与 ``1 + b S + a S^2``。共轭复根只取上半平面的一支。
    """
    factors: list[tuple[str, float, float]] = []
    _, finite = split_origin(roots)
    for r in finite:
        scale = max(1.0, abs(r))
        if abs(r.imag) <= TOL * scale:
            factors.append(("linear", float(-1.0 / r.real), 0.0))
        elif r.imag > 0.0:
            re = r.real if abs(r.real) > TOL * scale else 0.0
            a = 1.0 / (r.real * r.real + r.imag * r.imag)
            b = -2.0 * re * a
            factors.append(("quadratic", float(a), float(b)))
    return factors


def map_linear(a: float, omega: float, sigma: float) -> Factor:
    """映射一次因子 ``a S + 1``（论文 (4.2)/(4.3)）。"""
    if a >= 0.0:
        c = phi(omega * a, sigma)
    else:
        p = phi(omega * a, sigma)
        if p == 0.0:
            raise ValueError("phi(omega*a) = 0：该截止频率下映射出现除零")
        c = 1.0 / p
    poly = np.array([1.0, c])
    norm = float(poly.sum())
    if norm == 0.0:
        raise ValueError("映射因子在 z = 1 处为零，无法归一化")
    return Factor(poly, norm, 1.0)


def map_quad(a: float, b: float, omega: float, sigma: float) -> Factor:
    """映射二次因子 ``a S^2 + b S + 1``（论文 (4.4)/(4.5)/(4.6)）。"""
    if a <= 0.0:
        raise ValueError("二次因子的 S^2 系数必须为正")
    root_a = np.sqrt(a)
    x = omega * root_a
    y = b / (2.0 * root_a)
    p1, p2 = phi12(x, y, sigma)
    if b >= 0.0:
        c1, c2 = p1, p2
    else:
        if p2 == 0.0:
            raise ValueError("phi2 = 0：该截止频率下映射出现除零")
        c1, c2 = p1 / p2, 1.0 / p2
    poly = np.array([1.0, c1, c2])
    norm = float(poly.sum())
    if norm == 0.0:
        raise ValueError("映射因子在 z = 1 处为零，无法归一化")
    return Factor(poly, norm, 1.0)


def map_origin(omega: float) -> Factor:
    """映射原点零点 ``S -> (1 - z^-1) * omega/(2*pi)``（论文 (4.2)）。"""
    return Factor(np.array([1.0, -1.0]), 1.0, float(omega / (2.0 * np.pi)))


def prewarped_omega(fc: float, fs: float, sigma: float) -> float:
    """预畸变的 ``omega``：使数字参考频率精确落在 ``fc``（论文 (2.10)）。

    返回的是应代入 ``phi(omega a)`` 与 ``x = omega sqrt(a)`` 的 omega，而不是
    ``fs / fc``。它由 ``F(omega) = fc`` 解出：

        omega = sqrt(sigma^2 + pi^2 cot^2(pi / (fs/fc))).

    ``fc > fs/2`` 时按约定保持 ``omega = fs/fc``。
    """
    if fc > fs / 2.0:
        return fs / fc
    omega_raw = fs / fc
    ctg_val = 1.0 / np.tan(np.pi / omega_raw)
    return float(np.sqrt(sigma * sigma + (np.pi * ctg_val) ** 2))


def normalized_gain(zeros, poles, gain: float = 1.0) -> float:
    """把 zpk 形式的增益换算成“原点零点 + 常数项 1 因子”形式的增益。

    若 ``H(S) = gain * prod(S - z) / prod(S - p)``，则

        H(S) = gain_norm * S^m * prod(1 - S/z_f) / prod(1 - S/p_f),

    其中 ``gain_norm = gain * prod(-z_f) / prod(-p_f)``。scipy 的
    ``cheby2(..., output='zpk')`` 给出的 ``gain`` 必须经过这一步，才能直接
    作为 :func:`design` 的 ``gain``；否则 Chebyshev II 的通带 / 阻带电平会
    整体偏掉（例如低通的 ``gain = 1e-4``）。
    """
    _, finite_zeros = split_origin(zeros)
    n_pole_origin, finite_poles = split_origin(poles)
    if n_pole_origin:
        raise ValueError("s = 0 处的极点（积分器）不在本映射的支持范围内")
    num = np.prod(-finite_zeros) if finite_zeros.size else 1.0
    den = np.prod(-finite_poles) if finite_poles.size else 1.0
    ratio = np.real_if_close(num / den)
    if np.iscomplexobj(ratio) and abs(np.imag(ratio)) > 1e-9:
        raise ValueError("zpk 增益换算得到复数，检查共轭对是否完整")
    return float(gain * np.real(ratio))


def design(zeros, poles, fs: float, f_ref: float, sigma: float,
           prewarp: bool = False, gain: float = 1.0,
           omega: float | None = None) -> Design:
    """把 ``H(S) = gain * prod(1 - S/z)/prod(1 - S/p)`` 映射为 z 域因子。

    ``zeros`` / ``poles`` 是归一化变量 ``S = s/(2 pi f_ref)`` 下的根。
    原点零点单独计数并按 ``S -> (1-z^-1) omega/(2 pi)`` 映射；分子的阶数亏欠
    按退化因子补齐。``prewarp=True`` 时用 :func:`prewarped_omega` 代替
    ``fs/f_ref``，使数字参考频率精确落在 ``f_ref``。

    ``omega`` 给定时直接采用该设计比值（此时 ``f_ref`` 只用于标注原型，
    ``prewarp`` 不可再用）。需要自行解出 ``omega`` 的双点匹配见
    :mod:`two_point`。
    """
    if omega is None:
        omega = prewarped_omega(f_ref, fs, sigma) if prewarp else fs / f_ref
    elif prewarp:
        raise ValueError("omega 与 prewarp 不能同时给出：预畸变本身就是在求 omega")
    n_zero_origin, finite_zeros = split_origin(zeros)
    n_pole_origin, finite_poles = split_origin(poles)
    if n_pole_origin:
        raise ValueError("s = 0 处的极点（积分器）不在本映射的支持范围内")

    num: list[Factor] = []
    for kind, a, b in factorize(finite_zeros):
        num.append(map_linear(a, omega, sigma) if kind == "linear"
                   else map_quad(a, b, omega, sigma))
    num.extend(map_origin(omega) for _ in range(n_zero_origin))

    deficit = len(finite_poles) - (len(finite_zeros) + n_zero_origin)
    if deficit < 0:
        raise ValueError("分子次数高于分母，不是因果有理函数")
    num.extend(map_linear(0.0, omega, sigma) for _ in range(deficit))

    den: list[Factor] = []
    for kind, a, b in factorize(finite_poles):
        den.append(map_linear(a, omega, sigma) if kind == "linear"
                   else map_quad(a, b, omega, sigma))
    return Design(tuple(num), tuple(den), float(gain))


TINY = 1e-300
"""对数求值时避免 log(0) 的下限；不影响可听范围。"""


def _factor_logmag(fac: Factor, z_inv: np.ndarray) -> np.ndarray:
    value = fac.poly[0] + fac.poly[1] * z_inv
    if fac.poly.size == 3:
        value = value + fac.poly[2] * z_inv * z_inv
    return np.log(np.maximum(np.abs(fac.gain * value / fac.norm), TINY))


def digital_log_magnitude(design: Design, f, fs: float) -> np.ndarray:
    """数字滤波器的自然对数幅度；按因子累加，避免高频高阶下的溢出。"""
    w = 2.0 * np.pi * np.asarray(f, dtype=float) / fs
    z_inv = np.exp(-1j * w)
    logmag = np.full_like(w, np.log(max(abs(design.gain), TINY)))
    for fac in design.num:
        logmag += _factor_logmag(fac, z_inv)
    for fac in design.den:
        logmag -= _factor_logmag(fac, z_inv)
    return logmag


def digital_magnitude(design: Design, f, fs: float) -> np.ndarray:
    return np.exp(digital_log_magnitude(design, f, fs))


def analog_log_magnitude(zeros, poles, f, f_ref: float,
                         gain: float = 1.0) -> np.ndarray:
    """连续原型在频率 f（Hz）处的自然对数幅度，``S = j f / f_ref``。

    原点零点贡献 ``|S|^m``；``s = 0`` 处的极点不在支持范围内。
    """
    n_zero_origin, finite_zeros = split_origin(zeros)
    n_pole_origin, finite_poles = split_origin(poles)
    if n_pole_origin:
        raise ValueError("s = 0 处的极点（积分器）不在本映射的支持范围内")

    s = 1j * np.asarray(f, dtype=float) / f_ref
    logmag = np.full_like(np.abs(s), np.log(max(abs(gain), TINY)))
    if n_zero_origin:
        logmag += n_zero_origin * np.log(np.maximum(np.abs(s), TINY))
    for r in finite_zeros:
        logmag += np.log(np.maximum(np.abs(1.0 - s / r), TINY))
    for r in finite_poles:
        logmag -= np.log(np.maximum(np.abs(1.0 - s / r), TINY))
    return logmag


def analog_magnitude(zeros, poles, f, f_ref: float,
                     gain: float = 1.0) -> np.ndarray:
    return np.exp(analog_log_magnitude(zeros, poles, f, f_ref, gain))


def digital_zeros_poles(design: Design) -> tuple[np.ndarray, np.ndarray]:
    """把因子列表还原为 z 平面零点 / 极点（仅用于诊断）。"""

    def collect(factors: tuple[Factor, ...]) -> np.ndarray:
        out: list[complex] = []
        for fac in factors:
            if fac.poly.size == 2:
                out.append(-fac.poly[1] / fac.poly[0])
            else:
                out.extend(np.roots(fac.poly))
        return np.asarray(out, dtype=complex)

    return collect(design.num), collect(design.den)
