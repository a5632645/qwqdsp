"""Vicanek `match_biquad.hpp` 的 Python 移植（对照用）。

与 Ivantsov 的根本区别：**它不是频率变换**。

- 极点用脉冲响应不变：``a1 = -2 exp(-zeta wc) cos(sqrt(1-zeta^2) wc)``、
  ``a2 = exp(-2 zeta wc)``——极点角度就是 ``wc``（不畸变），半径 ``exp(-zeta wc)``。
- 分子（零点）再按"幅度平方在若干频点上等于模拟原型"拟合；``|H|^2`` 是
  ``phi = (1-s, s, 4s(1-s))``（``s = sin^2(w/2)``）的二次型，所以"3 个频点"
  就是关于分子三个系数的 3x3 线性方程组。

于是没有全局映射：不同滤波器类型（低通 / 高通 / shelf）拟合的频点、目标值都不同，
等价映射曲线也就各不相同；低通/高通只匹配 DC 与 ``wc``（Nyquist 不管），shelf 在
可行时匹配 DC / Nyquist / ``wc`` 三点，不可行时（极点频率 ``wc*sqrt_G`` 越过 Nyquist）
退化为重分布匹配点。

本文件只移植幅度相关部分，数值与 C++ 头文件一致。
"""
from __future__ import annotations

import numpy as np

PI = np.pi


def impulse_invariant(wc: float, Q: float) -> tuple[float, float]:
    """脉冲响应不变的极点（对应 C++ ``ImpluseInvarant``）。"""
    zeta = 1.0 / (2.0 * Q)
    e = np.exp(-zeta * wc)
    if zeta <= 1.0:
        a1 = -2.0 * e * np.cos(np.sqrt(1.0 - zeta * zeta) * wc)
    else:
        a1 = -2.0 * e * np.cosh(np.sqrt(zeta * zeta - 1.0) * wc)
    return float(a1), float(e * e)


def get_phi(w):
    s = np.sin(np.asarray(w, dtype=float) / 2.0) ** 2
    return np.stack([1.0 - s, s, 4.0 * (1.0 - s) * s], axis=-1)


def get_A(a1: float, a2: float) -> np.ndarray:
    return np.array([(1.0 + a1 + a2) ** 2, (1.0 - a1 + a2) ** 2, -4.0 * a2])


def solve_b(a1: float, a2: float, B0: float, B1: float, B2: float):
    """由三个约束解出 b0/b1/b2（对应 C++ ``Solveb``，系数归一 ``a0 = 1``）。"""
    W = (np.sqrt(max(B0, 0.0)) + np.sqrt(max(B1, 0.0))) / 2.0
    b0 = (W + np.sqrt(max(W * W + B2, 0.0))) / 2.0
    b1 = (np.sqrt(max(B0, 0.0)) - np.sqrt(max(B1, 0.0))) / 2.0
    b2 = -B2 / (4.0 * b0) if b0 != 0.0 else 0.0
    return np.array([b0, b1, b2]), np.array([1.0, a1, a2])


def magnitude(b, a, w) -> np.ndarray:
    zi = np.exp(-1j * np.asarray(w, dtype=float))
    return np.abs((b[0] + b[1] * zi + b[2] * zi * zi) /
                  (a[0] + a[1] * zi + a[2] * zi * zi))


def lowpass(wc: float, Q: float):
    """2 极点低通：DC = 1，``|H(wc)| = Q``（Nyquist 不匹配）。"""
    a1, a2 = impulse_invariant(wc, Q)
    A = get_A(a1, a2)
    phi = get_phi(wc)
    R = Q * Q * float(np.dot(phi, A))
    B0 = A[0]
    B1 = max((R - B0 * phi[0]) / phi[1], 0.0)
    d1 = 1.0 + a1 + a2
    b0 = (d1 + np.sqrt(B1)) / 2.0
    return np.array([b0, d1 - b0, 0.0]), np.array([1.0, a1, a2])


def highpass(wc: float, Q: float):
    """2 极点高通：DC = 0，``|H(wc)| = Q``（Nyquist 不匹配）。"""
    a1, a2 = impulse_invariant(wc, Q)
    A = get_A(a1, a2)
    phi = get_phi(wc)
    R = max(float(np.dot(phi, A)), 0.0)
    b0 = np.sqrt(R) * Q / (4.0 * phi[1])
    return np.array([b0, -2.0 * b0, b0]), np.array([1.0, a1, a2])


def shelf_target(db: float, wc: float, Q: float, kind: str):
    """2 极点 shelf 的模拟原型 ``|H|^2``（对应 ``analog_responce.hpp``）。"""
    A = 10.0 ** (db / 40.0)
    sa = np.sqrt(A)

    def target(w):
        s = 1j * np.asarray(w, dtype=float) / wc
        s2 = s * s
        if kind == "high":
            h = A * (A * s2 + sa / Q * s + 1.0) / (s2 + sa / Q * s + A)
        else:
            h = A * (s2 + sa / Q * s + A) / (A * s2 + sa / Q * s + 1.0)
        h = np.abs(h) ** 2
        return float(h) if np.ndim(h) == 0 else h

    return target


def fit_numerator(a1: float, a2: float, wc: float, target):
    """三点拟合分子（对应 C++ ``FitNumerator``，含不可行时的退化分支）。"""
    A = get_A(a1, a2)
    phi = get_phi(wc)
    B0 = A[0] * target(0.0)
    B1 = A[1] * target(PI)
    R1 = target(wc) * float(np.dot(phi, A))
    B2 = (R1 - B0 * phi[0] - B1 * phi[1]) / phi[2]

    def feasible(B) -> bool:
        if not (B[0] > 0.0 and B[1] > 0.0):
            return False
        W = (np.sqrt(B[0]) + np.sqrt(B[1])) / 2.0
        return W * W + B[2] > 0.0

    if feasible([B0, B1, B2]):
        return solve_b(a1, a2, B0, B1, B2)

    w_mid, w_high = 0.5 * wc, wc
    for _ in range(40):
        ps = np.array([get_phi(0.0), get_phi(w_mid), get_phi(w_high)])
        rhs = np.array([target(0.0) * float(np.dot(ps[0], A)),
                        target(w_mid) * float(np.dot(ps[1], A)),
                        target(w_high) * float(np.dot(ps[2], A))])
        try:
            cand = np.linalg.solve(ps, rhs)
        except np.linalg.LinAlgError:
            cand = None
        if cand is not None and feasible(list(cand)):
            return solve_b(a1, a2, *cand)
        w_mid, w_high = 0.5 * (w_mid + w_high), 0.5 * (w_high + PI)
    return solve_b(a1, a2, B0, B1, B2)


def highshelf(wc: float, Q: float, db: float):
    sa = 10.0 ** (db / 80.0)
    a1, a2 = impulse_invariant(wc * sa, Q)
    return fit_numerator(a1, a2, wc, shelf_target(db, wc, Q, "high"))


def lowshelf(wc: float, Q: float, db: float):
    sa = 10.0 ** (db / 80.0)
    a1, a2 = impulse_invariant(wc / sa, Q)
    return fit_numerator(a1, a2, wc, shelf_target(db, wc, Q, "low"))


def shelf_roots(db: float, Q: float) -> tuple[np.ndarray, np.ndarray]:
    """shelf 原型的零点 / 极点（归一化 ``r = w/wc``），供等价映射使用。"""
    A = 10.0 ** (db / 40.0)
    sa = np.sqrt(A)
    zeta = 1.0 / (2.0 * Q)
    # 先解 r = w/wc 平面（分母 r^2 - j(sa/Q) r - A = 0），再乘 j 换到 S = j r
    poles = 1j * np.array([0.5 * (1j * sa / Q + np.sqrt(-(sa / Q) ** 2 + 4.0 * A)),
                           0.5 * (1j * sa / Q - np.sqrt(-(sa / Q) ** 2 + 4.0 * A))])
    zeros = 1j * np.array([0.5 * (1j * sa / (A * Q) + np.sqrt(-(sa / (A * Q)) ** 2 + 4.0 / A)),
                           0.5 * (1j * sa / (A * Q) - np.sqrt(-(sa / (A * Q)) ** 2 + 4.0 / A))])
    return zeros, poles
