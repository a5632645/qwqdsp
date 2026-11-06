"""三点功率匹配的公共工具。

- 功率函数的 phi 二次型 ``(1-s, s, 4s(1-s))``；
- 由分子功率系数 ``(B0,B1,B2)`` 反解 ``b0/b1/b2``；
- 稳定分母的 ``(p,q)`` 参数化与双线性极点；
- 以 Nelder–Mead 做等波纹（最小最大 ``|ΔdB|``）搜索。

低通/高通各自只需给出模拟功率与三点目标，再用 :func:`best_denominator`
选分母的两个自由值。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

PI = np.pi
"""圆周率。"""

FS = 48000.0
"""实验采样率。"""

FC_LIST = (20.0, 117.6, 691.4, 4065.0, 23900.0)
"""与其余 match_filter 实验一致的 5 个截止频率。"""

R2_SEEDS = (-0.9, -0.8, -0.75, -0.696, -0.6, -0.4, -0.2, 0.0, 0.3, 0.6)
"""第二极点的初值：低 fc 下低通/高通的最优第二极点约 -0.696 / -0.556。"""


def phi(w):
    """``(1 - s, s, 4s(1 - s))``，``s = sin^2(w/2)``。"""
    s = np.sin(np.asarray(w, dtype=float) / 2.0) ** 2
    return np.stack([1.0 - s, s, 4.0 * s * (1.0 - s)], axis=-1)


def biquad_power(b, a, w):
    """``|H(e^jw)|^2``；``a = (a1, a2)``（``a0 = 1``）。"""
    z = np.exp(-1j * np.asarray(w, dtype=float))
    num = b[0] + b[1] * z + b[2] * z * z
    den = 1.0 + a[0] * z + a[1] * z * z
    return np.abs(num / den) ** 2


@dataclass(frozen=True)
class Match:
    """三点匹配得到的二阶节；``a = (a1, a2)``，不含 ``a0 = 1``。"""

    b: np.ndarray
    a: np.ndarray
    p: float
    q: float
    beta: np.ndarray
    alpha: np.ndarray

    def power(self, w):
        return biquad_power(self.b, self.a, w)

    def poles(self):
        return np.roots([1.0, self.a[0], self.a[1]])

    def zeros(self):
        return np.roots([self.b[0], self.b[1], self.b[2]])


def alpha_of(p, q):
    """分母的 phi 二次型系数 ``(A0, A1, A2)``。"""
    return np.array([p * p, q * q, 4.0 - 2.0 * (p + q)])


def solve_b(beta):
    """由 ``(B0, B1, B2)`` 反解 ``b0/b1/b2``（Vicanek ``Solveb``）。"""
    W = (np.sqrt(beta[0]) + np.sqrt(beta[1])) / 2.0
    b0 = (W + np.sqrt(W * W + beta[2])) / 2.0
    return np.array([b0, (np.sqrt(beta[0]) - np.sqrt(beta[1])) / 2.0,
                     -beta[2] / (4.0 * b0)])


def feasible(beta):
    """分子能否还原成实数 ``b``；允许 ``B0 = 0``（高通在直流处为零）。"""
    W = (np.sqrt(beta[0]) + np.sqrt(beta[1])) / 2.0
    return bool(beta[0] >= 0.0 and beta[1] >= 0.0 and W * W + beta[2] > 0.0)


def bilinear_pole(fc, fs):
    """一阶双线性变换的极点位置（分母 ``1+a1 z^-1`` 时 ``z = -a1``）。"""
    wc = 2.0 * PI * fc / fs
    k = 1.0 / np.tan(wc / 2.0)
    return (k - 1.0) / (k + 1.0)


def bilinear_poles(fc, fs, Q):
    """二阶模拟极点对的双线性像（复共轭对或一对实极点）。"""
    zeta = 1.0 / (2.0 * Q)
    wc = 2.0 * PI * fc / fs
    t = np.tan(wc / 2.0)
    if zeta < 1.0:
        mu = np.sqrt(1.0 - zeta * zeta)
        z = (1.0 + (-zeta + 1j * mu) * t) / (1.0 - (-zeta + 1j * mu) * t)
        return np.array([z, np.conj(z)])
    if zeta == 1.0:
        z = (1.0 - zeta * t) / (1.0 + zeta * t)
        return np.array([z, z])
    r = np.sqrt(zeta * zeta - 1.0)
    z1 = (1.0 + (-zeta + r) * t) / (1.0 - (-zeta + r) * t)
    z2 = (1.0 + (-zeta - r) * t) / (1.0 - (-zeta - r) * t)
    return np.array([z1, z2])


def refine(objective, p0, q0, maxfev=2500):
    """Nelder–Mead 局部最小化 ``objective(p, q)``；不可行返回 ``inf``。"""
    from scipy.optimize import minimize

    def wrapped(x):
        value = objective(x[0], x[1])
        return 1e3 if not np.isfinite(value) else value

    res = minimize(wrapped, np.array([p0, q0], dtype=float), method="Nelder-Mead",
                   options={"xatol": 1e-11, "fatol": 1e-11,
                            "maxiter": 2000, "maxfev": maxfev})
    return res.x, float(res.fun)


def best_denominator(fc, fs, max_error, r1=None):
    """等波纹意义下最优的分母两个自由值 ``(p, q)``。

    ``max_error(p, q, n)`` 返回全带 ``max|ΔdB|``（不可行时 ``inf``）。以 ``r1``
    （默认双线性极点）为第一个极点、扫描第二极点生成初值，粗网格选优后在细
    网格上收敛。shelf 的 ``r1`` 是模拟极点 ``wc*sqrt(g)`` 的双线性像。
    """
    if r1 is None:
        r1 = bilinear_pole(fc, fs)
    best = None
    for r2 in R2_SEEDS:
        p0 = (1.0 - r1) * (1.0 - r2)
        q0 = (1.0 + r1) * (1.0 + r2)
        x, value = refine(lambda p, q: max_error(p, q, 2001), p0, q0)
        if best is None or value < best[1]:
            best = (x, value)
    x, _ = refine(lambda p, q: max_error(p, q, 4001), best[0][0], best[0][1])
    return float(x[0]), float(x[1])
