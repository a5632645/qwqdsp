"""High-shelf 模拟原型：复刻 ``iir_design_shelf.hpp`` 的构造。

原型以 ``(1) rad/s`` 为中点，输出 ``H(S)`` 的零点 / 极点（常数项归一化后的
zpk 根），且 ``H(0) = 1``。与 §4 的一般形式配合时，直接传这些根、
``gain = 1`` 即可；高频端电平由根的位置决定。

- Butterworth：特征函数 ``F = w^N``，零点 / 极点互为倒易，``H(inf) = 10^(boost/20)``。
- Chebyshev II：``H`` 的极点取“特征 ``F/b``”的自然切比雪夫 II 原型的极点，
  零点取“特征 ``bF``”的原型的极点；两组原型的**传输零点相同**，在对消后消失。
  频率尺度 ``m`` 把中点放到 1 rad/s。

本文件只实现 ``even_modify = False`` 的 Chebyshev II。scipy 的 ``cheb2ap`` 没有
偶数阶修正原型；修正版需要 C++ 头文件里的修改特征函数。未修正时高频端落在
棚带低谷：``H(inf) = boost - ripple`` dB（衰减时为 ``boost + ripple``），而棚带
内的极值仍在有限频率处达到 ``boost`` dB。这是模拟原型本身的性质，不是离散化
误差。

参考：Vadim Zavalishin, *The Art of VA Filter Design*, 2.1.2, Chapter 10；
以及本仓库的 ``qwqdsp/include/qwqdsp/filter/iir_design_shelf.hpp``。
"""
from __future__ import annotations

import numpy as np
from scipy import signal


def butterworth_shelf(num_pairs: int, boost_db: float) -> tuple[np.ndarray, np.ndarray]:
    """Butterworth high-shelf 原型（阶数 ``2 * num_pairs``）。

    零点在 ``p / s``、极点在 ``p * s``，``s = (10^(boost/20))^(1/N)``，
    ``p`` 是单位圆 Butterworth 极点。返回共轭完整的根。
    """
    if num_pairs < 1:
        raise ValueError("num_pairs 必须 >= 1")
    if abs(boost_db) < 1e-9:
        raise ValueError("boost_db = 0 时 shelf 退化为直通，不生成原型")
    n = 2 * num_pairs
    amp = 10.0 ** (boost_db / 20.0)
    b = np.sqrt(amp)
    s = b ** (1.0 / n)
    zeros: list[complex] = []
    poles: list[complex] = []
    for k in range(1, num_pairs + 1):
        phi = (2.0 * k - 1.0) * np.pi / (2.0 * n)
        p = complex(-np.sin(phi), np.cos(phi))
        zeros.extend([p / s, np.conj(p / s)])
        poles.extend([p * s, np.conj(p * s)])
    return np.asarray(zeros), np.asarray(poles)


def _rs_of_eps(eps: float) -> float:
    """特征函数参数 -> Chebyshev II 的阻带衰减（负 dB）。"""
    return 10.0 * np.log10(eps * eps / (1.0 + eps * eps))


def chebyshev2_shelf(num_pairs: int, boost_db: float,
                     ripple_db: float) -> tuple[np.ndarray, np.ndarray]:
    """Chebyshev II high-shelf 原型（``even_modify = False``）。

    参数与头文件一致：``num_pairs`` 是极点对数（阶数 ``2*num_pairs``），
    ``ripple_db`` 是棚带等波纹深度（dB，> 0）。``boost_db`` 可正可负。

    未修正 Chebyshev II 的高频端是 ``boost - sign(boost) * ripple_db`` dB；
    棚带内的极值仍到 ``boost`` dB。返回共轭完整的零点 / 极点，``H(0) = 1``。
    """
    if num_pairs < 1:
        raise ValueError("num_pairs 必须 >= 1")
    if abs(boost_db) < 1e-9:
        raise ValueError("boost_db = 0 时 shelf 退化为直通，不生成原型")
    if ripple_db <= 0.0:
        raise ValueError("ripple_db 必须 > 0")

    n = 2 * num_pairs
    beta2 = 10.0 ** (abs(boost_db) / 20.0)       # 幅度棚高
    beta4 = beta2 * beta2
    ripple = float(np.clip(ripple_db, 1e-3, abs(boost_db) * 0.999))
    rho = 10.0 ** (-ripple / 10.0)
    eps2_sq = (1.0 - rho) / (rho * beta4 - 1.0)
    if not (eps2_sq > 0.0):
        raise ValueError("eps2^2 <= 0：boost / ripple 组合无效")
    eps2 = np.sqrt(eps2_sq)
    eps1 = beta2 * eps2

    # 频率尺度：让中点落在 (1) rad/s
    r = 1.0 / (eps2 * np.sqrt(beta2))
    u = (np.cosh(np.arccosh(r) / n) if r >= 1.0
         else np.cos(np.arccos(min(r, 1.0)) / n))
    m = u

    rs1 = -_rs_of_eps(eps1)
    rs2 = -_rs_of_eps(eps2)
    _, p1, _ = signal.cheb2ap(n, rs1)
    _, p2, _ = signal.cheb2ap(n, rs2)
    p1 = np.asarray(p1, dtype=complex) * m
    p2 = np.asarray(p2, dtype=complex) * m

    if boost_db >= 0.0:
        return p2, p1
    return p1, p2
