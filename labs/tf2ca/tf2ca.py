# -*- coding: utf-8 -*-
"""
tf2ca.py
========

把**奇数阶**数字 IIR 低通滤波器 ``H(z) = B(z)/A(z)`` 分解为两条全通链之和：

    H(z) = 1/2 * ( A0(z) + A1(z) )

同时给出功率互补的高通：

    H_c(z) = 1/2 * ( A0(z) - A1(z) )

这正是 MATLAB ``tf2ca`` / ``ca2tf`` 所做的事，也是 EMQF（Elliptic Minimal
Q-Factor）滤波器「一个滤波器加一个加法」实现的理论基础。

数学背景
--------
设 A0 = N0/D0、A1 = N1/D1 为全通（N_i 是 D_i 的系数反序，即 N_i = D_i^#），
代入并通分：

    H = (N0*D1 + N1*D0) / (2*D0*D1)   =>   A = D0*D1,  B = (N0*D1 + N1*D0)/2

因此可分解的必要条件是：

1. 阶数 n = deg A 为**奇数**（D0/D1 阶数一奇一偶，保证 B(-1)=0，即低通）；
2. 分子 B 为**镜像对称多项式**（palindromic，``b[k] == b[n-k]``）；
3. 归一化到 H(1)=1。
   ``butter`` / ``cheby1`` / ``cheby2`` / ``ellip`` 的输出天然满足 2、3。

令 B_c 为互补高通的分子，则

    B + B_c = N0*D1 ,   B - B_c = N1*D0

两者相乘即得

    B_c^2 = B^2 - A*A^#          （A^# = A 的镜像多项式）

B_c 满足**反镜像对称** ``c[k] == -c[n-k]``（奇数阶，故 B_c(1)=0）。

极点分组
--------
对 A 的极点 p（互异、单根）：

    p 属于 D0  <=>  B_c(p) = +B(p)
    p 属于 D1  <=>  B_c(p) = -B(p)

等价地 ``D0 = A / gcd(A, B+B_c)``、``D1 = A / gcd(A, B-B_c)``。

工程上对 butter/cheby1/cheby2/ellip 四族，规则可以简化为：

    把 A 的所有极点按模 |p| 升序排列（1 个实极点 + (n-1)/2 个共轭对），
    逐个（共轭对整体）交替分配到两条链。

实极点所在的链阶数为奇数。**注意不是按幅角交替**（见 README 中 cheby2 的反例）。

本模块提供两种实现：

- :func:`split_by_radius` —— 上面的工程规则，数值稳健，适合四族经典设计，
  与 MATLAB ``tf2ca`` 发布样例、EMQF ``apellip_du`` 的 beta 排序交替完全一致。
- :func:`split_by_bc`    —— 精确的代数判据，对任意可分解滤波器成立；但窄带
  高阶时 ``B^2 - A*A^#`` 会发生消位误差，见 README。
"""
from __future__ import annotations

import numpy as np

# ------------------------------------------------------------
# 极点分组为「一阶 / 二阶」节点
# ------------------------------------------------------------


def sections(a: np.ndarray) -> list[dict]:
    """把分母 A 的极点整理成实极点(一阶)与共轭对(二阶)两种节点。

    返回列表，每项为 ``{"poles": [...], "poly": [...]}``：

    - ``poles``：节点的极点（实极点 1 个，共轭对 2 个）；
    - ``poly`` ：该节点的分母多项式系数（降幂，首项为 1）。

    使用多项式系数而不是极点本身做卷积，可以避免反复 roots/poly 往返。
    """
    poles = np.roots(np.asarray(a, dtype=float))
    real = sorted((p.real for p in poles if abs(p.imag) < 1e-7))
    pairs = sorted(
        (p for p in poles if p.imag > 1e-7),
        key=lambda p: abs(p),
    )
    out: list[dict] = []
    for p in real:
        out.append({"poles": [p], "poly": np.array([1.0, -p])})
    for p in pairs:
        out.append(
            {
                "poles": [p, p.conjugate()],
                "poly": np.poly([p, p.conjugate()]).real,
            }
        )
    return out


def allpass_from_denoms(denoms: list[np.ndarray]) -> np.ndarray:
    """把若干全通节点的分母多项式串接（卷积）成整条链的分母 D。"""
    d = np.array([1.0])
    for poly in denoms:
        d = np.convolve(d, poly)
    return d


# ------------------------------------------------------------
# 方法一：按极点模长交替（工程规则）
# ------------------------------------------------------------


def split_by_radius(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """按极点模 |p| 升序交替分配节点，返回 (D0, D1)。

    与 MATLAB ``tf2ca``、EMQF ``apellip_du`` 的结果一致（分支可互换）。
    """
    nodes = sections(a)
    nodes.sort(key=lambda s: abs(s["poles"][0]))
    d = [np.array([1.0]), np.array([1.0])]
    for k, node in enumerate(nodes):
        d[k % 2] = np.convolve(d[k % 2], node["poly"])
    return d[0], d[1]


# ------------------------------------------------------------
# 方法二：精确代数判据
# ------------------------------------------------------------


def bc_antisymmetric_sqrt(b: np.ndarray, a: np.ndarray) -> np.ndarray:
    """求反镜像分子 B_c，满足 ``B_c^2 = B^2 - A*A^#``。

    设 ``c`` 为 B_c 的降幂系数，则 ``c[k] = -c[n-k]``。对 ``k <= (n-1)/2``
    的系数方程 ``(c*c)[k] == R[k]`` 是**下三角**的，可逐项递推：

        c[0] = sqrt(R[0])
        c[k] = ( R[k] - sum_{i=1}^{k-1} c[i]*c[k-i] ) / (2*c[0])

    其余系数由反镜像对称补出，且 ``k > (n-1)/2`` 的方程自动成立。
    """
    b = np.asarray(b, dtype=float)
    a = np.asarray(a, dtype=float)
    n = a.size - 1
    if n % 2 != 1:
        raise ValueError("阶数必须为奇数")
    r = np.convolve(b, b) - np.convolve(a, a[::-1])
    m = (n - 1) // 2
    c = np.zeros(n + 1)
    c[0] = np.sqrt(r[0])
    for k in range(1, m + 1):
        acc = sum(c[i] * c[k - i] for i in range(1, k))
        c[k] = (r[k] - acc) / (2.0 * c[0])
    for k in range(n - m, n + 1):
        c[k] = -c[n - k]
    return c


def split_by_bc(b: np.ndarray, a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """精确判据：用 B_c(p) 与 B(p) 的符号决定每个极点节点归哪条链。

    对任意可分解的低通都成立；但窄带/高阶时 ``B^2 - A*A^#`` 存在消位误差，
    精度不如 :func:`split_by_radius`。
    """
    b = np.asarray(b, dtype=float)
    a = np.asarray(a, dtype=float)
    c = bc_antisymmetric_sqrt(b, a)
    d = [np.array([1.0]), np.array([1.0])]
    for node in sections(a):
        p = node["poles"][0]
        bc_p = np.polyval(c, p)
        b_p = np.polyval(b, p)
        idx = 0 if abs(bc_p - b_p) < abs(bc_p + b_p) else 1
        d[idx] = np.convolve(d[idx], node["poly"])
    return d[0], d[1]


# ------------------------------------------------------------
# 主入口 / 校验
# ------------------------------------------------------------


def tf2ca(b: np.ndarray, a: np.ndarray, method: str = "radius") -> tuple[np.ndarray, np.ndarray]:
    """把低通 ``b/a`` 分解为两条全通链的分母 (D0, D1)。

    :param method: ``"radius"``（默认，四族经典设计）或 ``"bc"``（通用精确）。
    """
    if method == "radius":
        return split_by_radius(a)
    if method == "bc":
        return split_by_bc(b, a)
    raise ValueError(f"未知方法: {method}")


def reconstruct(d0: np.ndarray, d1: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """由 (D0, D1) 还原 (B, A, B_c)。

    其中 ``B = (N0*D1 + N1*D0)/2``、``B_c = (N0*D1 - N1*D0)/2``、``A = D0*D1``。
    """
    n0, n1 = np.asarray(d0)[::-1], np.asarray(d1)[::-1]
    b = (np.convolve(n0, d1) + np.convolve(n1, d0)) / 2.0
    bc = (np.convolve(n0, d1) - np.convolve(n1, d0)) / 2.0
    return b, np.convolve(d0, d1), bc


def zfreq(poly: np.ndarray, w: np.ndarray) -> np.ndarray:
    """求 ``Σ_k poly[k] * z^-k`` 在 ``z = exp(jw)`` 上的值。

    降幂系数 ``poly`` 表示 ``Σ_k poly[k] z^-k``（与 scipy ``freqz`` 一致），
    所以要先反序再代入 ``z^-1``。
    """
    zinv = np.exp(-1j * np.asarray(w, dtype=float))
    return np.polyval(np.asarray(poly)[::-1], zinv)


def tf_response(b: np.ndarray, a: np.ndarray, w: np.ndarray) -> np.ndarray:
    """``H(z) = B(z)/A(z)`` 在 ``z = exp(jw)`` 上的复响应。"""
    return zfreq(b, w) / zfreq(a, w)


def allpass_response(d: np.ndarray, w: np.ndarray) -> np.ndarray:
    """全通链 ``A(z) = N(z)/D(z)`` 在 ``z = exp(jw)`` 上的复响应（N 为 D 的系数反序）。"""
    return zfreq(np.asarray(d)[::-1], w) / zfreq(d, w)


def response(d0: np.ndarray, d1: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """返回并行全通结构的低通 / 高通复响应 (H, H_c)。

    H = (A0+A1)/2 为低通，H_c = (A0-A1)/2 为功率互补高通，|H|^2 + |H_c|^2 = 1。
    """
    a0 = allpass_response(d0, w)
    a1 = allpass_response(d1, w)
    return (a0 + a1) / 2.0, (a0 - a1) / 2.0


def to_biquads(d: np.ndarray) -> list[tuple[int, np.ndarray]]:
    """把整条链的分母 D 拆成 1 阶 / 2 阶全通节点，便于逐节点实现。

    返回 ``[(order, [d0, d1, ...]), ...]``，每个节点的分子为其分母的反序。
    """
    out: list[tuple[int, np.ndarray]] = []
    for node in sections(d):
        poly = node["poly"]
        out.append((poly.size - 1, poly))
    return out


# ------------------------------------------------------------
# 设计辅助与误差度量
# ------------------------------------------------------------


def design(kind: str, order: int, cutoff: float, rp: float = 0.2, rs: float = 40.0):
    """用 scipy 设计低通，``cutoff`` 以 Nyquist(=1) 归一。返回 (b, a)。"""
    from scipy import signal

    if kind == "butter":
        return signal.butter(order, cutoff, "low")
    if kind == "cheby1":
        return signal.cheby1(order, rp, cutoff, "low")
    if kind == "cheby2":
        return signal.cheby2(order, rs, cutoff, "low")
    if kind == "ellip":
        return signal.ellip(order, rp, rs, cutoff, "low")
    raise ValueError(f"未知类型: {kind}")


def max_response_error(b, a, d0, d1, num: int = 4000) -> float:
    """并行全通结构与原传递函数在频响上的最大偏差（复数，含相位）。"""
    w = np.linspace(1e-6, np.pi - 1e-6, num)
    ref = tf_response(b, a, w)
    got, _ = response(d0, d1, w)
    return float(np.max(np.abs(got - ref)))


def allpass_flatness(d0, d1, num: int = 4000) -> float:
    """两条链 |A_i(e^jw)| 偏离 1 的最大量（应约等于机器精度）。"""
    w = np.linspace(1e-6, np.pi - 1e-6, num)
    return float(
        max(
            np.max(np.abs(np.abs(allpass_response(d0, w)) - 1.0)),
            np.max(np.abs(np.abs(allpass_response(d1, w)) - 1.0)),
        )
    )
