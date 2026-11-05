"""一、二阶高通 / 低通 / 带通 / 带阻与一阶 shelf：同时精确匹配 fc 与 Nyquist。

论文固定 ``omega = fs/fc`` 与自由参数 ``sigma``，在 Nyquist 处匹配模拟原型
``fs/sigma`` 的值（(2.5)），再用 (2.10) 或 (2.11) 修截止频率：两件事只能做
一件。把设计比值 ``omega`` 与 ``sigma`` 都当作未知量，两个条件联立即得唯一解
——``sigma`` 不再自由，且恒不大于 2：

    F(omega) = fc                     (2.9) -> 截止频率
    |Hd(-1)| = |Ha(j fs/2)|           (2.5) -> Nyquist

两种类型的 ``(omega, sigma)`` 相同。截止条件对高通和低通都是

    cot^2(theta_c / 2) = (omega^2 - sigma^2) / pi^2,

即 (2.9) 的 ``cos(2 pi/u) = (omega^2 - pi^2 - sigma^2) / (omega^2 + pi^2 - sigma^2)``。
Nyquist 处高通的幅度是 ``omega / sqrt(omega^2 + sigma^2)``、低通是
``sigma / sqrt(omega^2 + sigma^2)``，两者都要等于模拟原型的
``u / sqrt(u^2 + 4)``，于是

    sigma = 2 omega / u,     u = fs / fc,

代入截止条件即得

    omega = pi cot(pi/u) / sqrt(1 - 4/u^2)
    sigma = 2 omega / u

一阶 shelf 共用同一组解。半增益点（幅度 ``sqrt(g)``）的条件化简后仍是
``cot^2(theta_c/2) = (omega^2 - sigma^2)/pi^2``；Nyquist 处的高 shelf 幅度是
``sqrt(omega^2 g + sigma^2)/sqrt(omega^2/g + sigma^2)``，与模拟原型在 ``fs/2``
处相等同样给出 ``sigma = 2 omega / u``——``g`` 在两个条件里都对消。DC 增益由
因子归一化自动精确（高 shelf 为 1、低 shelf 为 ``g``），所以 shelf 实际精确
匹配三点（DC、中点、Nyquist）；模拟原型的 ``g`` 是 ``s -> infinity`` 的渐近值，
数字端只逼近它。

同一组解对任意阶都成立。在 ``theta_c = 2 pi F(omega)/fs`` 处，映射的每个一次
因子（含 ``a = 0`` 的亏欠因子与原点零点）的数字幅度等于对应模拟因子在 ``S = j``
处的幅度乘以 ``r = omega / sqrt(pi^2 + omega^2 - sigma^2)``，每个二次因子则乘以
``r^2``——与系数无关。分子与分母的总次数相同（亏欠因子补齐），``r`` 的幂次相消，
于是 ``|Hd(fc)| = |Ha(fc)|`` 与阶数、类型无关。二阶只是把同一组
``(omega, sigma)`` 代入 (3.8)/Table 3.1，不需要新的解。

``u -> 2``（``fc -> fs/2``）时 ``omega, sigma -> 0``，极点趋于 ``z = -1``，
解退化；``u -> infinity`` 时 ``omega -> u``、``sigma -> 2``，即论文
``sigma = 2`` 的无预畸变设计。``fc >= fs/2`` 无解。

代价随类型不同：

- 高通：低频斜率被缩到 ``omega/u`` 倍，即模拟原型以下额外衰减
  ``20 log10(omega/u)`` dB（``u = 2.4`` 时 3.95 dB，``u = 24`` 时 0.02 dB）。
- 低通：DC 增益由因子归一化自动为 1，双点解因此同时精确匹配 DC、截止频率、
  Nyquist 三点；残差集中在高频端（``u = 2.4`` 时最大 0.56 dB 在过渡带，
  ``u = 2400`` 时 0.97 dB 出现在模拟已衰减到 -59 dB 处）。
- shelf：DC 与中点、Nyquist 都精确，残差集中在拐点附近，约按 boost 线性增长
  （``u = 2.4``、±12 dB 时 1.14 dB，±6 dB 时 0.58 dB；``u >= 24`` 时不超过
  0.06 dB）。
- 二阶高通：低频斜率被缩到 ``(omega/u)^2``，即 ``40 log10(omega/u)`` dB——比
  一阶高一倍（``u = 2.4`` 时 7.90 dB，``u = 3.2`` 时 3.02 dB）。二阶低通没有这
  一项（DC 增益归一化），残差仍在高频端（``u = 2.4`` 时 0.61 dB，``u = 9.6``
  时 1.85 dB，ζ = 1/sqrt(2)）。
- 二阶带通 / 带阻：中心点（带通峰值 1、带阻零点 0）与 Nyquist 都精确。带通有
  一个原点零点，低频侧因此同高通一样多衰减 ``20 log10(omega/u)`` dB；带阻的
  残差集中在陷波两侧的裙边（``|Ha| >= 0.5`` 内 ``u = 2.4`` 时 2.9 dB，
  ``u >= 24`` 时不超过 0.05 dB）。

``sigma`` 被压到论文幅度最优区间 ``[2, pi sqrt(2/3)]`` 之下。
"""
from __future__ import annotations

import numpy as np

import ivantsovy as I


def omega_sigma(fc: float, fs: float) -> tuple[float, float]:
    """双点匹配的 ``(omega, sigma)``；要求 ``0 < fc < fs/2``。"""
    if not 0.0 < fc < fs / 2.0:
        raise ValueError("双点匹配要求 0 < fc < fs/2：映射的截止频率达不到 fs/2")
    u = fs / fc
    omega = np.pi / np.tan(np.pi / u) / np.sqrt(1.0 - 4.0 / u ** 2)
    return float(omega), float(2.0 * omega / u)


def highpass(fc: float, fs: float) -> I.Design:
    """一阶高通 ``S/(S+1)`` 的双点匹配离散化。"""
    omega, sigma = omega_sigma(fc, fs)
    return I.design(np.array([0.0]), np.array([-1.0]), fs, fc, sigma, omega=omega)


def lowpass(fc: float, fs: float) -> I.Design:
    """一阶低通 ``1/(S+1)`` 的双点匹配离散化（DC 增益自动为 1）。"""
    omega, sigma = omega_sigma(fc, fs)
    return I.design(np.array([]), np.array([-1.0]), fs, fc, sigma, omega=omega)


def shelf_prototype(kind: str, boost_db: float) -> tuple[np.ndarray, np.ndarray, float]:
    """一阶 shelf 原型的 ``(zeros, poles, gain)``，中点在 ``fc``（归一化 ``S``）。

    高 shelf ``(1+sqrt(g) S)/(1+S/sqrt(g))``：``H(0)=1``、``H(inf)=g``；
    低 shelf 是其倒数乘 ``g``：``H(0)=g``、``H(inf)=1``。两者在 ``S=j`` 处
    的幅度都是 ``sqrt(g)``。
    """
    if boost_db == 0.0:
        raise ValueError("boost_db = 0 时 shelf 退化为直通")
    g = 10.0 ** (boost_db / 20.0)
    root = np.sqrt(g)
    if kind == "high":
        return np.array([-1.0 / root]), np.array([-root]), 1.0
    if kind == "low":
        return np.array([-root]), np.array([-1.0 / root]), g
    raise ValueError(f"未知 shelf 类型 {kind!r}")


def high_shelf(fc: float, fs: float, boost_db: float) -> I.Design:
    """一阶高 shelf；中点增益 ``sqrt(g)`` 精确落在 ``fc``。"""
    zeros, poles, gain = shelf_prototype("high", boost_db)
    omega, sigma = omega_sigma(fc, fs)
    return I.design(zeros, poles, fs, fc, sigma, omega=omega, gain=gain)


def low_shelf(fc: float, fs: float, boost_db: float) -> I.Design:
    """一阶低 shelf；中点增益 ``sqrt(g)`` 精确落在 ``fc``。"""
    zeros, poles, gain = shelf_prototype("low", boost_db)
    omega, sigma = omega_sigma(fc, fs)
    return I.design(zeros, poles, fs, fc, sigma, omega=omega, gain=gain)


def shelf_prototype2(kind: str, boost_db: float,
                     zeta: float) -> tuple[np.ndarray, np.ndarray, float]:
    """二阶 shelf 原型：``S^2 + 2 zeta S + 1`` 的根乘以 ``g^{±1/4}``。

    高 shelf 取 ``(1 + b S + a S^2)`` 的分子零点为 ``g^{-1/4} p``、极点为
    ``g^{1/4} p``，于是 ``H(0)=1``、``H(inf)=g``、``S=j`` 处幅度 ``sqrt(g)``；
    低 shelf 交换两者并乘 ``g``。``g`` 的指数是 ``1/4`` 而不是一阶的 ``1/2``，
    因为分子分母各升了一阶、每对根只承担一半的幅频跨度。
    """
    if boost_db == 0.0:
        raise ValueError("boost_db = 0 时 shelf 退化为直通")
    g = 10.0 ** (boost_db / 20.0)
    root = g ** 0.25
    core = second_order_poles(zeta)
    if kind == "high":
        return core / root, core * root, 1.0
    if kind == "low":
        return core * root, core / root, g
    raise ValueError(f"未知 shelf 类型 {kind!r}")


def high_shelf2(fc: float, fs: float, boost_db: float, zeta: float) -> I.Design:
    """二阶高 shelf；中点增益 ``sqrt(g)`` 精确落在 ``fc``。``zeta = 1/(2Q)``。"""
    zeros, poles, gain = shelf_prototype2("high", boost_db, zeta)
    omega, sigma = omega_sigma(fc, fs)
    return I.design(zeros, poles, fs, fc, sigma, omega=omega, gain=gain)


def low_shelf2(fc: float, fs: float, boost_db: float, zeta: float) -> I.Design:
    """二阶低 shelf；中点增益 ``sqrt(g)`` 精确落在 ``fc``。``zeta = 1/(2Q)``。"""
    zeros, poles, gain = shelf_prototype2("low", boost_db, zeta)
    omega, sigma = omega_sigma(fc, fs)
    return I.design(zeros, poles, fs, fc, sigma, omega=omega, gain=gain)


def second_order_poles(zeta: float) -> np.ndarray:
    """归一化二阶原型 ``S^2 + 2 zeta S + 1`` 的极点（共轭对或实根对）。"""
    if not zeta > 0.0:
        raise ValueError("zeta 必须为正，否则原型不稳定")
    if zeta < 1.0:
        r = np.sqrt(1.0 - zeta ** 2)
        return np.array([-zeta + 1j * r, -zeta - 1j * r])
    r = np.sqrt(zeta ** 2 - 1.0)
    return np.array([-zeta + r, -zeta - r])


def highpass2(fc: float, fs: float, zeta: float) -> I.Design:
    """二阶高通 ``S^2/(S^2 + 2 zeta S + 1)`` 的双点匹配离散化。"""
    omega, sigma = omega_sigma(fc, fs)
    return I.design(np.array([0.0, 0.0]), second_order_poles(zeta), fs, fc,
                    sigma, omega=omega)


def lowpass2(fc: float, fs: float, zeta: float) -> I.Design:
    """二阶低通 ``1/(S^2 + 2 zeta S + 1)`` 的双点匹配离散化。"""
    omega, sigma = omega_sigma(fc, fs)
    return I.design(np.array([]), second_order_poles(zeta), fs, fc,
                    sigma, omega=omega)


def bandpass(fc: float, fs: float, zeta: float) -> I.Design:
    """二阶带通 ``2 zeta S/(S^2 + 2 zeta S + 1)``；峰值精确落在 ``fc``。"""
    omega, sigma = omega_sigma(fc, fs)
    return I.design(np.array([0.0]), second_order_poles(zeta), fs, fc,
                    sigma, omega=omega, gain=2.0 * zeta)


def bandstop(fc: float, fs: float, zeta: float) -> I.Design:
    """二阶带阻 ``(S^2 + 1)/(S^2 + 2 zeta S + 1)``；零点精确落在 ``fc``。"""
    omega, sigma = omega_sigma(fc, fs)
    return I.design(np.array([1j, -1j]), second_order_poles(zeta), fs, fc,
                    sigma, omega=omega)
