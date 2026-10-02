# -*- coding: utf-8 -*-
"""
arg_s_plane.py
==============

验证"按角度交替"的一种变体：把数字极点**先映射回 s 平面**（双线性反变换），再测 s 平面角

    s = 2*fs*(p - 1)/(p + 1)        fs = 2  ==>  s = 4*(p - 1)/(p + 1)
    theta(p) = arg(s)

这正是双线性变换下真正的"设计角"：它与模拟原型极点的角度只差一个**正实**比例
（``lp2lp`` 的频率缩放不改变角度），所以排序与设计序号一致。

数值路径**全程不碰多项式**：

- 极点：``iirfilter(..., output="zpk")`` 直接给；
- 参考响应：``freqz_zpk``；
- 每条全通链：极点已知，零点就是 ``1/conj(p)``，同样交给 ``freqz_zpk``。

（之前用 ``np.roots(a)`` + ``B/A`` 多项式求值时，窄带高阶例子的条件数误差有 1e-2 量级，
会被误读成"规则失败"。）

扫描网格：四族 × 奇数阶 3/5/7/9/11 × cutoff 0.1/0.2/0.3/0.35/0.5，rp = 0.2 dB、rs = 40 dB，
共 100 例。**只打印出问题的例子**：重建残差 ``max|dH| > 1e-8``，或分区与 ``|p|`` 规则不同。
共轭对取上半平面那一个极点算角度。
"""
from __future__ import annotations

import numpy as np
from scipy import signal as sg

FAMILIES = ("butter", "cheby1", "cheby2", "ellip")
ORDERS = (3, 5, 7, 9, 11)
CUTOFFS = (0.1, 0.2, 0.3, 0.35, 0.5)
RP = 0.2
RS = 40.0
TOL = 1e-8
FS = 2.0


def to_s(p: complex) -> complex:
    """双线性反变换：数字极点 -> s 平面极点（``2*fs`` 只是正实比例，不影响 arg）。"""
    return 2.0 * FS * (p - 1.0) / (p + 1.0)


def design_zpk(kind: str, order: int, cutoff: float, rp: float, rs: float):
    """scipy 直接给出数字滤波器的 zpk（``Wn`` 仍以 Nyquist = 1 归一）。"""
    return sg.iirfilter(order, cutoff, rp, rs, btype="low", ftype=kind, output="zpk")


def nodes_from_zpk(p: np.ndarray) -> list[complex]:
    """每个共轭对取上半平面那一个代表，实极点单独一个，作为一节。"""
    return [q for q in np.asarray(p) if q.imag >= 0.0]


def split_order(nodes: list[complex], key) -> list[int]:
    """按 key(node) 升序给出名次排列（名次奇偶 = 归哪条链）。"""
    return sorted(range(len(nodes)), key=lambda i: key(nodes[i]))


def chain_zpk(nodes: list[complex], order: list[int], chain: int):
    """一条全通链的 (零点, 极点)：零点 = 极点的共轭倒数。"""
    ps = []
    for k, i in enumerate(order):
        if k % 2 != chain:
            continue
        q = nodes[i]
        ps.append(q)
        if abs(q.imag) > 1e-12:
            ps.append(np.conj(q))
    ps = np.asarray(ps)
    return 1.0 / np.conj(ps), ps


def chain_response(nodes, order, chain, w):
    """该链的全通响应，用 zpk 因子形式求值。

    零点取极点共轭倒数之后，因子式与 ``tf2ca.allpass_response``（零点 = 极点倒数、
    ``A(1) = 1``）只差一个增益，这里按 z = 1 处的值归一，免得高阶链的增益因子
    （``prod(1/|p|)`` 量级）把数值淹没。
    """
    z, p = chain_zpk(nodes, order, chain)
    h = sg.freqz_zpk(z, p, 1.0, worN=w)[1]
    h0 = sg.freqz_zpk(z, p, 1.0, worN=[0.0])[1][0]
    return h / h0


def partition_key(nodes: list[complex], order: list[int]) -> tuple:
    """分区指纹：两条链各自的极点集合（允许互换）。"""
    chains = []
    for chain in (0, 1):
        _, ps = chain_zpk(nodes, order, chain)
        chains.append(tuple(np.round(np.sort_complex(ps), 9)))
    return tuple(sorted(chains))


def main() -> None:
    w = np.linspace(0, np.pi, 2000)
    n_case = 0
    n_only_this = 0
    n_both = 0
    worst = 0.0
    worst_r = 0.0
    print("== arg(s) 规则（先双线性反变换回 s 平面；极点取自 scipy zpk）：只列出有问题的例子 ==")
    for kind in FAMILIES:
        for order in ORDERS:
            for cut in CUTOFFS:
                z, p, k = design_zpk(kind, order, cut, RP, RS)
                nodes = nodes_from_zpk(p)
                ref = sg.freqz_zpk(z, p, k, worN=w)[1]

                ord1 = split_order(nodes, lambda q: np.angle(to_s(q)))
                got = (chain_response(nodes, ord1, 0, w)
                       + chain_response(nodes, ord1, 1, w)) / 2.0
                err = float(np.max(np.abs(got - ref)))

                ordr = split_order(nodes, lambda q: abs(q))
                got_r = (chain_response(nodes, ordr, 0, w)
                         + chain_response(nodes, ordr, 1, w)) / 2.0
                err_r = float(np.max(np.abs(got_r - ref)))
                same = partition_key(nodes, ord1) == partition_key(nodes, ordr)

                worst = max(worst, err)
                worst_r = max(worst_r, err_r)
                n_case += 1
                if err > TOL or not same:
                    if err_r > TOL:
                        tag = "两者都坏"
                        n_both += 1
                    else:
                        tag = "仅本规则"
                        n_only_this += 1
                    print(f"  {kind:7s} N={order:<3d} cutoff={cut:<6g} "
                          f"arg(s) max|dH|={err:.3e}   |p| max|dH|={err_r:.3e}   "
                          f"分区与|p|{'相同' if same else '不同'}   [{tag}]")
    print(f"\n共 {n_case} 例：本规则出问题 {n_only_this + n_both} 例"
          f"（其中仅本规则 {n_only_this} 例，与 |p| 同时出问题 {n_both} 例）")
    print(f"  最大残差：本规则 {worst:.3e}，|p| 规则 {worst_r:.3e}")


if __name__ == "__main__":
    main()
