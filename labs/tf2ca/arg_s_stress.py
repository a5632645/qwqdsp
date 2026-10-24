# -*- coding: utf-8 -*-
"""
arg_s_stress.py
===============

对"把数字极点映回 s 平面，按 ``arg(s)`` 升序交替"这条规则做**更严重 / 更细分 / 更极端**的
压力测试。映射与求值路径与 ``arg_s_plane.py`` 完全一致（全程 zpk，不碰多项式）：

    s = 2*fs*(p - 1)/(p + 1)        fs = 2
    s 平面角升序 -> 逐个节点交替分到两条链
    残差 = max|H_重建 - H_设计|，两侧都用 freqz_zpk

三组扫描：

A. **细分截止频率**：orders 3,5,…,31 × cutoff（对数刻度 5e-4 … 0.4999），rp=0.2/rs=40
B. **极端规格**：rp ∈ {0.001, 0.01, 0.1, 1, 3} × rs ∈ {20, 40, 80, 120, 160}（取 rs > rp）
   × orders {3,7,11,15,21} × cutoff {0.05, 0.2, 0.45}
C. **超高阶**：orders {31,41,51,61,71} × cutoff {0.1,0.3,0.45} × (rp,rs) ∈ {(0.2,40),(1,120)}

只打印有问题的例子（残差 > TOL，或分区与 ``|p|`` 规则不同），另给出每组统计，包括
**arg(s) 相邻间隔最小值**（deg）—— 用来判断失败到底是精度问题还是结构问题。

``--controls`` 会拿一条**已知会坏**的规则（数字域 ``arg(p)``）跑同一套件，作为阳性对照：
如果套件连它都报不出来，那"0 例出问题"就没有意义。
"""
from __future__ import annotations

import argparse
import collections
import warnings

import numpy as np
from scipy import signal as sg

FS = 2.0
TOL = 1e-8
NUMW = 1200

w = np.linspace(0, np.pi, NUMW)


def to_s(p: complex) -> complex:
    """双线性反变换：数字极点 -> s 平面极点。"""
    return 2.0 * FS * (p - 1.0) / (p + 1.0)


def build(kind: str, order: int, cut: float, rp: float, rs: float):
    """scipy 直接给 zpk；节点 = 上半平面极点（每对取一个代表）。"""
    z, p, k = sg.iirfilter(order, cut, rp, rs, btype="low", ftype=kind, output="zpk")
    return z, p, k, [q for q in np.asarray(p) if q.imag >= 0.0]


def split_order(nodes, key) -> list[int]:
    return sorted(range(len(nodes)), key=lambda i: key(nodes[i]))


def chain_zpk(nodes, order: list[int], chain: int):
    ps = []
    for n, i in enumerate(order):
        if n % 2 != chain:
            continue
        q = nodes[i]
        ps.append(q)
        if abs(q.imag) > 1e-12:
            ps.append(np.conj(q))
    ps = np.asarray(ps)
    return 1.0 / np.conj(ps), ps


def chain_resp(nodes, order, chain) -> np.ndarray:
    """全通链响应（零点 = 极点共轭倒数，按 z=1 归一到 A(1)=1）。"""
    z, p = chain_zpk(nodes, order, chain)
    h = sg.freqz_zpk(z, p, 1.0, worN=w)[1]
    h0 = sg.freqz_zpk(z, p, 1.0, worN=[0.0])[1][0]
    return h / h0


def rebuilt(nodes, order) -> np.ndarray:
    return (chain_resp(nodes, order, 0) + chain_resp(nodes, order, 1)) / 2.0


def partition_key(nodes, order) -> tuple:
    ch = []
    for c in (0, 1):
        _, ps = chain_zpk(nodes, order, c)
        ch.append(tuple(np.round(np.sort_complex(ps), 9)))
    return tuple(sorted(ch))


def min_gap_deg(nodes) -> float:
    """arg(s) 排序后相邻间隔的最小值（度）。越小越脆弱。"""
    a = np.sort([np.degrees(np.angle(to_s(q))) for q in nodes])
    return float(np.min(np.diff(a))) if len(a) > 1 else 180.0


def min_gap_abs(nodes) -> float:
    """|p| 排序后相邻间隔的最小值（绝对值）。"""
    a = np.sort([abs(q) for q in nodes])
    return float(np.min(np.diff(a))) if len(a) > 1 else 1.0


def run_case(kind, order, cut, rp, rs):
    z, p, k, nodes = build(kind, order, cut, rp, rs)
    ref = sg.freqz_zpk(z, p, k, worN=w)[1]
    o_s = split_order(nodes, lambda q: np.angle(to_s(q)))
    o_r = split_order(nodes, lambda q: abs(q))
    e_s = float(np.max(np.abs(rebuilt(nodes, o_s) - ref)))
    e_r = float(np.max(np.abs(rebuilt(nodes, o_r) - ref)))
    same = partition_key(nodes, o_s) == partition_key(nodes, o_r)
    return e_s, e_r, same, min_gap_deg(nodes), min_gap_abs(nodes)


def sweep(title: str, cases: list[tuple]) -> None:
    """两个规则对称报告：各自失败数、最大残差、最小排序间隔，外加分区是否分歧。"""
    print(f"\n=== {title}（{len(cases)} 例）===")
    nbad = {"arg(s)": 0, "|p|": 0}
    ndiff = 0
    worst = {"arg(s)": 0.0, "|p|": 0.0}
    tight = {"arg(s)": (180.0, None), "|p|": (1.0, None)}
    for kind, order, cut, rp, rs in cases:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                e_s, e_r, same, gap_s, gap_r = run_case(kind, order, cut, rp, rs)
        except Exception as exc:  # noqa: BLE001 - 极端规格可能直接设计失败
            nbad["arg(s)"] += 1
            nbad["|p|"] += 1
            print(f"  {kind:7s} N={order:<3d} cutoff={cut:<8g} rp={rp:<6g} rs={rs:<6g}  "
                  f"[设计失败: {type(exc).__name__}: {exc}]")
            continue
        for name, e, gap in (("arg(s)", e_s, gap_s), ("|p|", e_r, gap_r)):
            worst[name] = max(worst[name], e)
            if gap < tight[name][0]:
                tight[name] = (gap, (kind, order, cut, rp, rs))
        if e_s > TOL:
            nbad["arg(s)"] += 1
        if e_r > TOL:
            nbad["|p|"] += 1
        if e_s > TOL or e_r > TOL or not same:
            ndiff += 0 if same else 1
            print(f"  {kind:7s} N={order:<3d} cutoff={cut:<8g} rp={rp:<6g} rs={rs:<6g}  "
                  f"arg(s)={e_s:.3e}  |p|={e_r:.3e}  分区{'相同' if same else '不同'}  "
                  f"间隔 arg(s)={gap_s:.4g} deg / |p|={gap_r:.4g}")
    print(f"  失败：arg(s) {nbad['arg(s)']} 例，|p| {nbad['|p|']} 例；分区不同 {ndiff} 例")
    print(f"  最大残差：arg(s) {worst['arg(s)']:.3e}，|p| {worst['|p|']:.3e}")
    print(f"  最小间隔：arg(s) {tight['arg(s)'][0]:.4g} deg @ {tight['arg(s)'][1]}")
    print(f"             |p| {tight['|p|'][0]:.4g} @ {tight['|p|'][1]}")


def controls(cases_by_name: list[tuple[str, list[tuple]]]) -> None:
    """阳性对照：已知会坏的规则在同一套件上应该被大量报出来。"""
    rules = (("digital arg", lambda q: abs(np.angle(q))),)
    suites = [cases_by_name[0], cases_by_name[1], cases_by_name[2]]  # A/B/C 足够
    for name, key in rules:
        bad = collections.Counter()
        tot = collections.Counter()
        for _, cases in suites:
            for kind, order, cut, rp, rs in cases:
                tot[kind] += 1
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        z, p, k, nodes = build(kind, order, cut, rp, rs)
                    ref = sg.freqz_zpk(z, p, k, worN=w)[1]
                    o = split_order(nodes, key)
                    if np.max(np.abs(rebuilt(nodes, o) - ref)) > TOL:
                        bad[kind] += 1
                except Exception:  # noqa: BLE001
                    pass
        total = sum(tot.values())
        n_bad = sum(bad.values())
        detail = "  ".join(f"{k}:{bad[k]}/{tot[k]}" for k in ("butter", "cheby1", "cheby2", "ellip"))
        print(f"\n=== 阳性对照 {name}：失败 {n_bad}/{total} 例 ===")
        print(f"  {detail}")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="arg(s) 规则的极端压力测试")
    ap.add_argument("--controls", action="store_true",
                    help="额外跑阳性对照（数字域 arg(p)）")
    args = ap.parse_args(argv)

    cuts_a = [5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.15, 0.2, 0.25, 0.3,
              0.35, 0.4, 0.45, 0.49, 0.499, 0.4999]
    suite_a = [(k, n, c, 0.2, 40.0)
               for k in ("butter", "cheby1", "cheby2", "ellip")
               for n in range(3, 32, 2)
               for c in cuts_a]

    suite_b = [(k, n, c, rp, rs)
               for k in ("butter", "cheby1", "cheby2", "ellip")
               for rp in (0.001, 0.01, 0.1, 1.0, 3.0)
               for rs in (20.0, 40.0, 80.0, 120.0, 160.0)
               if rs > rp + 1e-9
               for n in (3, 7, 11, 15, 21)
               for c in (0.05, 0.2, 0.45)]

    suite_c = [(k, n, c, rp, rs)
               for k in ("butter", "cheby1", "cheby2", "ellip")
               for n in (31, 41, 51, 61, 71)
               for c in (0.1, 0.3, 0.45)
               for rp, rs in ((0.2, 40.0), (1.0, 120.0))]

    suite_d = [(k, n, c, 0.2, 40.0)
               for k in ("butter", "cheby1", "cheby2", "ellip")
               for n in (3, 7, 15, 31, 51)
               for c in (1e-4, 2e-4, 5e-4, 0.4999, 0.49999, 0.499999)]

    print(f"TOL = {TOL:g}   freqz 点数 = {NUMW}")
    sweep("A. 细分截止频率（orders 3..31 × 18 个 cutoff，rp=0.2/rs=40）", suite_a)
    sweep("B. 极端规格（rp × rs × orders × cutoff）", suite_b)
    sweep("C. 超高阶（orders 31..71）", suite_c)
    sweep("D. 贴边截止（1e-4 … 0.499999）", suite_d)

    cuts_e = np.linspace(1e-3, 0.4999, 400)
    suite_e = [(k, n, c, 0.2, 40.0)
               for k in ("butter", "cheby1", "cheby2", "ellip")
               for n in (7, 31)
               for c in cuts_e]
    sweep("E. 细扫 cutoff（400 个点 × 4 族 × orders 7/31）", suite_e)

    if args.controls:
        controls([("A", suite_a), ("B", suite_b), ("C", suite_c)])


if __name__ == "__main__":
    main()
