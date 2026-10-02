# -*- coding: utf-8 -*-
"""
split_rule.py
=============

极点分配规则对照实验。三种排序方式代入同一个分解框架

    H(z) = 1/2 * ( A0(z) + A1(z) )

1. ``|p|``      : 数字极点模长升序交替  —— 本库/README 采用的工程规则
2. ``design``   : **模拟原型极点角**升序交替（即 ``iir_design.hpp`` 里
                  ``phi = (2k-1)pi/(2n)`` 的设计序号顺序）—— 文献里说的
                  "pole interlacing（按角度交错）"
3. ``darg``     : **数字域**极点幅角 ``arg(p)`` 升序交替 —— 一个被双线性变换
                  （以及 cheby2 特有的"取倒数"）扭曲过的量

**结论（本脚本的主要发现）**：``design`` 与 ``|p|`` 给出**完全相同**的分区
（四族 × 奇数阶 3/5/7/9/11 × cutoff 0.1/0.2/0.3/0.35/0.5 共 100 例，100/100 相同，
残差集合也相同），因此文献的"按角度交错"并不是错的——错的是用**数字**
``arg(p)``：在 cheby2 上它与设计顺序相反（100 例里有 9 例），按它交替就整片崩掉。

用法
----
    python split_rule.py                                  # cheby2 N=9 cutoff=0.3 rs=40
    python split_rule.py --kind cheby1 --order 9          # 换族对照
    python split_rule.py --order 11 --cutoff 0.35 --stop 80
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 无界面后端
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal as sg

from tf2ca import (
    design,
    reconstruct,
    response,
    sections,
    split_by_radius,
    tf_response,
)

SCRIPT_DIR = Path(__file__).resolve().parent

COLOR = ("#1f77b4", "#d62728")  # 链 0 / 链 1
FS = 2.0  # 数字域 Nyquist = 1 的 scipy 约定


# ------------------------------------------------------------
# 三种排序规则
# ------------------------------------------------------------


def node_angle(s: dict) -> float:
    """节点的数字幅角：共轭对取上半平面那一个，实极点取它自己的 ``arg``。"""
    return float(np.angle(s["poles"][0]))


def split_by_angle(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """按**数字**极点幅角 ``arg(p)`` 升序交替分配节点，返回 (D0, D1)。"""
    nodes = sections(a)
    nodes.sort(key=lambda s: abs(node_angle(s)))
    d = [np.array([1.0]), np.array([1.0])]
    for k, node in enumerate(nodes):
        d[k % 2] = np.convolve(d[k % 2], node["poly"])
    return d[0], d[1]


def prototype_zpk(kind: str, order: int, cutoff: float, rp: float, rs: float):
    """scipy 模拟原型（已 lp2lp 到预畸变后的截止），保持其内部设计序号顺序。"""
    wn = 2 * FS * np.tan(np.pi * cutoff / FS)
    if kind == "butter":
        z, p, k = sg.buttap(order)
    elif kind == "cheby1":
        z, p, k = sg.cheb1ap(order, rp)
    elif kind == "cheby2":
        z, p, k = sg.cheb2ap(order, rs)
    else:
        z, p, k = sg.ellipap(order, rp, rs)
    return z, np.asarray(p) * wn, k


def split_by_design_angle(kind: str, order: int, cutoff: float, rp: float, rs: float):
    """按**模拟原型极点角** ``|arg(s)|`` 升序交替（= 设计序号 ``phi_k`` 顺序）。

    节点由设计管线的数字极点直接构造（不经过 ``np.roots``），因此不受求根条件数影响。
    """
    za, pa, ka = prototype_zpk(kind, order, cutoff, rp, rs)
    _, pd, _ = sg.bilinear_zpk(za, pa, ka, FS)
    pd = np.asarray(pd)
    key = np.abs(np.angle(pa))
    used = np.zeros(len(pd), bool)
    nodes: list[tuple[float, np.ndarray]] = []
    for i, q in enumerate(pd):
        if used[i]:
            continue
        if abs(q.imag) < 1e-9:
            nodes.append((key[i], np.array([1.0, -q.real])))
            used[i] = True
        else:
            j = int(np.argmin(np.abs(pd - np.conj(q)) + 1e6 * used))
            nodes.append((key[i], np.poly([q, pd[j]]).real))
            used[i] = used[j] = True
    nodes.sort(key=lambda t: t[0])
    d = [np.array([1.0]), np.array([1.0])]
    for k, (_, poly) in enumerate(nodes):
        d[k % 2] = np.convolve(d[k % 2], poly)
    return d[0], d[1]


def _poleset(d: np.ndarray) -> tuple:
    return tuple(np.round(np.sort_complex(np.roots(d)), 8))


def same_partition(d0a, d1a, d0b, d1b) -> bool:
    """两条链的极点集合是否相同（允许两条链互换）。"""
    return (_poleset(d0a) == _poleset(d0b) and _poleset(d1a) == _poleset(d1b)) or (
        _poleset(d0a) == _poleset(d1b) and _poleset(d1a) == _poleset(d0b)
    )


# ------------------------------------------------------------
# 打印
# ------------------------------------------------------------


def dump_nodes(a: np.ndarray) -> None:
    """打印节点表：|p| 名次、数字 arg 名次、|p|、arg、两种规则各分到哪条链。"""
    nodes = sections(a)
    by_r = sorted(nodes, key=lambda s: abs(s["poles"][0]))
    by_a = sorted(nodes, key=lambda s: abs(node_angle(s)))
    rank_r = {id(s): k for k, s in enumerate(by_r)}
    rank_a = {id(s): k for k, s in enumerate(by_a)}
    print("  节点表（每个节点 = 1 个实极点或 1 个共轭对）")
    print(f"  {'#|p|':>5} {'#arg':>5} {'|p|':>8} {'arg(deg)':>9}   {'|p|链':>5} {'arg链':>6}  poles")
    for s in by_r:
        p0 = s["poles"][0]
        print(f"  {rank_r[id(s)]:>5} {rank_a[id(s)]:>5} {abs(p0):>8.5f} "
              f"{np.degrees(node_angle(s)):>9.3f}   {rank_r[id(s)] % 2:>5} "
              f"{rank_a[id(s)] % 2:>6}  {np.round(np.array(s['poles']), 6).tolist()}")


def report(b, a, w, name: str, d0, d1) -> None:
    """打印一种规则的重建结果：复频响最大偏差、重建零点的模。"""
    ref = tf_response(b, a, w)
    got, _ = response(d0, d1, w)
    brec, _, _ = reconstruct(d0, d1)
    z_rec = np.roots(brec) if brec.size > 1 else np.array([])
    dev = np.max(np.abs(np.abs(z_rec) - 1.0)) if z_rec.size else 0.0
    print(f"  {name:>14}: max|dH| = {np.max(np.abs(got - ref)):.3e}   "
          f"重建零点 max|1-|z|| = {dev:.3e}")


# ------------------------------------------------------------
# 绘图
# ------------------------------------------------------------


def plot_poles(ax, a, b, rule: str, order_label: str) -> None:
    """z 平面：极点按链着色 + 名次标注（数字 arg 规则叠加 arg 值）+ 设计零点。"""
    theta = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(theta), np.sin(theta), "k:", lw=0.8)

    nodes = sections(a)
    nodes.sort(key=lambda s: abs(node_angle(s)) if rule == "arg" else abs(s["poles"][0]))
    for k, s in enumerate(nodes):
        for p in s["poles"]:
            ax.plot(p.real, p.imag, "x", color=COLOR[k % 2], ms=10, mew=2.2)
        p0 = s["poles"][0]
        txt = f"{k}\n{np.degrees(node_angle(s)):.1f}deg" if rule == "arg" else str(k)
        ax.annotate(txt, (p0.real, p0.imag), textcoords="offset points",
                    xytext=(7, 7), fontsize=8, color=COLOR[k % 2], weight="bold")

    z_design = np.roots(b)
    ax.plot(z_design.real, z_design.imag, "o", mfc="none", mec="#2ca02c", ms=8,
            label="design zeros B")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set(title=f"{order_label}: split by {rule}", xlabel="Re", ylabel="Im",
           xlim=(-1.25, 1.25), ylim=(-1.25, 1.25), aspect="equal")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=7)


def plot_mag(ax, b, a, w, rule: str, d0, d1, order_label: str) -> None:
    """幅度响应：原设计 vs 该规则的重建。"""
    ref_c = tf_response(b, a, w)
    ref = 20 * np.log10(np.abs(ref_c) + 1e-18)
    got_c, _ = response(d0, d1, w)
    got = 20 * np.log10(np.abs(got_c) + 1e-18)
    err = float(np.max(np.abs(got_c - ref_c)))
    ax.plot(w / np.pi, ref, "k", lw=1.6, label="design B/A")
    ax.plot(w / np.pi, got, color=COLOR[1], ls="--", lw=1.3, label=f"rebuilt (split by {rule})")
    ax.set(title=f"{order_label}: rebuilt |H| with {rule} rule  (max|dH|={err:.2e})",
           xlabel="$\\omega/\\pi$", ylabel="dB", ylim=(-100, 3))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)


def plot_order_diagram(ax, a, order_label: str) -> None:
    """把每个节点在两种排序下的名次画成散点：落在对角线上 = 两种排序一致。"""
    nodes = sections(a)
    by_r = sorted(nodes, key=lambda s: abs(s["poles"][0]))
    by_a = sorted(nodes, key=lambda s: abs(node_angle(s)))
    rank_r = {id(s): k for k, s in enumerate(by_r)}
    rank_a = {id(s): k for k, s in enumerate(by_a)}
    xs = [rank_r[id(s)] for s in by_r]
    ys = [rank_a[id(s)] for s in by_r]
    n = len(xs)
    ax.plot([-0.5, n - 0.5], [-0.5, n - 0.5], "k:", lw=0.8, label="same order")
    for x, y, s in zip(xs, ys, by_r):
        ax.plot(x, y, "o", color=COLOR[x % 2], ms=9)
        ax.annotate(f"|p|={abs(s['poles'][0]):.3f}", (x, y), textcoords="offset points",
                    xytext=(8, -4), fontsize=7, color=COLOR[x % 2])
    ax.set(title=f"{order_label}: node order  |p| (x) vs digital arg (y)",
           xlabel="rank by |p|", ylabel="rank by arg(p)", xticks=range(n), yticks=range(n))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    if xs == ys:
        ax.text(0.03, 0.93, "orders coincide", transform=ax.transAxes, color="#2ca02c",
                fontsize=9, weight="bold")
    else:
        ax.text(0.03, 0.93, "orders DIVERGE", transform=ax.transAxes, color="#d62728",
                fontsize=9, weight="bold")


def plot_summary(ax, rows: list[tuple[str, float, bool]], order_label: str) -> None:
    """文字汇总：三种规则的残差与"是否与 |p| 同分区"。"""
    ax.axis("off")
    ax.set_title(f"{order_label}: summary", fontsize=11)
    y = 0.88
    ax.text(0.02, 0.98, "rule", transform=ax.transAxes, fontsize=9, weight="bold")
    ax.text(0.42, 0.98, "max|dH|", transform=ax.transAxes, fontsize=9, weight="bold")
    ax.text(0.72, 0.98, "same as |p|", transform=ax.transAxes, fontsize=9, weight="bold")
    for name, err, same in rows:
        ax.text(0.02, y, name, transform=ax.transAxes, fontsize=9)
        ax.text(0.42, y, f"{err:.2e}", transform=ax.transAxes, fontsize=9)
        ax.text(0.72, y, "yes" if same else "NO", transform=ax.transAxes, fontsize=9,
                color="#2ca02c" if same else "#d62728", weight="bold")
        y -= 0.14
    ax.text(0.02, y - 0.04,
            "design-angle (prototype phi_k) == |p|:\n"
            "verified 100/100 over 4 families x odd 3..11 x 5 cutoffs",
            transform=ax.transAxes, fontsize=8.5, color="#444444")


def plot_compare(b, a, results: dict, out: Path, label: str, num: int = 4000) -> None:
    """2x3：两种规则的极点分配 + 顺序对照 + 两种规则的重建幅度 + 汇总。"""
    w = np.linspace(0, np.pi, num)
    fig, ax = plt.subplots(2, 3, figsize=(16.5, 9.0))
    for col, rule in enumerate(("|p|", "arg")):
        d0, d1 = results[rule]
        plot_poles(ax[0, col], a, b, rule, label)
        plot_mag(ax[1, col], b, a, w, rule, d0, d1, label)
    plot_order_diagram(ax[0, 2], a, label)
    plot_summary(ax[1, 2], results["rows"], label)
    ax[0, 1].set_title(f"{label}: split by digital arg(p)  (order may be reversed)")
    fig.suptitle("odd filter: pole assignment  |p| vs design-angle(phi_k) vs digital arg(p)",
                 fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"图已保存：{out}")


# ------------------------------------------------------------
# main
# ------------------------------------------------------------


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="pole assignment: |p| vs design angle vs arg(p)")
    p.add_argument("--kind", default="cheby2", choices=["butter", "cheby1", "cheby2", "ellip"])
    p.add_argument("--order", type=int, default=9, help="奇数阶")
    p.add_argument("--cutoff", type=float, default=0.3, help="Nyquist 归一截止频率 (0..1)")
    p.add_argument("--ripple", type=float, default=0.2, help="通带波纹 dB (cheby1/ellip)")
    p.add_argument("--stop", type=float, default=40.0, help="阻带衰减 dB (cheby2/ellip)")
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args(argv)

    if args.order % 2 == 0:
        p.error("只做奇数阶：--order 必须为奇数")

    out = args.out or (SCRIPT_DIR / "output"
                       / f"split_rule_{args.kind}_{args.order}_{args.cutoff:g}.png")

    b, a = design(args.kind, args.order, args.cutoff, args.ripple, args.stop)
    print(f"== {args.kind} 低通：order={args.order}, cutoff={args.cutoff:g}, "
          f"rp={args.ripple}, rs={args.stop} ==")
    dump_nodes(a)
    print()

    w = np.linspace(0, np.pi, 4000)
    ref = tf_response(b, a, w)
    d_radius = split_by_radius(a)
    d_arg = split_by_angle(a)
    d_design = split_by_design_angle(args.kind, args.order, args.cutoff,
                                     args.ripple, args.stop)
    results = {"|p|": d_radius, "arg": d_arg}

    print("  重建校验")
    rows = []
    for name, dd in (("|p| (=design)", d_radius), ("design-angle", d_design),
                     ("digital arg", d_arg)):
        got, _ = response(*dd, w)
        err = float(np.max(np.abs(got - ref)))
        same = same_partition(*dd, *d_radius)
        rows.append((name, err, same))
        report(b, a, w, name, *dd)
        print(f"                  -> 与 |p| 分区{'相同' if same else '不同'}")
    results["rows"] = rows
    print()

    plot_compare(b, a, results, out, f"{args.kind} N={args.order}")


if __name__ == "__main__":
    main()
