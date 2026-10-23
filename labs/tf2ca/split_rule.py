"""极点分组规则的对照与压力测试。

三种排序：

    |p|       数字极点模长升序交替（四族经典设计的工程规则）
    arg(s)    双线性反变换回 s 平面后的极点角（与 |p| 分区相同，允许两支互换）
    arg(p)    数字域极点角（反例：cheby2 整族失败）

单次运行画三种规则的极点分配与重建幅度；``--stress`` 跑四族多规格扫描，
``--controls`` 加上已知会坏的 ``arg(p)`` 作阳性对照。

用法::

    python split_rule.py
    python split_rule.py --kind cheby1 --order 9
    python split_rule.py --stress --controls
"""
from __future__ import annotations

import argparse
import collections
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 无界面后端
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal as sg

from tf2ca import design, reconstruct, response, sections, split_by_radius, tf_response

SCRIPT_DIR = Path(__file__).resolve().parent

COLOR = ("#1f77b4", "#d62728")  # 链 0 / 链 1
FS = 2.0  # 数字域 Nyquist = 1 的 scipy 约定


# ------------------------------------------------------------
# 三种排序规则
# ------------------------------------------------------------


def to_s(p: complex) -> complex:
    """双线性反变换：数字极点 -> s 平面极点（``2*fs`` 是正实比例，不影响 arg）。"""
    return 2.0 * FS * (p - 1.0) / (p + 1.0)


def node_arg_digital(s: dict) -> float:
    """节点的数字幅角：共轭对取上半平面那一个，实极点取自身。"""
    return float(np.angle(s["poles"][0]))


def node_arg_s(s: dict) -> float:
    """节点的 s 平面幅角。"""
    return float(np.angle(to_s(s["poles"][0])))


def split_by_key(a: np.ndarray, key) -> tuple[np.ndarray, np.ndarray]:
    """按 ``key(node)`` 升序交替分配节点，返回 (D0, D1)。"""
    nodes = sections(a)
    nodes.sort(key=key)
    d = [np.array([1.0]), np.array([1.0])]
    for k, node in enumerate(nodes):
        d[k % 2] = np.convolve(d[k % 2], node["poly"])
    return d[0], d[1]


def split_by_angle(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """按数字域 ``|arg(p)|`` 升序交替。"""
    return split_by_key(a, lambda s: abs(node_arg_digital(s)))


def split_by_s_angle(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """按 s 平面 ``arg(s)`` 升序交替。"""
    return split_by_key(a, node_arg_s)


def _poleset(d: np.ndarray) -> tuple:
    return tuple(np.round(np.sort_complex(np.roots(d)), 8))


def same_partition(d0a, d1a, d0b, d1b) -> bool:
    """两条链的极点集合是否相同（允许两条链互换）。"""
    return (_poleset(d0a) == _poleset(d0b) and _poleset(d1a) == _poleset(d1b)) or (
        _poleset(d0a) == _poleset(d1b) and _poleset(d1a) == _poleset(d0b)
    )


# ------------------------------------------------------------
# 单次运行：打印与绘图
# ------------------------------------------------------------


def dump_nodes(a: np.ndarray) -> None:
    """打印节点在三种排序下的名次与链归属。"""
    nodes = sections(a)
    by_r = sorted(nodes, key=lambda s: abs(s["poles"][0]))
    by_s = sorted(nodes, key=node_arg_s)
    by_a = sorted(nodes, key=lambda s: abs(node_arg_digital(s)))
    rank_r = {id(s): k for k, s in enumerate(by_r)}
    rank_s = {id(s): k for k, s in enumerate(by_s)}
    rank_a = {id(s): k for k, s in enumerate(by_a)}
    print("  节点表（每个节点 = 1 个实极点或 1 个共轭对）")
    print(f"  {'#|p|':>5} {'#arg(s)':>7} {'#arg(p)':>7} {'|p|':>8} {'arg(p)':>9} "
          f"{'arg(s)':>9}   chains(|p|,arg(s),arg(p))  poles")
    for s in by_r:
        p0 = s["poles"][0]
        print(f"  {rank_r[id(s)]:>5} {rank_s[id(s)]:>7} {rank_a[id(s)]:>7} "
              f"{abs(p0):>8.5f} {np.degrees(node_arg_digital(s)):>9.3f} "
              f"{np.degrees(node_arg_s(s)):>9.3f}   "
              f"{rank_r[id(s)] % 2} {rank_s[id(s)] % 2} {rank_a[id(s)] % 2}  "
              f"{np.round(np.array(s['poles']), 6).tolist()}")


def report(b, a, w, name: str, d0, d1) -> None:
    """打印一种规则的重建结果。"""
    ref = tf_response(b, a, w)
    got, _ = response(d0, d1, w)
    brec, _, _ = reconstruct(d0, d1)
    z_rec = np.roots(brec) if brec.size > 1 else np.array([])
    dev = np.max(np.abs(np.abs(z_rec) - 1.0)) if z_rec.size else 0.0
    print(f"  {name:>10}: max|dH| = {np.max(np.abs(got - ref)):.3e}   "
          f"重建零点 max|1-|z|| = {dev:.3e}")


def plot_poles(ax, a, b, rule: str, label: str) -> None:
    """z 平面：极点按链着色 + 名次标注。"""
    theta = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(theta), np.sin(theta), "k:", lw=0.8)

    key = {"|p|": lambda s: abs(s["poles"][0]),
           "arg(s)": node_arg_s,
           "arg(p)": lambda s: abs(node_arg_digital(s))}[rule]
    nodes = sorted(sections(a), key=key)
    for k, s in enumerate(nodes):
        for p in s["poles"]:
            ax.plot(p.real, p.imag, "x", color=COLOR[k % 2], ms=10, mew=2.2)
        p0 = s["poles"][0]
        txt = str(k) if rule == "|p|" else f"{k}\n{np.degrees(key(s)):.1f}deg"
        ax.annotate(txt, (p0.real, p0.imag), textcoords="offset points",
                    xytext=(7, 7), fontsize=8, color=COLOR[k % 2], weight="bold")

    z_design = np.roots(b)
    ax.plot(z_design.real, z_design.imag, "o", mfc="none", mec="#2ca02c", ms=8,
            label="design zeros B")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set(title=f"{label}: split by {rule}", xlabel="Re", ylabel="Im",
           xlim=(-1.25, 1.25), ylim=(-1.25, 1.25), aspect="equal")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=7)


def plot_mag(ax, b, a, w, rule: str, d0, d1, label: str) -> None:
    """幅度响应：原设计 vs 该规则的重建。"""
    ref_c = tf_response(b, a, w)
    got_c, _ = response(d0, d1, w)
    err = float(np.max(np.abs(got_c - ref_c)))
    ax.plot(w / np.pi, 20 * np.log10(np.abs(ref_c) + 1e-18), "k", lw=1.6, label="design B/A")
    ax.plot(w / np.pi, 20 * np.log10(np.abs(got_c) + 1e-18), color=COLOR[1], ls="--",
            lw=1.3, label=f"rebuilt ({rule})")
    ax.set(title=f"{label}: {rule}  (max|dH|={err:.2e})",
           xlabel="$\\omega/\\pi$", ylabel="dB", ylim=(-100, 3))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)


def plot_compare(b, a, results: dict, out: Path, label: str, num: int = 4000) -> None:
    """2x3：三种规则的极点分配与重建幅度。"""
    w = np.linspace(0, np.pi, num)
    fig, ax = plt.subplots(2, 3, figsize=(16.5, 9.0))
    for col, rule in enumerate(("|p|", "arg(s)", "arg(p)")):
        d0, d1 = results[rule]
        plot_poles(ax[0, col], a, b, rule, label)
        plot_mag(ax[1, col], b, a, w, rule, d0, d1, label)
    fig.suptitle("odd filter: pole assignment  |p| vs arg(s) vs arg(p)", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"图已保存：{out}")


# ------------------------------------------------------------
# 压力测试：zpk 路径，避免多项式求根的条件数
# ------------------------------------------------------------


def design_zpk(kind: str, order: int, cut: float, rp: float, rs: float):
    """scipy 直接给数字滤波器的 zpk（``cut`` 以 Nyquist = 1 归一）。"""
    return sg.iirfilter(order, cut, rp, rs, btype="low", ftype=kind, output="zpk")


def nodes_from_zpk(p: np.ndarray) -> list[complex]:
    """每个共轭对取上半平面那一个代表，实极点单独一个。"""
    return [q for q in np.asarray(p) if q.imag >= 0.0]


def split_order(nodes: list[complex], key) -> list[int]:
    """按 ``key(node)`` 升序给出名次排列（名次奇偶 = 归哪条链）。"""
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


def chain_response(nodes: list[complex], order: list[int], chain: int, w) -> np.ndarray:
    """全通链响应（零点 = 极点共轭倒数，按 z=1 归一到 A(1)=1）。"""
    z, p = chain_zpk(nodes, order, chain)
    h = sg.freqz_zpk(z, p, 1.0, worN=w)[1]
    h0 = sg.freqz_zpk(z, p, 1.0, worN=[0.0])[1][0]
    return h / h0


def rebuilt(nodes: list[complex], order: list[int], w) -> np.ndarray:
    return (chain_response(nodes, order, 0, w) + chain_response(nodes, order, 1, w)) / 2.0


def partition_key(nodes: list[complex], order: list[int]) -> tuple:
    chains = []
    for chain in (0, 1):
        _, ps = chain_zpk(nodes, order, chain)
        chains.append(tuple(np.round(np.sort_complex(ps), 9)))
    return tuple(sorted(chains))


def min_gap_deg(nodes: list[complex]) -> float:
    """``arg(s)`` 排序后相邻间隔的最小值（度）。"""
    a = np.sort([np.degrees(np.angle(to_s(q))) for q in nodes])
    return float(np.min(np.diff(a))) if len(a) > 1 else 180.0


def min_gap_abs(nodes: list[complex]) -> float:
    """``|p|`` 排序后相邻间隔的最小值。"""
    a = np.sort([abs(q) for q in nodes])
    return float(np.min(np.diff(a))) if len(a) > 1 else 1.0


def run_case(kind, order, cut, rp, rs, w):
    """同一设计上比较 arg(s) 与 |p|：残差、分区是否一致、排序裕度。"""
    z, p, k = design_zpk(kind, order, cut, rp, rs)
    nodes = nodes_from_zpk(p)
    ref = sg.freqz_zpk(z, p, k, worN=w)[1]
    o_s = split_order(nodes, lambda q: np.angle(to_s(q)))
    o_r = split_order(nodes, lambda q: abs(q))
    e_s = float(np.max(np.abs(rebuilt(nodes, o_s, w) - ref)))
    e_r = float(np.max(np.abs(rebuilt(nodes, o_r, w) - ref)))
    same = partition_key(nodes, o_s) == partition_key(nodes, o_r)
    return e_s, e_r, same, min_gap_deg(nodes), min_gap_abs(nodes)


def sweep(title: str, cases: list[tuple], w, tol: float) -> None:
    """报告两种规则的失败数、最大残差与最小排序间隔。"""
    print(f"\n=== {title}（{len(cases)} 例）===")
    nbad = {"arg(s)": 0, "|p|": 0}
    ndiff = 0
    worst = {"arg(s)": 0.0, "|p|": 0.0}
    tight = {"arg(s)": (180.0, None), "|p|": (1.0, None)}
    for kind, order, cut, rp, rs in cases:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                e_s, e_r, same, gap_s, gap_r = run_case(kind, order, cut, rp, rs, w)
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
        nbad["arg(s)"] += e_s > tol
        nbad["|p|"] += e_r > tol
        if not same:
            ndiff += 1
        if e_s > tol or e_r > tol or not same:
            print(f"  {kind:7s} N={order:<3d} cutoff={cut:<8g} rp={rp:<6g} rs={rs:<6g}  "
                  f"arg(s)={e_s:.3e}  |p|={e_r:.3e}  分区{'相同' if same else '不同'}  "
                  f"间隔 arg(s)={gap_s:.4g} deg / |p|={gap_r:.4g}")
    print(f"  失败：arg(s) {nbad['arg(s)']} 例，|p| {nbad['|p|']} 例；分区不同 {ndiff} 例")
    print(f"  最大残差：arg(s) {worst['arg(s)']:.3e}，|p| {worst['|p|']:.3e}")
    print(f"  最小间隔：arg(s) {tight['arg(s)'][0]:.4g} deg @ {tight['arg(s)'][1]}")
    print(f"             |p| {tight['|p|'][0]:.4g} @ {tight['|p|'][1]}")


def controls(cases_by_name: list[tuple[str, list[tuple]]], w, tol: float) -> None:
    """阳性对照：已知会坏的 ``arg(p)`` 在同一套件上应被大量报出。"""
    bad = collections.Counter()
    tot = collections.Counter()
    for _, cases in cases_by_name:
        for kind, order, cut, rp, rs in cases:
            tot[kind] += 1
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    z, p, k = design_zpk(kind, order, cut, rp, rs)
                ref = sg.freqz_zpk(z, p, k, worN=w)[1]
                nodes = nodes_from_zpk(p)
                order_d = split_order(nodes, lambda q: abs(np.angle(q)))
                if np.max(np.abs(rebuilt(nodes, order_d, w) - ref)) > tol:
                    bad[kind] += 1
            except Exception:  # noqa: BLE001
                pass
    total = sum(tot.values())
    n_bad = sum(bad.values())
    detail = "  ".join(f"{k}:{bad[k]}/{tot[k]}" for k in ("butter", "cheby1", "cheby2", "ellip"))
    print(f"\n=== 阳性对照 arg(p)：失败 {n_bad}/{total} 例 ===")
    print(f"  {detail}")


def run_stress(tol: float, with_controls: bool) -> None:
    """四组批量扫描（A-E），比较 arg(s) 与 |p|。"""
    w = np.linspace(0, np.pi, 1200)
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
    cuts_e = np.linspace(1e-3, 0.4999, 400)
    suite_e = [(k, n, c, 0.2, 40.0)
               for k in ("butter", "cheby1", "cheby2", "ellip")
               for n in (7, 31)
               for c in cuts_e]

    print(f"TOL = {tol:g}   freqz 点数 = {len(w)}")
    sweep("A. 细分截止频率（orders 3..31 × 18 个 cutoff，rp=0.2/rs=40）", suite_a, w, tol)
    sweep("B. 极端规格（rp × rs × orders × cutoff）", suite_b, w, tol)
    sweep("C. 超高阶（orders 31..71）", suite_c, w, tol)
    sweep("D. 贴边截止（1e-4 … 0.499999）", suite_d, w, tol)
    sweep("E. 细扫 cutoff（400 个点 × 4 族 × orders 7/31）", suite_e, w, tol)
    if with_controls:
        controls([("A", suite_a), ("B", suite_b), ("C", suite_c)], w, tol)


# ------------------------------------------------------------
# main
# ------------------------------------------------------------


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="pole assignment: |p| vs arg(s) vs arg(p)")
    p.add_argument("--kind", default="cheby2", choices=["butter", "cheby1", "cheby2", "ellip"])
    p.add_argument("--order", type=int, default=9, help="奇数阶")
    p.add_argument("--cutoff", type=float, default=0.3, help="Nyquist 归一截止频率 (0..1)")
    p.add_argument("--ripple", type=float, default=0.2, help="通带波纹 dB (cheby1/ellip)")
    p.add_argument("--stop", type=float, default=40.0, help="阻带衰减 dB (cheby2/ellip)")
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--stress", action="store_true", help="跑四族多规格压力扫描")
    p.add_argument("--controls", action="store_true", help="压力扫描时加 arg(p) 阳性对照")
    p.add_argument("--tol", type=float, default=1e-8, help="压力扫描的残差容差")
    args = p.parse_args(argv)

    if args.stress:
        run_stress(args.tol, args.controls)
        return

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
    results = {
        "|p|": split_by_radius(a),
        "arg(s)": split_by_s_angle(a),
        "arg(p)": split_by_angle(a),
    }
    print("  重建校验")
    for name, dd in results.items():
        report(b, a, w, name, *dd)
        print(f"                  -> 与 |p| 分区"
              f"{'相同' if same_partition(*dd, *results['|p|']) else '不同'}")
    print()
    plot_compare(b, a, results, out, f"{args.kind} N={args.order}")


if __name__ == "__main__":
    main()
