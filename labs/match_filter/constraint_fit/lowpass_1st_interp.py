"""一阶低通：5 点有理插值 + BRASIL 区间等化。

数字二阶节的功率函数在 ``x = sin^2(w/2)`` 上是 2/2 有理函数：

    P(x) = (n0 + n1 x + n2 x^2) / (1 + d1 x + d2 x^2)

五个未知数，所以在 5 个节点上令 ``P = Pa`` 就是解一个 5x5 线性方程组。
分母不再由额外约束“锁定”，而是插值本身的产物。

有理逼近的 Chebyshev 交替定理（零缺陷）给出：n 阶最佳逼近在 ``2n+1`` 个
节点上插值，误差在 ``2n+2`` 个点等幅交替。这里 n = 2，即 5 个节点、6 个
交替极值。BRASIL（Hofreither, RICAM 2020-37）把这条定理变成迭代算法：
交替调整 5 个节点把区间切成的 6 段长度，让每段的最大 |ΔdB| 相等。
它不动极点，只动插值节点；线性收敛，实测 600 轮足够。

注意首末区间必须把 ``x = 0`` 与 ``x = 1`` 计入误差最大值——目标在
``x -> 1`` 处有平方根型边界奇点，Nyquist 端的相对误差通常是全带最大。

结果（fs = 48 kHz，max|ΔdB|）：20 / 117.6 / 691.4 / 4065 / 23900 Hz 分别为
0.0795 / 0.0793 / 0.0777 / 0.0676 / 0.0259 dB。作为对照，固定节点
（Chebyshev 根）为 16.3 / 3.45 / 0.47 / 0.46 / 0.23 dB。

运行：``python qwqdsp/labs/match_filter/constraint_fit/lowpass_1st_interp.py``
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "best_pole_fit"))
import power_match as PM  # noqa: E402

PI = PM.PI
"""圆周率。"""

FS = PM.FS
"""实验采样率。"""

FC_LIST = PM.FC_LIST
"""与其余实验一致的 5 个截止频率。"""

CHEB_NODES = 0.5 * (1.0 - np.cos((2.0 * np.arange(5) + 1.0) * PI / 10.0))
"""5 个 Chebyshev 根节点（在 x 上），作为 BRASIL 的初值。"""


# --------------------------------------------------------------------- 插值
def _match_from_nodes(fc, xnodes, fs=FS):
    """在给定节点上精确插值；失败返回 ``(None, 原因)``。"""
    x = np.sort(np.asarray(xnodes, dtype=float))
    if x.size != 5 or np.min(np.diff(x)) <= 1e-12:
        return None, "nodes"
    wc = 2.0 * PI * fc / fs
    w = 2.0 * np.arcsin(np.sqrt(x))
    target = 1.0 / (1.0 + (w / wc) ** 2)
    M = np.stack([np.ones_like(x), x, x * x, -target * x, -target * x * x], axis=1)
    scale = np.linalg.norm(M, axis=0)
    try:
        n0, n1, n2, d1, d2 = np.linalg.solve(M / scale, target) / scale
    except np.linalg.LinAlgError:
        return None, "singular"

    # 由 x 基分母反推 (p, q)：d1 = 1 + r^2 + 16/p^2 - 8(1+r)/p - 1
    r2 = d1 + 1.0 + d2
    if not r2 > 0.0:
        return None, "r2 <= 0"
    r = np.sqrt(r2)
    a1 = d1 + 2.0
    A = a1 - 1.0 - r * r
    B = 8.0 * (1.0 + r)
    disc = B * B + 64.0 * A
    if disc < 0.0:
        return None, "no pole"
    roots = ([(-B + s * np.sqrt(disc)) / (2.0 * A) for s in (1.0, -1.0)]
             if abs(A) > 1e-15 else [16.0 / B])
    p_bilinear = 2.0 / (1.0 / np.tan(wc / 2.0) + 1.0)
    stable = [p for p in roots if p > 0.0 and p * (1.0 + r) < 4.0]
    if not stable:
        return None, "unstable"
    p = min(stable, key=lambda z: abs(z - p_bilinear))
    q = r * p

    a0 = p * p
    beta = a0 * np.array([n0, n0 + n1 + n2, -n2 / 4.0])
    if not PM.feasible(beta):
        return None, "infeasible"
    a = np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])
    return PM.Match(PM.solve_b(beta), a, float(p), float(q), beta, PM.alpha_of(p, q)), "ok"


def error_db(match, x, fc, fs=FS):
    """在 ``x = sin^2(w/2)`` 上的 ``ΔdB``。"""
    x = np.asarray(x, dtype=float)
    ph = np.stack([1.0 - x, x, 4.0 * x * (1.0 - x)], axis=-1)
    hd2 = (ph @ match.beta) / (ph @ match.alpha)
    w = 2.0 * np.arcsin(np.sqrt(x))
    ha2 = 1.0 / (1.0 + (w / (2.0 * PI * fc / fs)) ** 2)
    return 10.0 * np.log10(np.maximum(hd2, 1e-300) / np.maximum(ha2, 1e-300))


# --------------------------------------------------------------------- BRASIL
def design(fc, fs=FS, max_iter=600, samples=300, step_max=0.15, kappa=0.15):
    """BRASIL 区间等化；返回 ``(match, nodes, max|ΔdB|)``。"""
    nodes = CHEB_NODES.copy()
    best = (np.inf, None, None)
    for _ in range(max_iter):
        match, why = _match_from_nodes(fc, nodes, fs)
        if match is None:
            raise RuntimeError(f"fc = {fc:g} Hz：{why}")
        edges = np.concatenate(([0.0], nodes, [1.0]))
        errs = np.empty(6)
        for i in range(6):
            g = np.linspace(edges[i], edges[i + 1], samples + 2)
            if i > 0:
                g = g[1:]
            if i < 5:
                g = g[:-1]
            errs[i] = np.max(np.abs(error_db(match, g, fc, fs)))
        top = float(errs.max())
        if top < best[0]:
            best = (top, match, nodes.copy())
        mean = errs.mean()
        dev = errs - mean
        dmax = np.abs(dev).max()
        if dmax <= 0.0:
            break
        if errs.min() > 0.0 and top / errs.min() - 1.0 < 1e-4:
            break
        delta = dev / dmax
        step = min(step_max, kappa / (dmax / mean)) if mean > 0.0 else step_max
        lengths = (1.0 + step) ** (-delta) * np.diff(edges)
        lengths /= lengths.sum()
        nodes = np.concatenate(([0.0], np.cumsum(lengths)))[1:-1]

    _, match, nodes = best
    err = float(np.max(np.abs(error_db(match, np.linspace(0.0, 1.0, 200001), fc, fs))))
    return match, nodes, err


# --------------------------------------------------------------------- 实验
def main():
    print(f"fs = {FS:g} Hz；5 点有理插值 + BRASIL 区间等化")
    print(f"{'fc':>8}{'max|dB|':>10}   节点（x = sin^2(w/2)）")
    rows = []
    for fc in FC_LIST:
        match, nodes, err = design(fc)
        rows.append((fc, match))
        print(f"{fc:8.1f}{err:10.4f}   " + " ".join(f"{v:.5f}" for v in nodes))
    print()

    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    _figure(rows, out)
    print(f"图：{out / 'lowpass_1st_interp.png'}")


def _figure(rows, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    fig, (ax_r, ax_e) = plt.subplots(
        2, 1, figsize=(10.5, 7.4), gridspec_kw={"height_ratios": (2.0, 1.2),
                                                "hspace": 0.28})
    f = np.geomspace(10.0, FS / 2.0, 4001)
    w = 2.0 * PI * f / FS
    for j, (fc, match) in enumerate(rows):
        color = f"C{j}"
        hd2 = match.power(w)
        ha2 = 1.0 / (1.0 + (w / (2.0 * PI * fc / FS)) ** 2)
        ax_r.semilogx(f, 10.0 * np.log10(np.maximum(ha2, 1e-300)),
                      color=color, lw=3.0, ls="--", alpha=0.45)
        ax_r.semilogx(f, 10.0 * np.log10(np.maximum(hd2, 1e-300)),
                      color=color, lw=1.3, label=f"{fc:g} Hz")
        err = error_db(match, np.linspace(0.0, 1.0, 200001), fc)
        ww = 2.0 * np.arcsin(np.sqrt(np.linspace(0.0, 1.0, 200001)))
        ax_e.plot(ww / PI * FS / 2.0, np.maximum(np.abs(err), 1e-12),
                  color=color, lw=1.2)

    ax_r.set_xscale("log")
    ax_r.set_xlim(10.0, FS / 2.0)
    ax_r.set_ylim(-90.0, 5.0)
    ax_r.set_ylabel("幅度 / dB")
    ax_r.grid(True, which="both", alpha=0.25)
    ax_r.legend(fontsize=8, title="fc", ncol=2)

    ax_e.set_xlim(0.0, FS / 2.0)
    ax_e.set_yscale("log")
    ax_e.set_ylim(1e-6, 1e2)
    ax_e.set_xlabel("频率 / Hz")
    ax_e.set_ylabel("|ΔdB|")
    ax_e.grid(True, which="both", alpha=0.25)

    ax_r.set_title("一阶低通: 无硬约束+BRASIL 5点有理插值")
    fig.subplots_adjust(left=0.085, right=0.98, top=0.92, bottom=0.09, hspace=0.28)
    fig.savefig(out / "lowpass_1st_interp.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
