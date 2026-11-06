"""一阶高通：4 点有理插值 + BRASIL 区间等化。

一阶模拟高通的功率函数

    Pa(w) = r^2 / (1 + r^2),    r = w / wc

在 ``x = sin^2(w/2)`` 下仍是 2/2 有理函数，但高通在直流处必须为零：
``P(0) = 0``，所以分子只能从一次项开始，

    P(x) = x (n1 + n2 x) / (1 + d1 x + d2 x^2)。

只剩 4 个自由参数，因此在 4 个节点上令 ``P = Pa`` 就是 4x4 线性方程组。
按有理逼近的交替定理，4 参数的最佳逼近在 5 个点等幅交替；BRASIL 相应地把
区间切成 5 段做长度等化（4 节点、5 区间）。

直流处数字与模拟功率都为 0，dB 误差取极限：

    lim P/P a = n1 * wc^2 / 4     (x -> 0)

这个极限计入全带误差。

结果（fs = 48 kHz，max|ΔdB|）：20 / 117.6 / 691.4 / 4065 / 23900 Hz 分别为
1.1e-7 / 3.7e-6 / 1.3e-4 / 3.8e-3 / 4.8e-2 dB。作为对照，三点匹配 + 等波纹
极点为 1.5e-7 / 5.3e-6 / 2.2e-4 / 6.2e-3 / 2.1e-1 dB。

运行：``python qwqdsp/labs/match_filter/constraint_fit/highpass_1st_interp.py``
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

CHEB_NODES = 0.5 * (1.0 - np.cos((2.0 * np.arange(4) + 1.0) * PI / 8.0))
"""4 个 Chebyshev 根节点（在 x 上），作为 BRASIL 的初值。"""


# --------------------------------------------------------------------- 插值
def _match_from_nodes(fc, xnodes, fs=FS):
    """在给定 4 个节点上精确插值；失败返回 ``(None, 原因)``。"""
    x = np.sort(np.asarray(xnodes, dtype=float))
    if x.size != 4 or np.min(np.diff(x)) <= 1e-12:
        return None, "nodes"
    wc = 2.0 * PI * fc / fs
    w = 2.0 * np.arcsin(np.sqrt(x))
    r = w / wc
    target = r * r / (1.0 + r * r)
    M = np.stack([x, x * x, -target * x, -target * x * x], axis=1)
    scale = np.linalg.norm(M, axis=0)
    try:
        n1, n2, d1, d2 = np.linalg.solve(M / scale, target) / scale
    except np.linalg.LinAlgError:
        return None, "singular"

    # 由 x 基分母反推 (p, q)
    r2 = d1 + 1.0 + d2
    if not r2 > 0.0:
        return None, "r2 <= 0"
    rr = np.sqrt(r2)
    a1 = d1 + 2.0
    A = a1 - 1.0 - rr * rr
    B = 8.0 * (1.0 + rr)
    disc = B * B + 64.0 * A
    if disc < 0.0:
        return None, "no pole"
    roots = ([(-B + s * np.sqrt(disc)) / (2.0 * A) for s in (1.0, -1.0)]
             if abs(A) > 1e-15 else [16.0 / B])
    p_bilinear = 2.0 / (1.0 / np.tan(wc / 2.0) + 1.0)
    stable = [p for p in roots if p > 0.0 and p * (1.0 + rr) < 4.0]
    if not stable:
        return None, "unstable"
    p = min(stable, key=lambda z: abs(z - p_bilinear))
    q = rr * p

    a0 = p * p
    beta = a0 * np.array([0.0, n1 + n2, -n2 / 4.0])
    if not PM.feasible(beta):
        return None, "infeasible"
    a = np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])
    return PM.Match(PM.solve_b(beta), a, float(p), float(q), beta, PM.alpha_of(p, q)), "ok"


def error_db(match, x, fc, fs=FS):
    """在 ``x = sin^2(w/2)`` 上的 ``ΔdB``；x = 0 取极限。"""
    x = np.asarray(x, dtype=float)
    wc = 2.0 * PI * fc / fs
    ph = np.stack([1.0 - x, x, 4.0 * x * (1.0 - x)], axis=-1)
    hd2 = (ph @ match.beta) / (ph @ match.alpha)
    w = 2.0 * np.arcsin(np.sqrt(np.clip(x, 0.0, 1.0)))
    r = w / wc
    ha2 = r * r / (1.0 + r * r)
    out = np.zeros_like(x)
    pos = x > 0.0
    out[pos] = 10.0 * np.log10(np.maximum(hd2[pos], 1e-300)
                               / np.maximum(ha2[pos], 1e-300))
    if np.any(~pos):
        n1 = (match.beta[1] + 4.0 * match.beta[2]) / (match.p * match.p)
        out[~pos] = 10.0 * np.log10(max(n1 * wc * wc / 4.0, 1e-300))
    return out


# --------------------------------------------------------------------- BRASIL
def design(fc, fs=FS, max_iter=800, samples=300, step_max=0.15, kappa=0.15):
    """BRASIL 区间等化（5 段）；返回 ``(match, nodes, max|ΔdB|)``。"""
    nodes = CHEB_NODES.copy()
    best = (np.inf, None, None)
    for _ in range(max_iter):
        match, why = _match_from_nodes(fc, nodes, fs)
        if match is None:
            raise RuntimeError(f"fc = {fc:g} Hz：{why}")
        edges = np.concatenate(([0.0], nodes, [1.0]))
        errs = np.empty(5)
        for i in range(5):
            g = np.linspace(edges[i], edges[i + 1], samples + 2)
            if i > 0:
                g = g[1:]
            if i < 4:
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
        if errs.min() > 0.0 and top / errs.min() - 1.0 < 1e-6:
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
    print(f"fs = {FS:g} Hz；4 点有理插值 + BRASIL 区间等化（高通）")
    print(f"{'fc':>8}{'max|dB|':>12}   节点（x = sin^2(w/2)）")
    rows = []
    for fc in FC_LIST:
        match, nodes, err = design(fc)
        rows.append((fc, match))
        print(f"{fc:8.1f}{err:12.3e}   " + " ".join(f"{v:.5f}" for v in nodes))
    print()

    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    _figure(rows, out)
    print(f"图：{out / 'highpass_1st_interp.png'}")


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
        r = w / (2.0 * PI * fc / FS)
        ha2 = r * r / (1.0 + r * r)
        ax_r.semilogx(f, 10.0 * np.log10(np.maximum(ha2, 1e-300)),
                      color=color, lw=3.0, ls="--", alpha=0.45)
        ax_r.semilogx(f, 10.0 * np.log10(np.maximum(hd2, 1e-300)),
                      color=color, lw=1.3, label=f"{fc:g} Hz")
        xs = np.linspace(0.0, 1.0, 200001)
        err = error_db(match, xs, fc)
        ww = 2.0 * np.arcsin(np.sqrt(xs))
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
    ax_e.set_ylim(1e-8, 1e0)
    ax_e.set_xlabel("频率 / Hz")
    ax_e.set_ylabel("|ΔdB|")
    ax_e.grid(True, which="both", alpha=0.25)

    ax_r.set_title("一阶高通: P(0)=0+BRASIL 4点有理插值")
    fig.subplots_adjust(left=0.085, right=0.98, top=0.92, bottom=0.09, hspace=0.28)
    fig.savefig(out / "highpass_1st_interp.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
