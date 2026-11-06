"""二阶带通：4 点有理插值 + BRASIL 区间等化。

二阶模拟带通（归一化峰值 1）

    Pa(w) = (r^2/Q^2) / ((1 - r^2)^2 + (r/Q)^2),    r = w / wc

在 ``x = sin^2(w/2)`` 下仍为 2/2 有理函数。带通在直流处为零，所以强制

    P(0) = 0  =>  n0 = 0,

分子只剩一次项起：

    P(x) = x (n1 + n2 x) / (1 + d1 x + d2 x^2)

即只约束了一个零点（z = 1 处），其余 4 个参数由 4 个节点插值确定；最佳逼近
在 5 个点等幅交替，BRASIL 切成 5 段做长度等化。

直流处数字与模拟功率都为零，dB 误差取极限

    lim P/P a = n1 * wc^2 * Q^2 / 4     (x -> 0)

Nyquist 不约束：目标 Pa(pi) 是有限值，不是零。

对 Q = 0.5 / 0.707 / 1 / 2 各扫 5 个 fc。结果（fs = 48 kHz，max|ΔdB|）
约 0.32–0.48 dB（fc <= 4 kHz），23.9 kHz 为 0.003–0.019 dB。

运行：``python qwqdsp/labs/match_filter/constraint_fit/bandpass_2nd_interp.py``
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

Q_LIST = (0.5, 1.0 / np.sqrt(2.0), 1.0, 2.0)
"""四个 Q：临界阻尼、Butterworth、1、以及有谐振峰的 2。"""


def analog_power(w, wc, Q):
    r = np.asarray(w, dtype=float) / wc
    return (r * r / (Q * Q)) / ((1.0 - r * r) ** 2 + (r / Q) ** 2)


def _match_from_nodes(fc, Q, xnodes, fs=FS):
    """在 4 个节点上精确插值；失败返回 ``(None, 原因)``。"""
    x = np.sort(np.asarray(xnodes, dtype=float))
    if x.size != 4 or np.min(np.diff(x)) <= 1e-12:
        return None, "nodes"
    wc = 2.0 * PI * fc / fs
    w = 2.0 * np.arcsin(np.sqrt(np.clip(x, 0.0, 1.0)))
    target = analog_power(w, wc, Q)
    M = np.stack([x, x * x, -target * x, -target * x * x], axis=1)
    scale = np.linalg.norm(M, axis=0)
    try:
        n1, n2, d1, d2 = np.linalg.solve(M / scale, target) / scale
    except np.linalg.LinAlgError:
        return None, "singular"
    if not n1 > 0.0:
        return None, "n1 <= 0"

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


def error_db(match, x, fc, Q, fs=FS):
    """在 ``x = sin^2(w/2)`` 上的 ``ΔdB``；x = 0 取极限。"""
    x = np.asarray(x, dtype=float)
    wc = 2.0 * PI * fc / fs
    a0 = match.p * match.p
    n1 = (match.beta[1] + 4.0 * match.beta[2]) / a0
    n2 = -4.0 * match.beta[2] / a0
    ph = np.stack([1.0 - x, x, 4.0 * x * (1.0 - x)], axis=-1)
    hd2 = a0 * (n1 * x + n2 * x * x) / (ph @ match.alpha)
    w = 2.0 * np.arcsin(np.sqrt(np.clip(x, 0.0, 1.0)))
    ha2 = analog_power(w, wc, Q)
    out = np.zeros_like(x)
    pos = x > 0.0
    out[pos] = 10.0 * np.log10(np.maximum(hd2[pos], 1e-300)
                               / np.maximum(ha2[pos], 1e-300))
    if np.any(~pos):
        out[~pos] = 10.0 * np.log10(max(n1 * wc * wc * Q * Q / 4.0, 1e-300))
    return out


def _inits(fc, fs=FS):
    xf = min(max(np.sin(PI * fc / fs) ** 2, 1e-9), 0.9995)
    return ([xf, 1e-5, 1e-3, 0.5], [xf * 1.2, 1e-5, 1e-2, 0.5],
            [1e-5, 1e-3, 0.05, 0.5], [1e-3, 0.01, 0.1, 0.5],
            [0.01, 0.1, 0.5, 0.9], [xf, 1e-4, 0.01, 0.3])


def design(fc, Q, fs=FS, max_iter=500, step_max=0.15, kappa=0.15):
    """多初值 + 回溯的 BRASIL；返回 ``(match, nodes, max|ΔdB|)``。"""
    wgrid = np.linspace(0.0, PI, 4001)
    xgrid = np.sin(wgrid / 2.0) ** 2
    best = (np.inf, None, None)
    for nodes0 in _inits(fc, fs):
        nodes = np.sort(np.asarray(nodes0, dtype=float))
        match, why = _match_from_nodes(fc, Q, nodes, fs)
        if match is None:
            continue
        local_best = np.inf
        stall = 0
        for _ in range(max_iter):
            edges = np.concatenate(([0.0], nodes, [1.0]))
            E = np.abs(error_db(match, xgrid, fc, Q, fs))
            idx = np.searchsorted(xgrid, edges)
            errs = np.empty(5)
            for i in range(5):
                lo, hi = idx[i], min(idx[i + 1] + 1, xgrid.size)
                errs[i] = E[lo:hi].max() if hi > lo else 0.0
            top = float(errs.max())
            if top < best[0]:
                best = (top, match, nodes.copy())
            if top < local_best - 1e-9:
                local_best = top
                stall = 0
            else:
                stall += 1
                if stall > 50:
                    break
            mean = errs.mean()
            dev = errs - mean
            dmax = np.abs(dev).max()
            if dmax <= 0.0 or (errs.min() > 0.0 and top / errs.min() - 1.0 < 1e-6):
                break
            base = min(step_max, kappa / (dmax / mean)) if mean > 0.0 else step_max
            for bt in range(7):
                lengths = (1.0 + base * 0.5 ** bt) ** (-dev / dmax) * np.diff(edges)
                lengths /= lengths.sum()
                cand = np.concatenate(([0.0], np.cumsum(lengths)))[1:-1]
                mc, _ = _match_from_nodes(fc, Q, cand, fs)
                if mc is not None:
                    nodes, match = cand, mc
                    break
            else:
                break

    _, match, nodes = best
    if match is None:
        raise RuntimeError(f"fc = {fc:g} Hz, Q = {Q:g}：没有可行解")
    xf_grid = np.sin(np.linspace(0.0, PI, 40001) / 2.0) ** 2
    err = float(np.max(np.abs(error_db(match, xf_grid, fc, Q, fs))))
    return match, nodes, err


def main():
    print(f"fs = {FS:g} Hz；二阶带通，4 点有理插值 + BRASIL")
    print(f"{'Q':>8}" + "".join(f"{fc:>12.1f}" for fc in FC_LIST))
    rows = {}
    for Q in Q_LIST:
        cells = []
        for fc in FC_LIST:
            match, nodes, err = design(fc, Q)
            rows[(Q, fc)] = match
            cells.append(f"{err:12.4f}")
        print(f"{Q:8.4f}" + "".join(cells))
    print()

    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    _figure(rows, out)
    print(f"图：{out / 'bandpass_2nd_interp.png'}")


def _figure(rows, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(4, 2, figsize=(12.5, 14.0),
                             gridspec_kw={"wspace": 0.22, "hspace": 0.35})
    f = np.geomspace(10.0, FS / 2.0, 3501)
    w = 2.0 * PI * f / FS
    xs = np.sin(np.linspace(0.0, PI, 20001) / 2.0) ** 2
    ww = 2.0 * np.arcsin(np.sqrt(xs))
    for r, Q in enumerate(Q_LIST):
        ax_r, ax_e = axes[r]
        for j, fc in enumerate(FC_LIST):
            match = rows[(Q, fc)]
            color = f"C{j}"
            hd2 = match.power(w)
            ha2 = analog_power(w, 2.0 * PI * fc / FS, Q)
            ax_r.semilogx(f, 10.0 * np.log10(np.maximum(ha2, 1e-300)),
                          color=color, lw=3.0, ls="--", alpha=0.45)
            ax_r.semilogx(f, 10.0 * np.log10(np.maximum(hd2, 1e-300)),
                          color=color, lw=1.2, label=f"{fc:g} Hz")
            err = error_db(match, xs, fc, Q)
            ax_e.plot(ww / PI * FS / 2.0, np.maximum(np.abs(err), 1e-12),
                      color=color, lw=1.1)
        ax_r.set_xscale("log")
        ax_r.set_xlim(10.0, FS / 2.0)
        ax_r.set_ylim(-110.0, 25.0)
        ax_r.set_ylabel(f"Q = {Q:g}\n幅度 / dB")
        ax_r.grid(True, which="both", alpha=0.25)
        ax_e.set_xlim(0.0, FS / 2.0)
        ax_e.set_yscale("log")
        ax_e.set_ylim(1e-8, 1e1)
        ax_e.set_ylabel("|ΔdB|")
        ax_e.grid(True, which="both", alpha=0.25)
    axes[0][0].legend(fontsize=7, title="fc", ncol=2)
    axes[-1][0].set_xlabel("频率 / Hz")
    axes[-1][1].set_xlabel("频率 / Hz")
    fig.suptitle("二阶带通: P(0)=0+BRASIL 4点有理插值", y=0.995)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.965, bottom=0.05)
    fig.savefig(out / "bandpass_2nd_interp.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
