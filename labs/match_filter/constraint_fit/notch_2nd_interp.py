"""二阶 notch：结构约束 + BRASIL 区间等化。

模拟 notch 的功率（``r = w/wc``）

    Pa(w) = (1 - r^2)^2 / ((1 - r^2)^2 + (r/Q)^2)

硬约束：

    P(0) = 1,   P(x_fc) = 0

功率在 ``x = sin^2(w/2)`` 下是 2/2 有理函数。缺口要求分子在 ``x_fc``
处为零；要保证功率非负，这个零点必须是二重的。于是分子被
``P(0)=1`` 完全定死：

    N(x) = (1 - x/xc)^2,    xc = sin^2(wc/2)

只剩分母 ``1 + d1 x + d2 x^2`` 的 2 个自由度。用 2 个自由节点插值，
BRASIL 调整节点位置，使 3 段的 max|ΔdB| 等幅；``x_fc`` 在段内部，
不切段。缺口处的 0/0 取解析极限（两边的二阶曲率之比）。

分子落在可行性边界 ``W^2 + B2 = 0`` 上，因此还原分子时不能要求
严格 ``> 0``，而是允许边界并在 ``rad`` 略负时截断。

结果（fs = 48 kHz，max|ΔdB|）：fc <= 4 kHz 为 0.0000–0.103，
23.9 kHz 为 1.467–1.695；全表 0–1.70。

运行：``python qwqdsp/labs/match_filter/constraint_fit/notch_2nd_interp.py``
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "best_pole_fit"))
import brasil as BRASIL  # noqa: E402
import power_match as PM  # noqa: E402

PI = PM.PI
FS = PM.FS
FC_LIST = PM.FC_LIST
Q_LIST = (0.5, 1.0 / np.sqrt(2.0), 1.0, 2.0)


def analog_power(w, wc, Q):
    r = np.asarray(w, dtype=float) / wc
    a = (1.0 - r * r) ** 2
    return a / (a + (r / Q) ** 2)


def _to_match(n1, n2, d1, d2, wc):
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
    beta = a0 * np.array([1.0, 1.0 + n1 + n2, -n2 / 4.0])
    # notch 的分子落在可行性边界 W^2 + B2 = 0 上，不能要求严格 > 0
    if not (beta[0] > 0.0 and beta[1] > 0.0):
        return None, "infeasible"
    W = (np.sqrt(beta[0]) + np.sqrt(beta[1])) / 2.0
    rad = W * W + beta[2]
    if rad < -1e-8 * max(W * W, abs(beta[2]), 1e-300):
        return None, "infeasible"
    b0 = (W + np.sqrt(max(rad, 0.0))) / 2.0
    b1 = (np.sqrt(beta[0]) - np.sqrt(beta[1])) / 2.0
    b2 = -beta[2] / (4.0 * b0)
    a = np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])
    return PM.Match(np.array([b0, b1, b2]), a, float(p), float(q), beta,
                    PM.alpha_of(p, q)), "ok"


def _match(fc, Q, nodes, fs=FS):
    wc = 2.0 * PI * fc / fs
    xc = np.sin(wc / 2.0) ** 2
    xs = np.sort(np.asarray(nodes, float))
    if (xs.size != 2 or np.min(np.diff(xs)) <= 1e-12
            or np.any(np.abs(xs - xc) < 1e-12)):
        return None, "nodes"
    rs = np.asarray(2.0 * np.arcsin(np.sqrt(xs)), float)
    vals = analog_power(rs, wc, Q)
    N = (xs - xc) ** 2 / xc ** 2
    M = np.stack([xs, xs * xs], axis=1)
    scale = np.linalg.norm(M, axis=0)
    try:
        d1, d2 = np.linalg.solve(M / scale, N / vals - 1.0) / scale
    except np.linalg.LinAlgError:
        return None, "singular"
    xx = np.linspace(0.0, 1.0, 1001)
    if np.min(1.0 + d1 * xx + d2 * xx * xx) <= 0.0:
        return None, "D <= 0"
    return _to_match(-2.0 / xc, 1.0 / xc ** 2, d1, d2, wc)


def error_db(match, x, fc, Q, fs=FS):
    x = np.asarray(x, dtype=float)
    wc = 2.0 * PI * fc / fs
    xc = np.sin(wc / 2.0) ** 2
    a_scale = match.alpha[0]
    d1 = (match.alpha[1] + 4.0 * match.alpha[2]) / a_scale - 1.0
    d2 = -4.0 * match.alpha[2] / a_scale
    n2 = -4.0 * match.beta[2] / match.beta[0]
    n1 = match.beta[1] / match.beta[0] - 1.0 - n2
    D = 1.0 + d1 * x + d2 * x * x
    N = 1.0 + n1 * x + n2 * x * x
    w = 2.0 * np.arcsin(np.sqrt(np.clip(x, 0.0, 1.0)))
    r = w / wc
    dx = x - xc
    tol = 1e-13 * max(xc, 1e-300)
    dwdx = np.empty_like(x)
    near = np.abs(dx) <= tol
    dwdx[near] = 1.0 / np.sqrt(xc * (1.0 - xc))
    dwdx[~near] = (w[~near] - wc) / dx[~near]
    a_over = ((dwdx * (1.0 + r)) / wc) ** 2
    ratio = (1.0 / xc ** 2) / a_over * (
        ((1.0 - r * r) ** 2 + (r / Q) ** 2) / D)
    return 10.0 * np.log10(np.maximum(ratio, 1e-300))


def _measure_grid(fc, fs=FS):
    xc = np.sin(PI * fc / fs) ** 2
    xs = np.linspace(0.0, 1.0, 4001)
    width = max(4.0 * xc * (1.0 - xc), 1e-12)
    loc = xc + np.linspace(-1.0, 1.0, 201) * min(width, 0.25)
    return np.unique(np.clip(np.concatenate((xs, loc, [xc])), 0.0, 1.0))


def design(fc, Q, fs=FS):
    xc = np.sin(PI * fc / fs) ** 2
    xgrid = _measure_grid(fc, fs)

    def measure(nodes):
        if np.any(np.abs(nodes - xc) < 1e-12):
            return None
        match, _ = _match(fc, Q, nodes, fs)
        if match is None:
            return None
        return BRASIL.segment_maxima(
            xgrid, np.abs(error_db(match, xgrid, fc, Q, fs)), nodes)

    near = max(2.0 * xc, 1e-8)
    sets = ((0.05, 0.5), (0.1, 0.9), (0.01, 0.7), (0.3, 0.99),
            (0.5, 0.95), (max(0.5 * xc, 1e-8), near),
            (max(0.1 * xc, 1e-9), 0.5), (0.01, max(10.0 * xc, 1e-8)),
            (0.99, 0.9999), (0.999, 0.99999))
    best = (np.inf, None)
    for t in sets:
        nodes0 = np.sort(np.clip(np.asarray(t, float), 1e-9, 1.0 - 1e-9))
        if (np.min(np.diff(nodes0)) <= 1e-12
                or np.any(np.abs(nodes0 - xc) < 1e-12)):
            continue
        nodes, _, err = BRASIL.equalize(nodes0, measure)
        if nodes is not None and err < best[0]:
            best = (err, nodes)
    if best[1] is None:
        return None, None, np.inf
    match, _ = _match(fc, Q, best[1], fs)
    if match is None:
        return None, None, np.inf
    return match, best[1], float(np.max(np.abs(
        error_db(match, xgrid, fc, Q, fs))))


def main():
    print(f"fs = {FS:g} Hz；二阶 notch：P(0)=1, P(fc)=0 + BRASIL")
    print(f"{'':>12}" + "".join(f"{fc:>12.1f}" for fc in FC_LIST))
    rows = {}
    for Q in Q_LIST:
        cells = []
        for fc in FC_LIST:
            match, _, err = design(fc, Q)
            rows[(Q, fc)] = match
            cells.append(f"{err:12.4f}")
        print(f"{f'Q={Q:g}':>12}" + "".join(cells))
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    _figure(rows, out)
    print(f"图：{out / 'notch_2nd_interp.png'}")


def _figure(rows, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(4, 2, figsize=(12.5, 14.0),
                             gridspec_kw={"wspace": 0.22, "hspace": 0.35})
    f = np.geomspace(10.0, FS / 2.0, 3001)
    w = 2.0 * PI * f / FS
    xs = np.sin(np.linspace(0.0, PI, 20001) / 2.0) ** 2
    ww = 2.0 * np.arcsin(np.sqrt(xs))
    for r, Q in enumerate(Q_LIST):
        ax_r, ax_e = axes[r]
        for j, fc in enumerate(FC_LIST):
            match = rows[(Q, fc)]
            if match is None:
                continue
            color = f"C{j}"
            ax_r.semilogx(
                f, 10.0 * np.log10(np.maximum(
                    analog_power(w, 2.0 * PI * fc / FS, Q), 1e-300)),
                color=color, lw=3.0, ls="--", alpha=0.45)
            ax_r.semilogx(
                f, 10.0 * np.log10(np.maximum(match.power(w), 1e-300)),
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
    fig.suptitle("二阶 notch: P(0)=1, P(fc)=0+BRASIL 2点有理插值", y=0.995)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.96, bottom=0.05)
    fig.savefig(out / "notch_2nd_interp.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
