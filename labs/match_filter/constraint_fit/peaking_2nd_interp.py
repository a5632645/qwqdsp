"""二阶 peaking（bell）：结构约束 + 自由节点插值。

模拟功率

    Pa(w) = ((1-r^2)^2 + (r A/Q)^2) / ((1-r^2)^2 + (r/(A Q))^2),
    A = 10^(db/40),  r = w/wc.

只保留两条用途必需的硬约束：

    P(0) = 1,   P(x_fc) = A^4   （中心增益）

Nyquist 不约束（远离钟形），其余 3 个自由度由自由节点插值。求解时把
P(0)=1 代入（``n0 = 1``），对剩余 4 个系数解 4×4，并按节点几何均值缩放
x，避免低 fc 时条件数退化。

结果（fs = 48 kHz，max|ΔdB|）：全档可解，0–0.18；4 kHz 从仓库的 0.30
降到 0.13。

运行：``python qwqdsp/labs/match_filter/constraint_fit/peaking_2nd_interp.py``
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "best_pole_fit"))
import power_match as PM  # noqa: E402

PI = PM.PI
FS = PM.FS
FC_LIST = PM.FC_LIST
Q_LIST = (0.5, 1.0 / np.sqrt(2.0), 1.0, 2.0)
DB_LIST = (12.0, 6.0, -6.0, -12.0)


def analog_power(w, wc, Q, db):
    r = np.asarray(w, dtype=float) / wc
    A = 10.0 ** (db / 40.0)
    return (((1.0 - r * r) ** 2 + (r * A / Q) ** 2)
            / ((1.0 - r * r) ** 2 + (r / (A * Q)) ** 2))


def _match(fc, Q, db, nodes, fs=FS):
    wc = 2.0 * PI * fc / fs
    xc = np.sin(wc / 2.0) ** 2
    t = 10.0 ** (db / 10.0)
    xs = np.concatenate(([xc], np.sort(np.asarray(nodes, float))))
    vals = np.array([t] + [
        float(analog_power(np.array([2.0 * np.arcsin(np.sqrt(x))]), wc, Q, db)[0])
        for x in np.sort(np.asarray(nodes, float))])
    ref = float(np.exp(np.mean(np.log(xs))))
    u = xs / ref
    M = np.stack([u, u * u, -vals * u, -vals * u * u], axis=1)
    scale = np.linalg.norm(M, axis=0)
    try:
        n1u, n2u, d1u, d2u = np.linalg.solve(M / scale, vals - 1.0) / scale
    except np.linalg.LinAlgError:
        return None, "singular"
    n1, n2 = n1u / ref, n2u / ref ** 2
    d1, d2 = d1u / ref, d2u / ref ** 2

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
    if not PM.feasible(beta):
        return None, "infeasible"
    a = np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])
    return PM.Match(PM.solve_b(beta), a, float(p), float(q), beta, PM.alpha_of(p, q)), "ok"


def error_db(match, x, fc, Q, db, fs=FS):
    x = np.asarray(x, dtype=float)
    wc = 2.0 * PI * fc / fs
    ph = np.stack([1.0 - x, x, 4.0 * x * (1.0 - x)], axis=-1)
    hd2 = (ph @ match.beta) / (ph @ match.alpha)
    ha2 = analog_power(2.0 * np.arcsin(np.sqrt(np.clip(x, 0.0, 1.0))), wc, Q, db)
    out = np.zeros_like(x)
    pos = x > 0.0
    out[pos] = 10.0 * np.log10(np.maximum(hd2[pos], 1e-300)
                               / np.maximum(ha2[pos], 1e-300))
    return out


def design(fc, Q, db, fs=FS):
    wc = 2.0 * PI * fc / fs
    xc = np.sin(wc / 2.0) ** 2
    wgrid = np.linspace(0.0, PI, 4001)
    xgrid = np.sin(wgrid / 2.0) ** 2

    def to_nodes(y):
        x = 1.0 / (1.0 + np.exp(-np.clip(np.asarray(y, float), -30.0, 30.0)))
        return np.sort(np.clip(x, 1e-9, 1.0 - 1e-9))

    def objective(y):
        nodes = to_nodes(y)
        if np.min(np.diff(nodes)) < 1e-9 or np.any(np.abs(nodes - xc) < 1e-12):
            return 1e3
        match, _ = _match(fc, Q, db, nodes, fs)
        if match is None:
            return 1e3
        return float(np.max(np.abs(error_db(match, xgrid, fc, Q, db, fs))))

    near = max(2.0 * xc, 1e-8)
    sets = ((0.05, 0.4, 0.9), (0.1, 0.5, 0.95), (0.01, 0.2, 0.7),
            (0.3, 0.7, 0.99), (0.5, 0.85, 0.99),
            (near, 0.1, 0.6), (near, 0.5, 0.95),
            (max(1.5 * xc, 1e-8), 0.01, 0.3))
    starts = []
    for t in sets:
        t = tuple(float(np.clip(v, 1e-12, 1.0 - 1e-9)) for v in t)
        if np.any(np.abs(np.array(t) - xc) < 1e-9):
            continue
        starts.append([np.log(v / (1.0 - v)) for v in t])

    best = (np.inf, None)
    for y0 in starts:
        res = minimize(lambda y: objective(y), np.array(y0, float), method="Nelder-Mead",
                       options={"xatol": 1e-9, "fatol": 1e-11,
                                "maxiter": 600, "maxfev": 600})
        if res.fun < best[0]:
            best = (float(res.fun), to_nodes(res.x))
    if best[1] is None or best[0] >= 1e3:
        return None, None, np.inf
    match, why = _match(fc, Q, db, best[1], fs)
    if match is None:
        return None, None, np.inf
    xf = np.sin(np.linspace(0.0, PI, 40001) / 2.0) ** 2
    return match, best[1], float(np.max(np.abs(error_db(match, xf, fc, Q, db, fs))))


def main():
    print(f"fs = {FS:g} Hz；二阶 peaking：结构约束 + 自由节点")
    print(f"{'':>14}" + "".join(f"{fc:>12.1f}" for fc in FC_LIST))
    rows = {}
    for Q in Q_LIST:
        for db in DB_LIST:
            cells = []
            for fc in FC_LIST:
                match, _, err = design(fc, Q, db)
                rows[(Q, db, fc)] = match
                cells.append(f"{err:12.4f}")
            print(f"{f'Q={Q:g} {db:+.0f}dB':>14}" + "".join(cells))
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    _figure(rows, out)
    print(f"图：{out / 'peaking_2nd_interp.png'}")


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
    for r, db in enumerate(DB_LIST):
        ax_r, ax_e = axes[r]
        for j, fc in enumerate(FC_LIST):
            match = rows[(1.0, db, fc)]
            if match is None:
                continue
            color = f"C{j}"
            ax_r.semilogx(f, 10.0 * np.log10(np.maximum(analog_power(w, 2.0 * PI * fc / FS, 1.0, db), 1e-300)),
                          color=color, lw=3.0, ls="--", alpha=0.45)
            ax_r.semilogx(f, 10.0 * np.log10(np.maximum(match.power(w), 1e-300)),
                          color=color, lw=1.2, label=f"{fc:g} Hz")
            err = error_db(match, xs, fc, 1.0, db)
            ax_e.plot(ww / PI * FS / 2.0, np.maximum(np.abs(err), 1e-12),
                      color=color, lw=1.1)
        ax_r.set_xscale("log")
        ax_r.set_xlim(10.0, FS / 2.0)
        ax_r.set_ylim(-90.0, 25.0)
        ax_r.set_ylabel(f"{db:+.0f} dB\n幅度 / dB")
        ax_r.grid(True, which="both", alpha=0.25)
        ax_e.set_xlim(0.0, FS / 2.0)
        ax_e.set_yscale("log")
        ax_e.set_ylim(1e-8, 1e1)
        ax_e.set_ylabel("|ΔdB|")
        ax_e.grid(True, which="both", alpha=0.25)
    axes[0][0].legend(fontsize=7, title="fc", ncol=2)
    axes[-1][0].set_xlabel("频率 / Hz")
    axes[-1][1].set_xlabel("频率 / Hz")
    fig.suptitle("二阶 peaking: P(0)=1, P(fc)=A⁴+BRASIL 3点有理插值（Q=1）", y=0.995)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.96, bottom=0.05)
    fig.savefig(out / "peaking_2nd_interp.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
