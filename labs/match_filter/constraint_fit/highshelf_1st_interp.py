"""一阶 high shelf：结构约束 + 自由节点插值。

模拟功率（``g = 10^(db/20)``，``r = w/wc``）

    Pa(w) = (1 + g r^2) / (1 + r^2 / g)

只保留两条用途必需的硬约束：

    P(0) = 1,   P(pi) = Pa(pi)   （直流与高频增益）

fc 不作为硬约束，其余 3 个自由度由自由节点插值。结果（fs = 48 kHz，
max|ΔdB|）：0–0.052，全部档均可解。

运行：``python qwqdsp/labs/match_filter/constraint_fit/highshelf_1st_interp.py``
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
DB_LIST = (12.0, 6.0, -6.0, -12.0)


def analog_power(w, wc, db):
    g = 10.0 ** (db / 20.0)
    r = np.asarray(w, dtype=float) / wc
    return (1.0 + g * r * r) / (1.0 + r * r / g)


def _match(fc, db, nodes, fs=FS):
    wc = 2.0 * PI * fc / fs
    t = float(analog_power(np.array([PI]), wc, db)[0])
    rows = [[1.0, 0.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, -t, -t]]
    rhs = [1.0, t]
    for x in nodes:
        f = float(analog_power(np.array([2.0 * np.arcsin(np.sqrt(x))]), wc, db)[0])
        rows.append([1.0, x, x * x, -f * x, -f * x * x]); rhs.append(f)

    M = np.array(rows)
    scale = np.linalg.norm(M, axis=0)
    try:
        n0, n1, n2, d1, d2 = np.linalg.solve(M / scale, np.array(rhs)) / scale
    except np.linalg.LinAlgError:
        return None, "singular"

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
    beta = a0 * np.array([n0, n0 + n1 + n2, -n2 / 4.0])
    if not PM.feasible(beta):
        return None, "infeasible"
    a = np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])
    return PM.Match(PM.solve_b(beta), a, float(p), float(q), beta, PM.alpha_of(p, q)), "ok"


def error_db(match, x, fc, db, fs=FS):
    x = np.asarray(x, dtype=float)
    wc = 2.0 * PI * fc / fs
    ph = np.stack([1.0 - x, x, 4.0 * x * (1.0 - x)], axis=-1)
    hd2 = (ph @ match.beta) / (ph @ match.alpha)
    ha2 = analog_power(2.0 * np.arcsin(np.sqrt(np.clip(x, 0.0, 1.0))), wc, db)
    out = np.zeros_like(x)
    pos = x > 0.0
    out[pos] = 10.0 * np.log10(np.maximum(hd2[pos], 1e-300)
                               / np.maximum(ha2[pos], 1e-300))
    return out


def design(fc, db, fs=FS):
    wgrid = np.linspace(0.0, PI, 4001)
    xgrid = np.sin(wgrid / 2.0) ** 2

    def to_nodes(y):
        x = 1.0 / (1.0 + np.exp(-np.clip(np.asarray(y, float), -30.0, 30.0)))
        return np.sort(np.clip(x, 1e-9, 1.0 - 1e-9))

    def objective(y):
        nodes = to_nodes(y)
        if np.min(np.diff(nodes)) < 1e-9:
            return 1e3
        match, _ = _match(fc, db, nodes, fs)
        if match is None:
            return 1e3
        return float(np.max(np.abs(error_db(match, xgrid, fc, db, fs))))

    sets = ((0.05, 0.4, 0.9), (0.1, 0.5, 0.95), (0.01, 0.2, 0.7),
            (0.3, 0.7, 0.99), (0.5, 0.85, 0.99))
    starts = []
    for t in sets:
        t = tuple(float(np.clip(v, 1e-4, 1.0 - 1e-4)) for v in t)
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
    match, why = _match(fc, db, best[1], fs)
    if match is None:
        return None, None, np.inf
    xf = np.sin(np.linspace(0.0, PI, 40001) / 2.0) ** 2
    return match, best[1], float(np.max(np.abs(error_db(match, xf, fc, db, fs))))


def main():
    print(f"fs = {FS:g} Hz；一阶 high shelf：结构约束 + 自由节点")
    print(f"{'':>10}" + "".join(f"{fc:>12.1f}" for fc in FC_LIST))
    rows = {}
    for db in DB_LIST:
        cells = []
        for fc in FC_LIST:
            match, _, err = design(fc, db)
            rows[(db, fc)] = match
            cells.append(f"{err:12.4f}")
        print(f"{f'{db:+.0f} dB':>10}" + "".join(cells))
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    _figure(rows, out)
    print(f"图：{out / 'highshelf_1st_interp.png'}")


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
            match = rows[(db, fc)]
            if match is None:
                continue
            color = f"C{j}"
            ax_r.semilogx(f, 10.0 * np.log10(np.maximum(analog_power(w, 2.0 * PI * fc / FS, db), 1e-300)),
                          color=color, lw=3.0, ls="--", alpha=0.45)
            ax_r.semilogx(f, 10.0 * np.log10(np.maximum(match.power(w), 1e-300)),
                          color=color, lw=1.2, label=f"{fc:g} Hz")
            err = error_db(match, xs, fc, db)
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
    fig.suptitle("一阶 high shelf: P(0)=1, P(π)=Pa(π)+BRASIL 3点有理插值", y=0.995)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.96, bottom=0.05)
    fig.savefig(out / "highshelf_1st_interp.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
