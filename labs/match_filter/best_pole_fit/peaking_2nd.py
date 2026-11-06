"""二阶模拟 peaking（bell）：三点匹配 + 最佳极点（boost × Q）。

模拟原型

    H(S) = (S^2 + (A/Q) S + 1) / (S^2 + S/(A Q) + 1),   A = 10^(db/40),

即极点 ``Qp = A Q``、零点 ``Qz = Q/A``，峰值功率 ``|H(j)|^2 = A^4 = 10^(db/10)``。

数字节按功率函数在 **0、fc、fs/2** 三点精确匹配，得到分子；分母两个自由值按
全带等波纹（最小最大 ``|ΔdB|``）优化。II 基线取 Vicanek 的 Peaking：II 极点
``Qp = A Q`` 加**两点（DC + fc）**分子——它的 Nyquist 不匹配，cut 档尤其明显，
这正是极点优化要补的缺口。极点优化把 II 解、3 点解与密集网格一起比较，保证
不劣于 II。

参数：4 个 boost（+12/+6/−6/−12 dB）× 4 个 Q（0.5/0.707/1/2）× 5 个 fc。

运行：``python qwqdsp/labs/match_filter/best_pole_fit/peaking_2nd.py``
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

import power_match as PM

PI = PM.PI
"""圆周率。"""

FS = PM.FS
"""实验采样率。"""

FC_LIST = PM.FC_LIST
"""与 match_filter 其余实验一致的 5 个截止频率。"""

BOOST_LIST = (12.0, 6.0, -6.0, -12.0)
"""四个峰值（dB）。"""

Q_LIST = (0.5, 1.0 / np.sqrt(2.0), 1.0, 2.0)
"""四个 Q。"""


# --------------------------------------------------------------------- 功率函数
def analog_power(w, wc, Q, db):
    """二阶模拟 peaking 的功率：DC/Nyquist 为 1，fc 处为 10^(db/10)。"""
    r = np.asarray(w, dtype=float) / wc
    A = 10.0 ** (db / 40.0)
    num = (1.0 - r * r) ** 2 + (r * A / Q) ** 2
    den = (1.0 - r * r) ** 2 + (r / (A * Q)) ** 2
    return num / den


def build_match(fc, fs, Q, db, p, q):
    """三点（0/fc/fs2）功率匹配的分子 + 给定分母。"""
    if not (p > 0.0 and q > 0.0 and p + q < 4.0):
        return None
    wc = 2.0 * PI * fc / fs
    alpha = PM.alpha_of(p, q)
    phi = PM.phi(np.array([0.0, wc, PI]))
    target = np.array([1.0, 10.0 ** (db / 10.0),
                       analog_power(PI, wc, Q, db)])
    beta = np.linalg.solve(phi, target * (phi @ alpha))
    if not PM.feasible(beta):
        return None
    a = np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])
    return PM.Match(PM.solve_b(beta), a, float(p), float(q), beta, alpha)


def build_vicanek(fc, fs, Q, db, p, q):
    """Vicanek Peaking：II 极点 + DC/fc 两点分子。"""
    if not (p > 0.0 and q > 0.0 and p + q < 4.0):
        return None
    wc = 2.0 * PI * fc / fs
    alpha = PM.alpha_of(p, q)
    phi = PM.phi(np.array([wc]))[0]
    gain2 = 10.0 ** (db / 10.0)
    r1 = gain2 * float(phi @ alpha)
    r2 = gain2 * (-alpha[0] + alpha[1] + 4.0 * (phi[0] - phi[1]) * alpha[2])
    beta0 = alpha[0]
    beta2 = (r1 - r2 * phi[1] - beta0) / (4.0 * phi[1] * phi[1])
    beta1 = r2 + beta0 + 4.0 * (phi[1] - phi[0]) * beta2
    beta = np.array([beta0, beta1, beta2])
    if not PM.feasible(beta):
        return None
    a = np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])
    return PM.Match(PM.solve_b(beta), a, float(p), float(q), beta, alpha)


def error_db(match, fc, fs, Q, db, n=20001):
    s = np.linspace(0.0, 1.0, n)
    ph = np.stack([1.0 - s, s, 4.0 * s * (1.0 - s)], axis=-1)
    hd2 = (ph @ match.beta) / (ph @ match.alpha)
    w = 2.0 * np.arcsin(np.sqrt(s))
    ha2 = analog_power(w, 2.0 * PI * fc / fs, Q, db)
    return 10.0 * np.log10(np.maximum(hd2, 1e-300) / np.maximum(ha2, 1e-300)), s, w


def max_error_db(fc, fs, Q, db, p, q, n=20001):
    match = build_match(fc, fs, Q, db, p, q)
    if match is None:
        return np.inf
    err, _, _ = error_db(match, fc, fs, Q, db, n)
    return float(np.max(np.abs(err)))


# --------------------------------------------------------------------- 极点
def vicanek_pq(fc, fs, Q, db):
    """II 极点对应的 (p,q)：原型极点 Qp = A Q。"""
    A = 10.0 ** (db / 40.0)
    zeta = 1.0 / (2.0 * A * Q)
    wc = 2.0 * PI * fc / fs
    e = np.exp(-zeta * wc)
    mu = np.sqrt(max(0.0, 1.0 - zeta * zeta))
    a1 = -2.0 * e * np.cos(mu * wc)
    a2 = e * e
    return 1.0 + a1 + a2, 1.0 - a1 + a2


def best_design(fc, fs, Q, db):
    from scipy.optimize import minimize

    p_ii, q_ii = vicanek_pq(fc, fs, Q, db)
    ii_match = build_vicanek(fc, fs, Q, db, p_ii, q_ii)
    ii_err = (float(np.max(np.abs(error_db(ii_match, fc, fs, Q, db, 1000001)[0])))
              if ii_match is not None else np.inf)

    candidates = []
    seed = build_match(fc, fs, Q, db, p_ii, q_ii)
    if seed is not None:
        candidates.append((float(np.max(np.abs(
            error_db(seed, fc, fs, Q, db, 4001)[0]))), p_ii, q_ii))
    for p in np.geomspace(1e-4, 3.9, 30):
        for q in np.linspace(1e-4, 3.9, 30):
            if p + q >= 4.0:
                continue
            value = max_error_db(fc, fs, Q, db, p, q, 2001)
            if np.isfinite(value):
                candidates.append((value, p, q))
    candidates.sort()

    best = (np.inf, None)
    for _, p, q in candidates[:3]:
        res = minimize(
            lambda x: max_error_db(fc, fs, Q, db, x[0], x[1], 8001),
            np.array([p, q]), method="Nelder-Mead",
            options={"xatol": 1e-12, "fatol": 1e-12,
                     "maxiter": 2000, "maxfev": 1500})
        if res.fun < best[0]:
            best = (float(res.fun), res.x)
    if best[1] is None:
        return ii_match

    match = build_match(fc, fs, Q, db, *best[1])
    ours = float(np.max(np.abs(error_db(match, fc, fs, Q, db, 1000001)[0])))
    if ii_match is not None and ii_err <= ours:
        return ii_match
    return match


def design_peaking(fc, fs, Q, db):
    match = best_design(fc, fs, Q, db)
    if match is None:
        raise RuntimeError("没有可行的极点")
    return match


# --------------------------------------------------------------------- 实验
def _format_pole(z):
    if abs(z.imag) < 1e-9:
        return f"{z.real:.6f}"
    return f"{z.real:.6f}{z.imag:+.6f}j"


def _report(fs):
    print(f"fs = {fs:g} Hz，全带 max|ΔdB|（三点；II = Vicanek 两点）")
    print(f"{'boost':>6}{'Q':>8}{'fc':>8}{'p':>13}{'q':>12}"
          f"{'z1':>16}{'z2':>16}{'ours':>9}{'II':>9}")
    rows = []
    for db in BOOST_LIST:
        for Q in Q_LIST:
            for fc in FC_LIST:
                match = design_peaking(fc, fs, Q, db)
                err, _, _ = error_db(match, fc, fs, Q, db, 1000001)
                p_ii, q_ii = vicanek_pq(fc, fs, Q, db)
                ii_match = build_vicanek(fc, fs, Q, db, p_ii, q_ii)
                if ii_match is not None:
                    ii_err = float(np.max(np.abs(
                        error_db(ii_match, fc, fs, Q, db, 1000001)[0])))
                else:
                    ii_err = np.inf
                poles = match.poles()
                rows.append({"db": db, "Q": Q, "fc": fc, "match": match,
                             "poles": poles, "ours": float(np.max(np.abs(err))),
                             "ii": ii_err, "ii_match": ii_match})
                print(f"{db:6.0f}{Q:8.4f}{fc:8.1f}{match.p:13.6g}{match.q:12.6g}"
                      f"{_format_pole(poles[0]):>16}{_format_pole(poles[1]):>16}"
                      f"{float(np.max(np.abs(err))):9.4f}{ii_err:9.4f}")
        print()
    return rows


def _figure_panels(rows, fs, db, out_dir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    out_dir.mkdir(exist_ok=True)

    fig = plt.figure(figsize=(12.5, 8.8))
    outer = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.22)
    f = np.geomspace(10.0, fs / 2.0, 4001)
    w = 2.0 * PI * f / fs
    for index, Q in enumerate(Q_LIST):
        cell = outer[index // 2, index % 2].subgridspec(
            2, 1, height_ratios=(2.1, 1.1), hspace=0.25)
        ax_r = fig.add_subplot(cell[0])
        ax_e = fig.add_subplot(cell[1])
        group = [row for row in rows if row["db"] == db and row["Q"] == Q]
        for j, row in enumerate(group):
            color = f"C{j}"
            ph = PM.phi(w)
            hd2 = (ph @ row["match"].beta) / (ph @ row["match"].alpha)
            ha2 = analog_power(w, 2.0 * PI * row["fc"] / fs, Q, db)
            ax_r.semilogx(f, 10.0 * np.log10(np.maximum(ha2, 1e-300)),
                          color=color, lw=3.4, ls="--", alpha=0.30)
            ax_r.semilogx(f, 10.0 * np.log10(np.maximum(hd2, 1e-300)),
                          color=color, lw=1.4, label=f"{row['fc']:g} Hz")
            err, _, ww = error_db(row["match"], row["fc"], fs, Q, db, 200001)
            ax_e.semilogy(ww / PI, np.maximum(np.abs(err), 1e-12),
                          color=color, lw=1.2)
        ax_r.set_title(f"Q = {Q:g}")
        ax_r.set_ylabel("幅度 / dB")
        ax_r.set_xlim(10.0, fs / 2.0)
        ax_r.set_ylim(-24.0, 18.0)
        ax_r.set_yticks([-24.0, -18.0, -12.0, -6.0, 0.0, 6.0, 12.0, 18.0])
        ax_r.grid(True, which="both", alpha=0.25)
        ax_r.legend(loc="lower left", fontsize=7)
        ax_r.set_xlabel("频率 / Hz")
        ax_e.set_ylim(1e-8, 1e2)
        ax_e.set_ylabel("|ΔdB|")
        ax_e.set_xlabel("f / (fs/2)")
        ax_e.grid(True, which="both", alpha=0.25)

    fig.suptitle(f"二阶 peaking（boost = {db:+g} dB）：四个 Q 子类"
                 "（实线：最佳二阶节；虚线：模拟原型）")
    path = out_dir / f"peaking_2nd_response_boost{db:+g}.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _pole_metrics(poles):
    z = poles[np.argmax(np.abs(poles.imag))]
    if abs(z.imag) > 1e-9:
        return abs(z), abs(float(np.angle(z)))
    z = poles[np.argmax(np.abs(poles))]
    return abs(z), (PI if z.real < 0.0 else 0.0)


def _figure_poles(rows, fs, out_dir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    out_dir.mkdir(exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.8))
    colors = ("C0", "C1", "C2", "C3")
    for ax, db in zip(axes.ravel(), BOOST_LIST):
        ax_r = ax
        ax_a = ax.twinx()
        handles = []
        for Q, color in zip(Q_LIST, colors):
            group = sorted((row for row in rows if row["db"] == db and row["Q"] == Q),
                           key=lambda row: row["fc"])
            fcs = [row["fc"] for row in group]
            radius, angle = [], []
            for row in group:
                r, theta = _pole_metrics(row["poles"])
                radius.append(r)
                angle.append(np.degrees(theta))
            line_r, = ax_r.semilogx(fcs, radius, "o-", color=color, lw=1.3,
                                    ms=4.5, label=f"Q={Q:g} 半径")
            line_a, = ax_a.semilogx(fcs, angle, "s--", color=color, lw=1.0,
                                    ms=4, mfc="none", label=f"Q={Q:g} 角度")
            handles += [line_r, line_a]
        ax_r.set_ylim(0.0, 1.05)
        ax_r.set_ylabel("极点半径 |z|")
        ax_a.set_ylim(0.0, 180.0)
        ax_a.set_yticks([0.0, 45.0, 90.0, 135.0, 180.0])
        ax_a.set_ylabel("极点角度 / deg")
        ax_r.set_xlabel("fc / Hz")
        ax_r.set_title(f"boost = {db:+g} dB")
        ax_r.grid(True, which="both", alpha=0.25)
        ax_r.legend(handles=handles, loc="center left", fontsize=6.5, ncol=2)

    fig.suptitle("最佳极点：半径（实线/实心）与角度（虚线/空心），颜色区分 Q")
    fig.tight_layout()
    path = out_dir / "peaking_2nd_poles.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main():
    rows = _report(FS)
    out_dir = Path(__file__).parent / "output"
    paths = [_figure_panels(rows, FS, db, out_dir) for db in BOOST_LIST]
    paths.append(_figure_poles(rows, FS, out_dir))
    for path in paths:
        print(f"图：{path}")


if __name__ == "__main__":
    main()
