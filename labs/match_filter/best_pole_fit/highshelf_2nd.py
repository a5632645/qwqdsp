"""二阶高 shelf（任意 Q）：二点（DC + fc）匹配 + 最佳极点。

模拟原型（``A = 10^(db/40)``、``S = j r``）：

    H(S) = A * (A S^2 + (sqrt(A)/Q) S + 1) / (S^2 + (sqrt(A)/Q) S + A),

极点 ``wn = wc * A^(1/2)``、零点 ``wn = wc / A^(1/2)``，两者 Q 相同；DC 为 1，
fc 处幅度 ``sqrt(g)``（功率 ``g = 10^(db/20)``），高频渐近 ``g``。

三点（0、fc、fs/2）匹配完全不用；种子取 II 极点与 ``(p, q)`` 网格。最终
设计精确匹配 0/fc 幅度与 fc 处斜率（Hermite；0 点斜率恒为 0），fs/2 不约束；
给定极点后分子是 3×3 线性解（闭式）。极点半径限制在 ``rho <= 0.999``，防止
贴单位圆的近对消解。Vicanek 的 ``HighShelf`` 作为 II 基线/兜底。23.9 kHz 的 Q=2 从 3.52 dB 降到 2.99 dB，Q=1 从 2.39 降到 0.97；
4 kHz 档 ≤0.107，691 Hz 档 ≤0.0032。

参数：4 个 boost（+12/+6/−6/−12 dB）× 4 个 Q（0.5/0.707/1/2）× 5 个 fc。

运行：``python qwqdsp/labs/match_filter/best_pole_fit/highshelf_2nd.py``
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

import power_match as PM

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "vicanek"))
import vicanek as V  # noqa: E402  （Vicanek 的 highshelf 作为 II 基线）

PI = PM.PI
"""圆周率。"""

FS = PM.FS
"""实验采样率。"""

FC_LIST = PM.FC_LIST
"""与 match_filter 其余实验一致的 5 个截止频率。"""

BOOST_LIST = (12.0, 6.0, -6.0, -12.0)
"""四个高频增益（dB）。"""

Q_LIST = (0.5, 1.0 / np.sqrt(2.0), 1.0, 2.0)
"""四个 Q。"""

RHO_MAX = 0.999
"""极点半径上限：防止极点贴单位圆（`p→0` 的近对消解）。"""


# --------------------------------------------------------------------- 功率函数
def analog_target(fc, fs, Q, db):
    """模拟高 shelf 的功率函数 ``|H|^2``。"""
    return V.shelf_target(db, 2.0 * PI * fc / fs, Q, "high")


def build_match_slope(fc, fs, Q, db, p, q):
    """0/fc 幅度 + fc 处功率斜率（fs/2 不约束）。

    数字功率的 0 点斜率恒为 0（``P(w)`` 是偶函数），与模拟一致；有效条件为
    ``P(0)``、``P(fc)``、``P'(fc)`` 三条，对 ``B`` 是线性的。``P'(fc)`` 等价于
    fc 处的 dB 斜率（因为 ``P(fc)`` 已相等）。给定极点后分子是闭式。
    """
    if not (p > 0.0 and q > 0.0 and p + q < 4.0):
        return None
    wc = 2.0 * PI * fc / fs
    alpha = PM.alpha_of(p, q)
    phi0 = PM.phi(np.array([0.0]))[0]
    phic = PM.phi(np.array([wc]))[0]
    sc = np.sin(wc / 2.0) ** 2
    dsc = np.sin(wc) / 2.0
    dphic = np.array([-1.0, 1.0, 4.0 * (1.0 - 2.0 * sc)]) * dsc
    target = analog_target(fc, fs, Q, db)
    h = 1e-6
    tp = (target(wc + h) - target(wc - h)) / (2.0 * h)
    t0 = target(0.0)
    tc = target(wc)
    ac = float(phic @ alpha)
    row_d = dphic * ac - phic * float(dphic @ alpha)
    rows = np.stack([phi0, phic, row_d])
    rhs = np.array([t0 * alpha[0], tc * ac, tp * ac * ac])
    beta = np.linalg.solve(rows, rhs)
    if not PM.feasible(beta):
        return None
    a = np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])
    if np.max(np.abs(np.roots([1.0, a[0], a[1]]))) > RHO_MAX:
        return None
    return PM.Match(PM.solve_b(beta), a, float(p), float(q), beta, alpha)


def max_error_db_slope(fc, fs, Q, db, p, q, n=20001):
    match = build_match_slope(fc, fs, Q, db, p, q)
    if match is None:
        return np.inf
    err, _, _ = error_db(match, fc, fs, Q, db, n)
    return float(np.max(np.abs(err)))


def match_from_ba(b, a):
    """由 ``(b, a)`` 构造 ``Match``（用于 Vicanek 的 II 解）。

    ``V.highshelf`` 返回的 ``a`` 含 ``a0 = 1``；``Match.a`` 只存 ``a1, a2``。
    """
    a = np.asarray(a, float)[1:]
    p = 1.0 + a[0] + a[1]
    q = 1.0 - a[0] + a[1]
    alpha = np.array([p * p, q * q, -4.0 * a[1]])
    beta = np.array([(b[0] + b[1] + b[2]) ** 2,
                     (b[0] - b[1] + b[2]) ** 2,
                     -4.0 * b[0] * b[2]])
    return PM.Match(np.asarray(b, float), a, float(p), float(q), beta, alpha)


def build_vicanek(fc, fs, Q, db):
    """Vicanek HighShelf：II 极点 + 目标拟合（含退化分支）。"""
    b, a = V.highshelf(2.0 * PI * fc / fs, Q, db)
    return match_from_ba(b, a)


def error_db(match, fc, fs, Q, db, n=20001):
    s = np.linspace(0.0, 1.0, n)
    ph = np.stack([1.0 - s, s, 4.0 * s * (1.0 - s)], axis=-1)
    hd2 = (ph @ match.beta) / (ph @ match.alpha)
    w = 2.0 * np.arcsin(np.sqrt(s))
    ha2 = analog_target(fc, fs, Q, db)(w)
    return 10.0 * np.log10(np.maximum(hd2, 1e-300) / np.maximum(ha2, 1e-300)), s, w


# --------------------------------------------------------------------- 极点
def best_design(fc, fs, Q, db):
    """0/fc 幅度 + fc 斜率（Hermite）+ 极点优化；fs/2 不约束，II 兜底。"""
    from scipy.optimize import minimize

    ii_match = build_vicanek(fc, fs, Q, db)
    ii_err = float(np.max(np.abs(error_db(ii_match, fc, fs, Q, db, 1000001)[0])))

    seeds = [(ii_match.p, ii_match.q)]
    grid = []
    for p in np.geomspace(1e-4, 3.9, 24):
        for q in np.linspace(1e-4, 3.9, 24):
            if p + q >= 4.0:
                continue
            value = max_error_db_slope(fc, fs, Q, db, p, q, 2001)
            if np.isfinite(value):
                grid.append((value, p, q))
    grid.sort()
    seeds += [(p, q) for _, p, q in grid[:3]]

    best = (np.inf, None)
    for seed in seeds:
        res = minimize(
            lambda x: max_error_db_slope(fc, fs, Q, db, x[0], x[1], 8001),
            np.array(seed), method="Nelder-Mead",
            options={"xatol": 1e-12, "fatol": 1e-12,
                     "maxiter": 2500, "maxfev": 2500})
        if res.fun < best[0]:
            best = (float(res.fun), res.x)
    if best[1] is None:
        return ii_match
    match = build_match_slope(fc, fs, Q, db, *best[1])
    if match is None:
        return ii_match

    ours = float(np.max(np.abs(error_db(match, fc, fs, Q, db, 1000001)[0])))
    if ii_err <= ours:
        return ii_match
    return match


def design_highshelf(fc, fs, Q, db):
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
    print(f"fs = {fs:g} Hz，全带 max|ΔdB|（0/fc + fc 斜率，Hermite；II = Vicanek）")
    print(f"{'boost':>6}{'Q':>8}{'fc':>8}{'p':>13}{'q':>12}"
          f"{'z1':>16}{'z2':>16}{'ours':>9}{'II':>9}")
    rows = []
    for db in BOOST_LIST:
        for Q in Q_LIST:
            for fc in FC_LIST:
                match = design_highshelf(fc, fs, Q, db)
                err, _, _ = error_db(match, fc, fs, Q, db, 1000001)
                ii_match = build_vicanek(fc, fs, Q, db)
                ii_err = float(np.max(np.abs(
                    error_db(ii_match, fc, fs, Q, db, 1000001)[0])))
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
            ha2 = analog_target(row["fc"], fs, Q, db)(w)
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
        ax_r.set_ylim(-30.0, 20.0)
        ax_r.set_yticks([-30.0, -24.0, -18.0, -12.0, -6.0, 0.0, 6.0, 12.0, 18.0])
        ax_r.grid(True, which="both", alpha=0.25)
        ax_r.legend(loc="lower left", fontsize=7)
        ax_r.set_xlabel("频率 / Hz")
        ax_e.set_ylim(1e-8, 1e2)
        ax_e.set_ylabel("|ΔdB|")
        ax_e.set_xlabel("f / (fs/2)")
        ax_e.grid(True, which="both", alpha=0.25)

    fig.suptitle(f"二阶高 shelf（boost = {db:+g} dB）：四个 Q 子类"
                 "（实线：最佳二阶节；虚线：模拟原型）")
    path = out_dir / f"highshelf_2nd_response_boost{db:+g}.png"
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
    path = out_dir / "highshelf_2nd_poles.png"
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
