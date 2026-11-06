"""二阶模拟带通的二阶数字匹配：三点功率匹配 + 极点优化。

模拟原型取**归一化**形式

    H(S) = (S/Q) / (S^2 + S/Q + 1),   S = s/wc,

峰值 ``|H(fc)|^2 = 1``，低频渐近 ``|H|^2 ~ (f/fc)^2 / Q^2``。功率

    |Ha|^2 = (r^2/Q^2) / ((1 - r^2)^2 + (r/Q)^2),   r = f/fc,

在 ``0``、``fc``、``fs/2`` 三点匹配：

    P(0) = 0,   P(fc) = 1,   P(pi) = (r_pi^2/Q^2)/((1-r_pi^2)^2 + (r_pi/Q)^2).

分子由三点线性方程组唯一确定，分母用 ``p = A(1)``、``q = A(-1)`` 参数化；
两个自由值按全带等波纹（最小最大 ``|ΔdB|``）优化。作为基线同时计算
Vicanek 的脉冲响应不变（II）极点，用同一三点分子。

未归一化形式 ``H = S/(S^2 + S/Q + 1)``（峰值 ``Q^2``）只差一个幅度因子
``Q``：把归一化结果乘 ``Q`` 即得未归一化，极点与误差完全相同。

运行：``python qwqdsp/labs/match_filter/best_pole_fit/bandpass_2nd.py``
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

Q_LIST = (0.5, 1.0 / np.sqrt(2.0), 1.0, 2.0)
"""四个 Q（0.5 为临界阻尼，0.707 为 Butterworth，2 有谐振峰）。"""


# --------------------------------------------------------------------- 功率函数
def analog_power(w, wc, Q):
    """二阶模拟带通的功率（归一化峰值 1）。"""
    r = np.asarray(w, dtype=float) / wc
    return (r * r / (Q * Q)) / ((1.0 - r * r) ** 2 + (r / Q) ** 2)


def build_match(fc, fs, Q, p, q):
    """给定分母两个自由值，做三点功率匹配。"""
    if not (p > 0.0 and q > 0.0 and p + q < 4.0):
        return None
    wc = 2.0 * PI * fc / fs
    alpha = PM.alpha_of(p, q)
    ph = PM.phi(np.array([0.0, wc, PI]))
    target = np.array([0.0, 1.0, analog_power(PI, wc, Q)])
    beta = np.linalg.solve(ph, target * (ph @ alpha))
    if not PM.feasible(beta):
        return None
    a = np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])
    return PM.Match(PM.solve_b(beta), a, float(p), float(q), beta, alpha)


def error_db(match, fc, fs, Q, n=20001):
    s = np.linspace(0.0, 1.0, n)
    ph = np.stack([1.0 - s, s, 4.0 * s * (1.0 - s)], axis=-1)
    hd2 = (ph @ match.beta) / (ph @ match.alpha)
    w = 2.0 * np.arcsin(np.sqrt(s))
    ha2 = analog_power(w, 2.0 * PI * fc / fs, Q)
    return 10.0 * np.log10(np.maximum(hd2, 1e-300) / np.maximum(ha2, 1e-300)), s, w


def max_error_db(fc, fs, Q, p, q, n=20001):
    match = build_match(fc, fs, Q, p, q)
    if match is None:
        return np.inf
    err, _, _ = error_db(match, fc, fs, Q, n)
    return float(np.max(np.abs(err)))


# --------------------------------------------------------------------- 极点
def vicanek_pq(fc, fs, Q):
    """Vicanek 的脉冲响应不变极点，写成 p/q 参数。"""
    zeta = 1.0 / (2.0 * Q)
    wc = 2.0 * PI * fc / fs
    e = np.exp(-zeta * wc)
    mu = np.sqrt(max(0.0, 1.0 - zeta * zeta))
    a1 = -2.0 * e * np.cos(mu * wc)
    a2 = e * e
    return 1.0 + a1 + a2, 1.0 - a1 + a2


def best_denominator(fc, fs, Q):
    from scipy.optimize import minimize

    p_ii, q_ii = vicanek_pq(fc, fs, Q)
    candidates = [(max_error_db(fc, fs, Q, p_ii, q_ii, 4001), p_ii, q_ii)]
    for p in np.geomspace(1e-4, 3.9, 36):
        for q in np.linspace(1e-4, 3.9, 36):
            if p + q >= 4.0:
                continue
            value = max_error_db(fc, fs, Q, p, q, 2001)
            if np.isfinite(value):
                candidates.append((value, p, q))
    candidates.sort()

    best = (np.inf, None)
    for _, p, q in candidates[:4]:
        res = minimize(
            lambda x: max_error_db(fc, fs, Q, x[0], x[1], 8001),
            np.array([p, q]), method="Nelder-Mead",
            options={"xatol": 1e-12, "fatol": 1e-12,
                     "maxiter": 2000, "maxfev": 1500})
        if res.fun < best[0]:
            best = (float(res.fun), res.x)
    if best[1] is None:
        return p_ii, q_ii

    matched = build_match(fc, fs, Q, *best[1])
    ours = float(np.max(np.abs(error_db(matched, fc, fs, Q, 1000001)[0])))
    iim = build_match(fc, fs, Q, p_ii, q_ii)
    theirs = float(np.max(np.abs(error_db(iim, fc, fs, Q, 1000001)[0])))
    if theirs <= ours:
        return p_ii, q_ii
    return float(best[1][0]), float(best[1][1])


def design_bandpass_2nd(fc, fs, Q):
    p, q = best_denominator(fc, fs, Q)
    match = build_match(fc, fs, Q, p, q)
    if match is None:
        raise RuntimeError("没有可行的极点")
    return match


# --------------------------------------------------------------------- 实验
def _format_pole(z):
    if abs(z.imag) < 1e-9:
        return f"{z.real:.6f}"
    return f"{z.real:.6f}{z.imag:+.6f}j"


def _report(fs):
    print(f"fs = {fs:g} Hz，全带 max|ΔdB|（0 / fc / fs/2 三点，峰值 1）")
    print(f"{'Q':>8}{'fc':>8}{'p':>13}{'q':>12}"
          f"{'z1':>18}{'z2':>18}{'ours':>9}{'II':>9}{'>=-60dB':>9}")
    print("z1/z2 是分母 1+a1 z^-1+a2 z^-2 的两个根（极点）："
          "a1=-(z1+z2), a2=z1*z2")

    rows = []
    for Q in Q_LIST:
        for fc in FC_LIST:
            match = design_bandpass_2nd(fc, fs, Q)
            err, _, w = error_db(match, fc, fs, Q, 1000001)
            ours = float(np.max(np.abs(err)))
            ha = analog_power(w, 2.0 * PI * fc / fs, Q)
            practical = float(np.max(np.abs(err[ha >= 1e-6])))
            p_ii, q_ii = vicanek_pq(fc, fs, Q)
            ii_match = build_match(fc, fs, Q, p_ii, q_ii)
            err_ii, _, _ = error_db(ii_match, fc, fs, Q, 1000001)
            theirs = float(np.max(np.abs(err_ii)))
            poles = match.poles()
            rows.append({"Q": Q, "fc": fc, "match": match, "poles": poles,
                         "ours": ours, "ii": theirs, "practical": practical,
                         "ii_match": ii_match})
            print(f"{Q:8.4f}{fc:8.1f}{match.p:13.6g}{match.q:12.6g}"
                  f"{_format_pole(poles[0]):>18}{_format_pole(poles[1]):>18}"
                  f"{ours:9.4f}{theirs:9.4f}{practical:9.4f}")
        print()
    return rows


def _figure_panels(rows, fs, out_dir):
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
        group = [row for row in rows if row["Q"] == Q]
        for j, row in enumerate(group):
            color = f"C{j}"
            ph = PM.phi(w)
            hd2 = (ph @ row["match"].beta) / (ph @ row["match"].alpha)
            ha2 = analog_power(w, 2.0 * PI * row["fc"] / fs, Q)
            ax_r.semilogx(f, 10.0 * np.log10(np.maximum(ha2, 1e-300)),
                          color=color, lw=3.4, ls="--", alpha=0.30)
            ax_r.semilogx(f, 10.0 * np.log10(np.maximum(hd2, 1e-300)),
                          color=color, lw=1.4, label=f"{row['fc']:g} Hz")
            err, _, ww = error_db(row["match"], row["fc"], fs, Q, 200001)
            ax_e.semilogy(ww / PI, np.maximum(np.abs(err), 1e-12),
                          color=color, lw=1.2)
        ax_r.set_title(f"Q = {Q:g}")
        ax_r.set_ylabel("幅度 / dB")
        ax_r.set_xlim(10.0, fs / 2.0)
        ax_r.set_ylim(-80.0, 15.0)
        ax_r.set_yticks([-80.0, -60.0, -40.0, -20.0, 0.0, 15.0])
        ax_r.grid(True, which="both", alpha=0.25)
        ax_r.legend(loc="lower left", fontsize=7)
        ax_r.set_xlabel("频率 / Hz")
        ax_e.set_ylim(1e-8, 1.0)
        ax_e.set_ylabel("|ΔdB|")
        ax_e.set_xlabel("f / (fs/2)")
        ax_e.grid(True, which="both", alpha=0.25)

    fig.suptitle("二阶模拟带通：四个 Q 子类（实线：最佳二阶节；虚线：模拟原型）")
    path = out_dir / "bandpass_2nd_response.png"
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
    for ax, Q in zip(axes.ravel(), Q_LIST):
        group = sorted((row for row in rows if row["Q"] == Q),
                       key=lambda row: row["fc"])
        fcs = [row["fc"] for row in group]
        radius, angle_deg = [], []
        for row in group:
            r, theta = _pole_metrics(row["poles"])
            radius.append(r)
            angle_deg.append(np.degrees(theta))
        ii = [_pole_metrics(row["ii_match"].poles()) for row in group]
        radius_ii = [item[0] for item in ii]
        angle_ii = [np.degrees(item[1]) for item in ii]

        line_r, = ax.semilogx(fcs, radius, "o-", color="C0", lw=1.4, ms=5,
                              label="半径 |z|（最佳）")
        line_rb, = ax.semilogx(fcs, radius_ii, "o--", color="C0", lw=1.0,
                               ms=4, mfc="none", alpha=0.75,
                               label="半径 |z|（II）")
        ax.set_ylim(0.0, 1.05)
        ax.set_ylabel("极点半径 |z|", color="C0")
        ax.tick_params(axis="y", colors="C0")
        ax.set_xlabel("fc / Hz")
        ax.set_title(f"Q = {Q:g}")
        ax.grid(True, which="both", alpha=0.25)

        ax_angle = ax.twinx()
        line_a, = ax_angle.semilogx(fcs, angle_deg, "s-", color="C3", lw=1.4,
                                    ms=5, label="角度 θ（最佳）")
        line_ab, = ax_angle.semilogx(fcs, angle_ii, "s--", color="C3", lw=1.0,
                                     ms=4, mfc="none", alpha=0.75,
                                     label="角度 θ（II）")
        ax_angle.set_ylim(0.0, 180.0)
        ax_angle.set_yticks([0.0, 45.0, 90.0, 135.0, 180.0])
        ax_angle.set_ylabel("极点角度 / deg", color="C3")
        ax_angle.tick_params(axis="y", colors="C3")
        ax.legend(handles=[line_r, line_rb, line_a, line_ab],
                  loc="center left", fontsize=7)

    fig.suptitle("最佳极点的半径与角度（实线/实心：最佳；虚线/空心：II）")
    fig.tight_layout()
    path = out_dir / "bandpass_2nd_poles.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main():
    rows = _report(FS)
    out_dir = Path(__file__).parent / "output"
    p1 = _figure_panels(rows, FS, out_dir)
    p2 = _figure_poles(rows, FS, out_dir)
    print(f"图：{p1}\n图：{p2}")


if __name__ == "__main__":
    main()
