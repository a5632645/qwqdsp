"""二阶模拟 notch：DC + fc 两点匹配 + 极点优化。

模拟原型 ``H(S) = (S^2 + 1)/(S^2 + S/Q + 1)``，``S = s/wc``：

    |Ha|^2 = (1 - r^2)^2 / ((1 - r^2)^2 + (r/Q)^2),   r = f/fc,

DC 与 Nyquist 都是 1，``fc`` 处为 0。三点（0、fc、fs/2）匹配不可实现：
缺口要求 fc 处功率为 0，三点解把 ``B2`` 推到可行性边界外。因此按 Vicanek
的 Notch 结构做两点：

    B0 = A0,   B2 = -A0/(4 s_c^2),   B1 = A0 + 4(2 s_c - 1) B2,

对应分子零点在单位圆 ``e^{±j wc}``、DC 归一为 1、fc 处严格为零。``B`` 正好
落在可行性边界 ``W^2 + B2 = 0`` 上（``b0 = b2``）。分母的两个自由值按全带
等波纹（最小最大 ``|ΔdB|``）优化；同时用同一分子计算 Vicanek 的脉冲响应
不变（II）极点作为基线。

结果：fc <= 4 kHz 全带峰值 1e-4 … 0.06 dB（II 0.000 … 0.17），23.9 kHz
档 1.7 dB（II 40.5 dB）。

运行：``python qwqdsp/labs/match_filter/best_pole_fit/notch_2nd.py``
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
"""四个 Q（0.5 为临界阻尼，0.707 为 Butterworth，2 有窄缺口）。"""


# --------------------------------------------------------------------- 功率函数
def analog_power(w, wc, Q):
    """二阶模拟 notch 的功率（DC/Nyquist 为 1，fc 为 0）。"""
    r = np.asarray(w, dtype=float) / wc
    a = (1.0 - r * r) ** 2
    return a / (a + (r / Q) ** 2)


def build_match(fc, fs, Q, p, q):
    """Vicanek Notch 的分子（单位圆零点）+ 给定分母。"""
    if not (p > 0.0 and q > 0.0 and p + q < 4.0):
        return None
    wc = 2.0 * PI * fc / fs
    sc = np.sin(wc / 2.0) ** 2
    alpha = PM.alpha_of(p, q)
    b0p = alpha[0]
    b2p = -alpha[0] / (4.0 * sc * sc)
    b1p = alpha[0] + 4.0 * (2.0 * sc - 1.0) * b2p
    if not (b0p > 0.0 and b1p > 0.0):
        return None
    beta = np.array([b0p, b1p, b2p])
    W = (np.sqrt(beta[0]) + np.sqrt(beta[1])) / 2.0
    rad = W * W + beta[2]
    if rad < -1e-8 * max(W * W, abs(beta[2]), 1e-300):
        return None
    b0 = (W + np.sqrt(max(rad, 0.0))) / 2.0
    b1 = (np.sqrt(beta[0]) - np.sqrt(beta[1])) / 2.0
    b2 = -beta[2] / (4.0 * b0)
    a = np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])
    return PM.Match(np.array([b0, b1, b2]), a, float(p), float(q), beta, alpha)


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


def design_notch_2nd(fc, fs, Q):
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
    print(f"fs = {fs:g} Hz，全带 max|ΔdB|（DC + fc 缺口两点）")
    print(f"{'Q':>8}{'fc':>8}{'p':>13}{'q':>12}"
          f"{'z1':>18}{'z2':>18}{'ours':>10}{'II':>10}")
    print("z1/z2 是分母 1+a1 z^-1+a2 z^-2 的两个根（极点）："
          "a1=-(z1+z2), a2=z1*z2")

    rows = []
    for Q in Q_LIST:
        for fc in FC_LIST:
            match = design_notch_2nd(fc, fs, Q)
            err, _, _ = error_db(match, fc, fs, Q, 1000001)
            ours = float(np.max(np.abs(err)))
            p_ii, q_ii = vicanek_pq(fc, fs, Q)
            ii_match = build_match(fc, fs, Q, p_ii, q_ii)
            err_ii, _, _ = error_db(ii_match, fc, fs, Q, 1000001)
            theirs = float(np.max(np.abs(err_ii)))
            poles = match.poles()
            rows.append({"Q": Q, "fc": fc, "match": match, "poles": poles,
                         "ours": ours, "ii": theirs, "ii_match": ii_match})
            print(f"{Q:8.4f}{fc:8.1f}{match.p:13.6g}{match.q:12.6g}"
                  f"{_format_pole(poles[0]):>18}{_format_pole(poles[1]):>18}"
                  f"{ours:10.4f}{theirs:10.4f}")
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
        ax_e.set_ylim(1e-8, 1e2)
        ax_e.set_ylabel("|ΔdB|")
        ax_e.set_xlabel("f / (fs/2)")
        ax_e.grid(True, which="both", alpha=0.25)

    fig.suptitle("二阶模拟 notch：四个 Q 子类（实线：最佳二阶节；虚线：模拟原型）")
    path = out_dir / "notch_2nd_response.png"
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
    path = out_dir / "notch_2nd_poles.png"
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
