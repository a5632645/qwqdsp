"""二阶模拟低通的二阶数字匹配：功率函数三点匹配 + 最佳极点。

模拟原型为 ``H(S) = 1/(S^2 + 2 zeta S + 1)``，``S = s/wc``、``zeta = 1/(2Q)``：

    |Ha|^2 = 1 / ((1 - r^2)^2 + (r/Q)^2),   r = f/fc.

在 ``0``、``fc``、``fs/2`` 三点匹配：

    P(0) = 1,   P(fc) = Q^2,   P(pi) = 1/((1 - r_pi^2)^2 + (r_pi/Q)^2).

分子由三点线性方程组唯一确定，分母用 ``p = A(1)``、``q = A(-1)`` 参数化；
两个自由值取等波纹（最小最大 ``|ΔdB|``）最优。与一阶不同，二阶模拟低通的
数字极点通常是复极点；初值取模拟极点的双线性像

    z = (1 + s T/2) / (1 - s T/2),   s T/2 = (-zeta ± j sqrt(1-zeta^2)) tan(wc/2).

对 Q = 0.5 / 0.707 / 1 / 2 各扫 5 个 fc，输出：

- 一张四个子类的图：每个 Q 一个子类，子类内是 5 个 fc 的幅频响应与误差；
- 一张 z 平面图：Q 是否改变最佳极点（空心 = 双线性极点，实心 = 最佳极点）。

注意：二阶低通在 Nyquist 的目标约为 ``T^2``，``B2`` 的抵消比一阶更严重；
最优 ``(p,q)`` 需要约 15 位有效数字，``b/a`` 则对 float32 舍入稳定。

运行：``python qwqdsp/labs/match_filter/best_pole_fit/lowpass_2nd.py``
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
    """二阶模拟低通的功率；``wc`` 为 fc 的角频率。"""
    r = np.asarray(w, dtype=float) / wc
    return 1.0 / ((1.0 - r * r) ** 2 + (r / Q) ** 2)


def build_match(fc, fs, Q, p, q):
    """给定分母两个自由值，做三点功率匹配；不可行时返回 ``None``。"""
    if not (p > 0.0 and q > 0.0 and p + q < 4.0):
        return None
    wc = 2.0 * PI * fc / fs
    alpha = PM.alpha_of(p, q)
    ph = PM.phi(np.array([0.0, wc, PI]))
    target = np.array([1.0, Q * Q, analog_power(PI, wc, Q)])
    beta = np.linalg.solve(ph, target * (ph @ alpha))
    if not PM.feasible(beta):
        return None
    a = np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])
    return PM.Match(PM.solve_b(beta), a, float(p), float(q), beta, alpha)


# --------------------------------------------------------------------- 误差
def error_db(match, fc, fs, Q, n=20001):
    """在 ``s = sin^2(w/2)`` 的均匀网格上返回 ``(ΔdB, s, w)``。"""
    s = np.linspace(0.0, 1.0, n)
    ph = np.stack([1.0 - s, s, 4.0 * s * (1.0 - s)], axis=-1)
    hd2 = (ph @ match.beta) / (ph @ match.alpha)
    w = 2.0 * np.arcsin(np.sqrt(s))
    ha2 = analog_power(w, 2.0 * PI * fc / fs, Q)
    return 10.0 * np.log10(np.maximum(hd2, 1e-300) / np.maximum(ha2, 1e-300)), s, w


def max_error_db(fc, fs, Q, p, q, n=20001):
    """不可行时返回 ``inf``。"""
    match = build_match(fc, fs, Q, p, q)
    if match is None:
        return np.inf
    err, _, _ = error_db(match, fc, fs, Q, n)
    return float(np.max(np.abs(err)))


# --------------------------------------------------------------------- 最佳极点
def _seed_pq(fc, fs, Q):
    poles = PM.bilinear_poles(fc, fs, Q)
    a1 = -float(np.real(poles.sum()))
    a2 = float(np.real(np.prod(poles)))
    return 1.0 + a1 + a2, 1.0 - a1 + a2


def best_denominator(fc, fs, Q):
    """等波纹意义下最优的分母两个自由值 ``(p, q)``。"""
    p0, q0 = _seed_pq(fc, fs, Q)
    best = None
    for dp in (0.0, 1e-4, 1e-3, 1e-2):
        for dq in (0.0, 1e-4, 1e-3, 1e-2):
            p1 = max(p0 * (1.0 + dp), 1e-14)
            q1 = max(q0 * (1.0 + dq), 1e-12)
            x, value = PM.refine(
                lambda p, q: max_error_db(fc, fs, Q, p, q, 2001), p1, q1,
                maxfev=600)
            if best is None or value < best[1]:
                best = (x, value)
    x, _ = PM.refine(lambda p, q: max_error_db(fc, fs, Q, p, q, 4001),
                     best[0][0], best[0][1])
    return float(x[0]), float(x[1])


def design_lowpass_2nd(fc, fs, Q):
    """二阶模拟低通的等波纹最佳二阶数字匹配。"""
    p, q = best_denominator(fc, fs, Q)
    match = build_match(fc, fs, Q, p, q)
    if match is None:
        raise RuntimeError("最佳分母落在可行域外")
    return match


# --------------------------------------------------------------------- 实验
def _format_pole(z):
    if abs(z.imag) < 1e-9:
        return f"{z.real:.6f}"
    return f"{z.real:.6f}{z.imag:+.6f}j"


def _report(fs):
    print(f"fs = {fs:g} Hz，全带 max|ΔdB|（0 <= f <= fs/2）")
    print(f"{'Q':>8}{'fc':>8}{'p':>15}{'q':>14}"
          f"{'z1':>18}{'z2':>18}{'full':>9}{'>= -60dB':>10}")
    print("z1/z2 是分母 1+a1 z^-1+a2 z^-2 的两个根（极点）："
          "a1=-(z1+z2), a2=z1*z2")

    rows = []
    for Q in Q_LIST:
        for fc in FC_LIST:
            match = design_lowpass_2nd(fc, fs, Q)
            err, _, w = error_db(match, fc, fs, Q, 1000001)
            ha = analog_power(w, 2.0 * PI * fc / fs, Q)
            full = float(np.max(np.abs(err)))
            practical = float(np.max(np.abs(err[ha >= 1e-6])))
            poles = match.poles()
            ph3 = PM.phi(np.array([0.0, 2.0 * PI * fc / fs, PI]))
            hd3 = (ph3 @ match.beta) / (ph3 @ match.alpha)
            target3 = np.array([1.0, Q * Q, analog_power(PI, 2.0 * PI * fc / fs, Q)])
            residual = float(np.max(np.abs(hd3 - target3)))
            rows.append({"Q": Q, "fc": fc, "match": match,
                         "poles": poles, "full": full, "practical": practical,
                         "residual": residual})
            print(f"{Q:8.4f}{fc:8.1f}{match.p:15.9g}{match.q:14.9g}"
                  f"{_format_pole(poles[0]):>18}{_format_pole(poles[1]):>18}"
                  f"{full:9.4f}{practical:10.4f}")
        print()
    print(f"\n三点匹配残差最大：{max(row['residual'] for row in rows):.2e}")
    return rows


def _figure_panels(rows, fs, out_dir):
    """四个 Q 子类；每个子类含 5 个 fc 的幅频响应与误差。"""
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
            ax_r.semilogx(f, 10.0 * np.log10(ha2), color=color, lw=3.4,
                          ls="--", alpha=0.30)
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
        ax_e.set_ylim(1e-6, 1.0)
        ax_e.set_ylabel("|ΔdB|")
        ax_e.set_xlabel("f / (fs/2)")
        ax_e.grid(True, which="both", alpha=0.25)

    fig.suptitle("二阶模拟低通：四个 Q 子类（实线：最佳二阶节；虚线：模拟原型）")
    path = out_dir / "lowpass_2nd_response.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _pole_metrics(poles):
    """复极点对取上半平面根；实极点取模较大的根（角度 0 或 pi）。"""
    z = poles[np.argmax(np.abs(poles.imag))]
    if abs(z.imag) > 1e-9:
        return abs(z), abs(float(np.angle(z)))
    z = poles[np.argmax(np.abs(poles))]
    return abs(z), (PI if z.real < 0.0 else 0.0)


def _figure_poles(rows, fs, out_dir):
    """每个 Q 一个子图：极点半径与极点角度随 fc 的变化。"""
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
        radius = []
        angle_deg = []
        for row in group:
            r, theta = _pole_metrics(row["poles"])
            radius.append(r)
            angle_deg.append(np.degrees(theta))
        bilinear = [_pole_metrics(PM.bilinear_poles(fc, fs, Q)) for fc in fcs]
        radius_bil = [item[0] for item in bilinear]
        angle_bil = [np.degrees(item[1]) for item in bilinear]

        line_r, = ax.semilogx(fcs, radius, "o-", color="C0", lw=1.4, ms=5,
                              label="半径 |z|（最佳）")
        line_rb, = ax.semilogx(fcs, radius_bil, "o--", color="C0", lw=1.0,
                               ms=4, mfc="none", alpha=0.75,
                               label="半径 |z|（双线性）")
        ax.set_ylim(0.0, 1.05)
        ax.set_ylabel("极点半径 |z|", color="C0")
        ax.tick_params(axis="y", colors="C0")
        ax.set_xlabel("fc / Hz")
        ax.set_title(f"Q = {Q:g}")
        ax.grid(True, which="both", alpha=0.25)

        ax_angle = ax.twinx()
        line_a, = ax_angle.semilogx(fcs, angle_deg, "s-", color="C3", lw=1.4,
                                    ms=5, label="角度 θ（最佳）")
        line_ab, = ax_angle.semilogx(fcs, angle_bil, "s--", color="C3", lw=1.0,
                                     ms=4, mfc="none", alpha=0.75,
                                     label="角度 θ（双线性）")
        ax_angle.set_ylim(0.0, 180.0)
        ax_angle.set_yticks([0.0, 45.0, 90.0, 135.0, 180.0])
        ax_angle.set_ylabel("极点角度 / deg", color="C3")
        ax_angle.tick_params(axis="y", colors="C3")

        ax.legend(handles=[line_r, line_rb, line_a, line_ab],
                  loc="center left", fontsize=7)

    fig.suptitle("最佳极点的半径与角度（实线/实心：最佳；虚线/空心：双线性）")
    fig.tight_layout()
    path = out_dir / "lowpass_2nd_poles.png"
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
