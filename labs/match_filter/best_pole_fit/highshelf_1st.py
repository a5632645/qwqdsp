"""一阶模拟高 shelf 的二阶数字匹配：功率函数三点匹配 + 最佳极点。

与低通/高通同一套框架（见 ``power_match``），模拟原型取

    H(S) = (1 + sqrt(g) S) / (1 + S/sqrt(g)),   S = s/wc,  g = 10^(db/20),

即 ``H(0) = 1``、``H(inf) = g``、``|H(fc)|^2 = g``。功率函数

    |Ha|^2 = (1 + g r^2) / (1 + r^2/g),   r = f/fc,

在 ``0``、``fc``、``fs/2`` 三点匹配：

    P(0) = 1,   P(fc) = g,   P(pi) = (1 + g r_pi^2)/(1 + r_pi^2/g).

分子由三点线性方程组唯一确定，分母用 ``p = A(1)``、``q = A(-1)`` 参数化；
两个自由值取等波纹（最小最大 ``|ΔdB|``）最优。最佳极点的初值取模拟极点
``wc*sqrt(g)`` 的双线性像。

对四个 boost（+12 / +6 / -6 / -12 dB）各扫 5 个 fc，输出：

- 一张四个子类的图：每个 boost 一个子类，子类内是 5 个 fc 的幅频响应与误差；
- 一张极点图：boost 是否改变最佳极点（第一极点随 ``wc*sqrt(g)`` 变，第二
  极点在小 fc 下几乎不变，约 -0.556）。

运行：``python qwqdsp/labs/match_filter/best_pole_fit/highshelf_1st.py``
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

DB_LIST = (12.0, 6.0, -6.0, -12.0)
"""四个 shelf boost（dB）。"""


# --------------------------------------------------------------------- 功率函数
def analog_power(w, wc, db):
    """一阶模拟高 shelf 的功率；``wc`` 为 fc 的角频率。"""
    g = 10.0 ** (db / 20.0)
    r = np.asarray(w, dtype=float) / wc
    return (1.0 + g * r * r) / (1.0 + r * r / g)


def build_match(fc, fs, db, p, q):
    """给定分母两个自由值，做三点功率匹配；不可行时返回 ``None``。"""
    if not (p > 0.0 and q > 0.0 and p + q < 4.0):
        return None
    g = 10.0 ** (db / 20.0)
    wc = 2.0 * PI * fc / fs
    alpha = PM.alpha_of(p, q)
    ph = PM.phi(np.array([0.0, wc, PI]))
    target = np.array([1.0, g, analog_power(PI, wc, db)])
    beta = np.linalg.solve(ph, target * (ph @ alpha))
    if not PM.feasible(beta):
        return None
    a = np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])
    return PM.Match(PM.solve_b(beta), a, float(p), float(q), beta, alpha)


# --------------------------------------------------------------------- 误差
def error_db(match, fc, fs, db, n=20001):
    """在 ``s = sin^2(w/2)`` 的均匀网格上返回 ``(ΔdB, s, w)``。"""
    s = np.linspace(0.0, 1.0, n)
    ph = np.stack([1.0 - s, s, 4.0 * s * (1.0 - s)], axis=-1)
    hd2 = (ph @ match.beta) / (ph @ match.alpha)
    w = 2.0 * np.arcsin(np.sqrt(s))
    ha2 = analog_power(w, 2.0 * PI * fc / fs, db)
    return 10.0 * np.log10(np.maximum(hd2, 1e-300) / ha2), s, w


def max_error_db(fc, fs, db, p, q, n=20001):
    """不可行时返回 ``inf``。"""
    match = build_match(fc, fs, db, p, q)
    if match is None:
        return np.inf
    err, _, _ = error_db(match, fc, fs, db, n)
    return float(np.max(np.abs(err)))


# --------------------------------------------------------------------- 最佳极点
def shelf_pole(fc, fs, db):
    """模拟极点 ``wc*sqrt(g)`` 的双线性像；作为第一极点的初值。"""
    g = 10.0 ** (db / 20.0)
    wc = 2.0 * PI * fc / fs
    a = np.sqrt(g) * np.tan(wc / 2.0)
    return (1.0 - a) / (1.0 + a)


def best_denominator(fc, fs, db):
    """等波纹意义下最优的分母两个自由值 ``(p, q)``。"""
    return PM.best_denominator(
        fc, fs, lambda p, q, n: max_error_db(fc, fs, db, p, q, n),
        r1=shelf_pole(fc, fs, db))


def design_highshelf_1st(fc, fs, db):
    """一阶模拟高 shelf 的等波纹最佳二阶数字匹配。"""
    p, q = best_denominator(fc, fs, db)
    match = build_match(fc, fs, db, p, q)
    if match is None:
        raise RuntimeError("最佳分母落在可行域外")
    return match


# --------------------------------------------------------------------- 实验
def _report(fs):
    print(f"fs = {fs:g} Hz，全带 max|ΔdB|（0 <= f <= fs/2）")
    header = (f"{'boost':>7}{'fc':>8}{'p':>13}{'q':>12}"
              f"{'z1':>10}{'z2':>10}{'maxdB':>11}")
    print(header)
    print("z1/z2 是分母 1+a1 z^-1+a2 z^-2 的两个根（极点）："
          "a1=-(z1+z2), a2=z1*z2")

    rows = []
    for db in DB_LIST:
        for fc in FC_LIST:
            match = design_highshelf_1st(fc, fs, db)
            err, _, _ = error_db(match, fc, fs, db, 400001)
            poles = np.sort(match.poles())
            best = float(np.max(np.abs(err)))
            rows.append({"db": db, "fc": fc, "match": match,
                         "poles": poles, "max_db": best})
            print(f"{db:+7.1f}{fc:8.1f}{match.p:13.9f}{match.q:12.9f}"
                  f"{poles[0]:10.6f}{poles[1]:10.6f}{best:11.3e}")
        print()
    return rows


def _figure_panels(rows, fs, out_dir):
    """四个 boost 子类；每个子类含 5 个 fc 的幅频响应与误差。"""
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
    for index, db in enumerate(DB_LIST):
        cell = outer[index // 2, index % 2].subgridspec(
            2, 1, height_ratios=(2.1, 1.1), hspace=0.25)
        ax_r = fig.add_subplot(cell[0])
        ax_e = fig.add_subplot(cell[1])
        group = [row for row in rows if row["db"] == db]
        for j, row in enumerate(group):
            color = f"C{j}"
            ph = PM.phi(w)
            hd2 = (ph @ row["match"].beta) / (ph @ row["match"].alpha)
            ha2 = analog_power(w, 2.0 * PI * row["fc"] / fs, db)
            ax_r.semilogx(f, 10.0 * np.log10(ha2), color=color, lw=3.4,
                          ls="--", alpha=0.30)
            ax_r.semilogx(f, 10.0 * np.log10(np.maximum(hd2, 1e-300)),
                          color=color, lw=1.4, label=f"{row['fc']:g} Hz")
            err, _, ww = error_db(row["match"], row["fc"], fs, db, 200001)
            ax_e.semilogy(ww / PI, np.maximum(np.abs(err), 1e-12),
                          color=color, lw=1.2)
        ax_r.set_title(f"boost = {db:+g} dB")
        ax_r.set_ylabel("幅度 / dB")
        ax_r.set_xlim(10.0, fs / 2.0)
        ax_r.set_ylim(min(0.0, db) - 1.5, max(0.0, db) + 1.5)
        ax_r.grid(True, which="both", alpha=0.25)
        ax_r.legend(loc="lower right", fontsize=7)
        ax_r.set_xlabel("频率 / Hz")
        ax_e.set_ylim(1e-8, 1.0)
        ax_e.set_ylabel("|ΔdB|")
        ax_e.set_xlabel("f / (fs/2)")
        ax_e.grid(True, which="both", alpha=0.25)

    fig.suptitle("高 shelf：四个 boost 子类（实线：最佳二阶节；虚线：模拟原型）")
    path = out_dir / "highshelf_1st_response.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _figure_poles(rows, fs, out_dir):
    """boost 对最佳极点的影响：第一/第二极点随 fc 的轨迹。"""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    out_dir.mkdir(exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.6))
    for index, db in enumerate(DB_LIST):
        group = sorted((row for row in rows if row["db"] == db),
                       key=lambda row: row["fc"])
        first, second, seeds = [], [], []
        for row in group:
            r1 = shelf_pole(row["fc"], fs, db)
            poles = row["poles"]
            i = int(np.argmin(np.abs(poles - r1)))
            first.append(poles[i])
            second.append(poles[1 - i])
            seeds.append(r1)
        color = f"C{index}"
        fcs = [row["fc"] for row in group]
        ax1.semilogx(fcs, seeds, color=color, lw=0.9, ls="--", alpha=0.7)
        ax1.semilogx(fcs, first, "o-", color=color, lw=1.3, ms=4,
                     label=f"{db:+g} dB")
        ax2.semilogx(fcs, second, "s-", color=color, lw=1.3, ms=4,
                     label=f"{db:+g} dB")
    ax1.axhline(0.0, color="0.6", lw=0.8, ls=":")
    ax1.set_xlabel("fc / Hz")
    ax1.set_ylabel("第 1 极点 z")
    ax1.set_ylim(-1.02, 1.02)
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(fontsize=8)
    ax1.set_title("第 1 极点（虚线：模拟极点 wc√g 的双线性像）")

    ax2.axhline(-0.556, color="0.6", lw=0.8, ls=":")
    ax2.set_xlabel("fc / Hz")
    ax2.set_ylabel("第 2 极点 z")
    ax2.set_ylim(-1.02, 1.02)
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(fontsize=8)
    ax2.set_title("第 2 极点（点线：z = -0.556）")
    fig.tight_layout()

    path = out_dir / "highshelf_1st_poles.png"
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
