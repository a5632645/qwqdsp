"""一阶模拟高通的二阶数字匹配：功率函数三点匹配 + 最佳极点。

与低通同一套框架（见 ``power_match``）：数字二阶节的功率函数在
``s = sin^2(w/2)`` 上是两个二次型之比，在 ``0``、``fc``、``fs/2`` 三点令它
等于一阶模拟高通

    |Ha|^2 = (w/wc)^2 / (1 + (w/wc)^2),

即 ``P(0) = 0``、``P(fc) = 1/2``、``P(pi) = U = (pi/wc)^2/(1+(pi/wc)^2)``。
分子由三点线性方程组唯一确定，分母用

    p = A(1) = 1 + a1 + a2,
    q = A(-1) = 1 - a1 + a2

参数化；稳定的两个自由值 ``p, q`` 取等波纹（最小最大 ``|ΔdB|``）最优。

## 与低通的差别

低通的模拟功率在 Nyquist 处趋于 0，而双线性/数字实现也在 Nyquist 趋于 0，
三点匹配在低 fc 时留下一条 0.12 dB 的等波纹；高通在高频趋于 1，数字与模拟
的功率在高频只差 ``1/(w/wc)^2`` 量级，因此三点匹配在低 fc 下几乎精确：

- fc = 20 Hz：全带 ``max|ΔdB|`` 约 ``1.5e-7``；
- fc = 4 kHz：约 ``6e-3`` dB；
- fc = 23.9 kHz（截止贴近 Nyquist）：约 ``0.213`` dB。

最佳第二极点也与低通不同：fc <= 4 kHz 时稳定在 ``z ~ -0.556``（低通是
``-0.696``），且它几乎与一个零点相消，使二阶节在低频段有效退化为一阶；
只保留双线性极点（分子仍三点匹配）时误差大 4–32 倍。

运行：``python qwqdsp/labs/match_filter/best_pole_fit/highpass_1st.py``
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


# --------------------------------------------------------------------- 功率函数
def analog_power(w, wc):
    """一阶模拟高通功率 ``(w/wc)^2 / (1 + (w/wc)^2)``。"""
    r = np.asarray(w, dtype=float) / wc
    return r * r / (1.0 + r * r)


def build_match(fc, fs, p, q):
    """给定分母两个自由值，做三点功率匹配；不可行时返回 ``None``。"""
    if not (p > 0.0 and q > 0.0 and p + q < 4.0):
        return None
    wc = 2.0 * PI * fc / fs
    alpha = PM.alpha_of(p, q)
    ph = PM.phi(np.array([0.0, wc, PI]))
    target = np.array([0.0, 0.5, analog_power(PI, wc)])
    beta = np.linalg.solve(ph, target * (ph @ alpha))
    if not PM.feasible(beta):
        return None
    a = np.array([(p - q) / 2.0, (p + q) / 2.0 - 1.0])
    return PM.Match(PM.solve_b(beta), a, float(p), float(q), beta, alpha)


# --------------------------------------------------------------------- 误差
def error_db(match, fc, fs, n=20001):
    """在 ``s = sin^2(w/2)`` 的均匀网格上返回 ``(ΔdB, s, w)``。

    ``s = 0`` 处数字与模拟功率都为 0，分子分母同时取下限，记作 0 dB。
    """
    s = np.linspace(0.0, 1.0, n)
    ph = np.stack([1.0 - s, s, 4.0 * s * (1.0 - s)], axis=-1)
    hd2 = (ph @ match.beta) / (ph @ match.alpha)
    w = 2.0 * np.arcsin(np.sqrt(s))
    ha2 = analog_power(w, 2.0 * PI * fc / fs)
    ratio = np.maximum(hd2, 1e-300) / np.maximum(ha2, 1e-300)
    return 10.0 * np.log10(ratio), s, w


def max_error_db(fc, fs, p, q, n=20001):
    """不可行时返回 ``inf``。"""
    match = build_match(fc, fs, p, q)
    if match is None:
        return np.inf
    err, _, _ = error_db(match, fc, fs, n)
    return float(np.max(np.abs(err)))


# --------------------------------------------------------------------- 最佳极点
def best_denominator(fc, fs=FS):
    """等波纹意义下最优的分母两个自由值 ``(p, q)``。"""
    return PM.best_denominator(fc, fs, lambda p, q, n: max_error_db(fc, fs, p, q, n))


def design_highpass_1st(fc, fs=FS):
    """一阶模拟高通的等波纹最佳二阶数字匹配。"""
    p, q = best_denominator(fc, fs)
    match = build_match(fc, fs, p, q)
    if match is None:
        raise RuntimeError("最佳分母落在可行域外")
    return match


# --------------------------------------------------------------------- 实验
def _report(fs):
    print(f"fs = {fs:g} Hz，全带 max|ΔdB|（0 <= f <= fs/2）")
    print(f"{'fc':>8}{'p':>13}{'q':>12}{'z1':>10}{'z2':>10}"
          f"{'best':>11}{'bilin':>11}")
    print("z1/z2 是分母 1+a1 z^-1+a2 z^-2 的两个根（极点）："
          "a1=-(z1+z2), a2=z1*z2")

    rows = []
    for fc in FC_LIST:
        match = design_highpass_1st(fc, fs)
        err, _, _ = error_db(match, fc, fs, 400001)
        poles = np.sort(match.poles())
        best = float(np.max(np.abs(err)))

        r1 = PM.bilinear_pole(fc, fs)
        bilin = max_error_db(fc, fs, 1.0 - r1, 1.0 + r1, 100001)

        rows.append((fc, match, poles, best, bilin))
        print(f"{fc:8.1f}{match.p:13.9f}{match.q:12.9f}{poles[0]:10.6f}"
              f"{poles[1]:10.6f}{best:11.3e}{bilin:11.3e}")

    print("\n三点匹配残差与等波纹极值（最佳二阶节）")
    for fc, match, _, _, _ in rows:
        wc = 2.0 * PI * fc / fs
        ph = PM.phi(np.array([0.0, wc, PI]))
        hd = (ph @ match.beta) / (ph @ match.alpha)
        err, _, _ = error_db(match, fc, fs, 1000001)
        idx = np.where(np.diff(np.sign(np.diff(err))) != 0)[0]
        top = idx[np.argsort(-np.abs(err[idx]))[:3]]
        target = np.array([0.0, 0.5, analog_power(PI, wc)])
        print(f"  fc={fc:7.1f}  三点 |H|^2 = {hd[0]:.1e} {hd[1]:.12f} {hd[2]:.12g}"
              f"  残差 {np.max(np.abs(hd - target)):.2e}"
              f"  等波纹 {np.array2string(err[top], formatter={'float_kind': lambda x: f'{x:.2e}'})}")
    return rows


def _figures(rows, fs, out_dir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    out_dir.mkdir(exist_ok=True)

    # 图 1：幅频响应（5 个 fc）+ 误差幅值
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12.5, 4.8))
    f = np.geomspace(10.0, fs / 2.0, 4001)
    w = 2.0 * PI * f / fs
    for i, (fc, m, _, _, _) in enumerate(rows):
        color = f"C{i}"
        ax0.semilogx(f, 10 * np.log10(np.maximum(analog_power(w, 2 * PI * fc / fs), 1e-300)),
                     color=color, lw=3.4, ls="--", alpha=0.30)
        ax0.semilogx(f, 10 * np.log10(np.maximum(m.power(w), 1e-300)),
                     color=color, lw=1.5, label=f"{fc:g} Hz")
    ax0.set_xlim(10.0, fs / 2.0)
    ax0.set_ylim(-80.0, 3.0)
    ax0.set_xlabel("频率 / Hz")
    ax0.set_ylabel("幅度 / dB")
    ax0.set_title("幅度响应（实线：最佳二阶节；虚线：模拟原型）")
    ax0.grid(True, which="both", alpha=0.25)
    ax0.legend(loc="lower right", fontsize=8)

    for fc, m, _, _, _ in rows:
        err, _, w = error_db(m, fc, fs, 200001)
        ax1.semilogy(w / PI, np.maximum(np.abs(err), 1e-12), lw=1.4,
                     label=f"{fc:g} Hz")
    ax1.set_xlim(0.0, 1.0)
    ax1.set_ylim(1e-8, 1.0)
    ax1.set_xlabel("f / (fs/2)")
    ax1.set_ylabel("|ΔdB|（数字 − 模拟）")
    ax1.set_title("最佳二阶节的全带误差幅值（对数）")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    path1 = out_dir / "highpass_1st_error.png"
    fig.savefig(path1, dpi=150)
    plt.close(fig)

    # 图 2：最佳极点随 fc 的轨迹
    fcs = np.geomspace(20.0, 23900.0, 15)
    lp, extra, bil = [], [], []
    for fc in fcs:
        poles = design_highpass_1st(fc, fs).poles()
        r1 = PM.bilinear_pole(fc, fs)
        i = int(np.argmin(np.abs(poles - r1)))
        lp.append(poles[i])
        extra.append(poles[1 - i])
        bil.append(r1)
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.semilogx(fcs, bil, color="C0", lw=1.2, ls="--",
                label="双线性极点  z = −a1_b")
    ax.semilogx(fcs, lp, "o-", color="C3", lw=1.2, ms=4,
                label="第 1 极点（离双线性极点最近）")
    ax.semilogx(fcs, extra, "s-", color="C1", lw=1.2, ms=4,
                label="第 2 极点（另一个根）")
    ax.axhline(-0.556, color="C1", lw=0.8, ls=":")
    ax.set_xlabel("fc / Hz")
    ax.set_ylabel("极点 z（分母根 p1, p2）")
    ax.set_ylim(-1.02, 1.02)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="center left", fontsize=8)
    ax.set_title("最佳极点随 fc 的轨迹")
    ax.text(0.02, 0.03,
            "极点 z 是 1 + a1 z^-1 + a2 z^-2 = 0 的根：\n"
            "a1 = -(z1 + z2),  a2 = z1 z2",
            transform=ax.transAxes, fontsize=7.5, va="bottom",
            bbox={"boxstyle": "round,pad=0.3", "fc": "white", "alpha": 0.85})
    fig.tight_layout()
    path2 = out_dir / "highpass_1st_poles.png"
    fig.savefig(path2, dpi=150)
    plt.close(fig)
    return path1, path2


def main():
    rows = _report(FS)
    out_dir = Path(__file__).parent / "output"
    p1, p2 = _figures(rows, FS, out_dir)
    print(f"\n图：{p1}\n图：{p2}")


if __name__ == "__main__":
    main()
