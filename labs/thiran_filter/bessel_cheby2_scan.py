# -*- coding: utf-8 -*-
"""
bessel_cheby2_scan.py
=====================

用 scipy 设计两类数字低通滤波器，扫描设计参数、只看群延迟：

- 数字贝塞尔（scipy.signal.bessel，双线性变换；Wn 为归一化临界频率，1 = Nyquist）
- 数字逆切比雪夫 / cheby2（scipy.signal.cheby2，rs = 阻带衰减，Wn = 阻带边沿）

四张图（输出在本脚本同目录的 output/ 下），每张都是「定阶数扫频率」与「定频率扫阶数」之一：

    bessel_scan_cutoff.png   贝塞尔：N 固定，扫 Wn
    bessel_scan_order.png    贝塞尔：Wn 固定，扫 N
    cheby2_scan_cutoff.png   逆切比雪夫：N、rs 固定，扫 Wn
    cheby2_scan_order.png    逆切比雪夫：Wn、rs 固定，扫 N

群延迟由零极点直接求（-d arg H/dw），并在自检里与 scipy.signal.group_delay 逐点比对，
避免高阶时多项式系数的病态影响。
"""

import math
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal as ss

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei"]
matplotlib.rcParams["axes.unicode_minus"] = False

OUT_DIR = Path(__file__).parent / "output"
FS = 2.0  # 数字设计采样率：Wn 归一化到 1 = Nyquist
W = np.linspace(0.0, np.pi, 3001)
WX = W / np.pi


# ------------------------------------------------------------
# 频域量
# ------------------------------------------------------------

def zpk_group_delay(z, p, k, w):
    """由零极点直接算群延迟：tau_g = -d arg H/dw（k 的相位与 w 无关，不参与）。"""
    e = np.exp(1j * np.asarray(w))[:, None]
    z = np.atleast_1d(z).astype(complex)
    p = np.atleast_1d(p).astype(complex)
    tau = np.zeros_like(np.asarray(w), dtype=float)
    for roots, sign in ((p, 1.0), (z, -1.0)):
        if roots.size:
            tau += sign * np.imag(1j * e / (e - roots[None, :])).sum(axis=1)
    return tau


def zpk_response(z, p, k, w):
    e = np.exp(1j * np.asarray(w))[:, None]
    z = np.atleast_1d(z).astype(complex)
    p = np.atleast_1d(p).astype(complex)
    num = np.prod(e - z[None, :], axis=1) if z.size else np.ones_like(e[:, 0])
    den = np.prod(e - p[None, :], axis=1) if p.size else np.ones_like(e[:, 0])
    return k * num / den


def zpk_mag_db(z, p, k, w):
    return 20.0 * np.log10(np.abs(zpk_response(z, p, k, w)))


def bessel_zpk(n, wn, norm="phase"):
    return ss.bessel(n, wn, "low", analog=False, fs=FS, norm=norm, output="zpk")


def bessel_analog_zpk(n, wn, norm="delay"):
    """与数字版等价的模拟原型：scipy 对数字设计做了预畸变 Wn -> 2*fs*tan(pi*Wn/fs)。"""
    return ss.bessel(n, 2.0 * FS * math.tan(math.pi * wn / FS), "low", analog=True, norm=norm, output="zpk")


def analog_group_delay(z, p, om):
    """模拟域群延迟（秒）：tau_a = Re[sum_p 1/(jw-p)] - Re[sum_z 1/(jw-z)]。"""
    jom = 1j * np.asarray(om)[:, None]
    tau = np.zeros_like(np.asarray(om), dtype=float)
    for roots, sign in ((p, 1.0), (z, -1.0)):
        roots = np.atleast_1d(roots).astype(complex)
        if roots.size:
            tau += sign * np.real(1.0 / (jom - roots[None, :])).sum(axis=1)
    return tau


def bilinear_group_delay(z, p, wn, w):
    """数字群延迟的解析预测：tau_d(w) = fs*sec^2(w/2)*tau_a(Omega), Omega = 2*fs*tan(w/2)。"""
    za, pa, _ = bessel_analog_zpk(len(np.atleast_1d(p)), wn)
    om = 2.0 * FS * np.tan(np.asarray(w) / 2.0)
    return FS / np.cos(np.asarray(w) / 2.0) ** 2 * analog_group_delay(za, pa, om)


def cheby2_zpk(n, rs, wn):
    return ss.cheby2(n, rs, wn, "low", analog=False, fs=FS, output="zpk")


# ------------------------------------------------------------
# 出图：两张曲线 + 一张"直流群延迟 vs 扫描参数"
# ------------------------------------------------------------

def draw(cases, title, xlabel_param, fname, clip=None, extra=None, note=None):
    """cases: [(标签, z, p, k, 参数)]；左图群延迟族，右图直流群延迟 vs 扫描参数。
    extra: 右图附加曲线 [(标签, [参数...], [tau0...])]，用于对比另一种归一化。
    note: 左图角上的说明文字。"""
    fig, (ax_g, ax_0) = plt.subplots(1, 2, figsize=(12.0, 4.6), width_ratios=[2.0, 1.0])
    cmap = plt.get_cmap("viridis")
    colors = cmap(np.linspace(0.05, 0.9, len(cases)))

    taus = []
    for (label, z, p, k, param), color in zip(cases, colors):
        tau = zpk_group_delay(z, p, k, W)
        taus.append(tau[0])
        ax_g.plot(WX, tau, color=color, lw=1.6, label=f"{xlabel_param}={param:g}")
        ax_0.plot([param], [tau[0]], "o", color=color, ms=5)

    params = [c[4] for c in cases]
    ax_0.plot(params, taus, color="0.4", lw=1.0, zorder=0, label="本图（见左）")
    for label, xs, ys in (extra or []):
        ax_0.plot(xs, ys, "s--", color="0.55", ms=4, lw=1.0, label=label)

    ax_g.set_title(f"{title}：群延迟")
    ax_g.set_xlabel("ω / π")
    ax_g.set_ylabel("群延迟 τ$_g$ / 样本")
    ax_g.grid(alpha=0.3)
    ax_g.legend(fontsize=8, ncol=2)
    if clip is not None:
        ax_g.set_ylim(*clip)
    if note:
        ax_g.text(0.02, 0.06, note, transform=ax_g.transAxes, fontsize=7.5, color="0.35",
                  bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none", "pad": 1.5})

    ax_0.set_title("直流群延迟 vs 扫描参数")
    ax_0.set_xlabel(xlabel_param)
    ax_0.set_ylabel("τ$_g$(0) / 样本")
    ax_0.grid(alpha=0.3)
    ax_0.set_xscale("log")
    ax_0.set_yscale("log")
    ax_0.set_xticks(params)
    ax_0.set_xticklabels([f"{p:g}" for p in params])
    if extra:
        ax_0.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(OUT_DIR / fname, dpi=110, bbox_inches="tight")
    plt.close(fig)


def bessel_scan_cutoff(n=4, wn_list=(0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8)):
    cases = []
    for wn in wn_list:
        z, p, k = bessel_zpk(n, wn, norm="delay")
        cases.append((f"N={n} Wn={wn:g}", z, p, k, wn))
    phase = [zpk_group_delay(*bessel_zpk(n, wn, norm="phase"), w=[0.0])[0] for wn in wn_list]
    extra = [("scipy 默认 norm='phase'", list(wn_list), phase)]
    draw(cases, f"数字贝塞尔（scipy.signal.bessel，norm='delay'，N={n} 固定，扫 Wn）",
         "Wn", "bessel_scan_cutoff.png", extra=extra,
         note="τ$_g$(0) = 1/(2tan(πWn/2))，与阶数无关；奈奎斯特附近的抬升来自双线性频率拉伸")
    return wn_list


def bessel_scan_order(wn=0.2, orders=(1, 2, 3, 4, 5, 6, 8, 10)):
    cases = []
    for n in orders:
        z, p, k = bessel_zpk(n, wn, norm="delay")
        cases.append((f"n={n}", z, p, k, n))
    phase = [zpk_group_delay(*bessel_zpk(n, wn, norm="phase"), w=[0.0])[0] for n in orders]
    extra = [("scipy 默认 norm='phase'", list(orders), phase)]
    draw(cases, f"数字贝塞尔（norm='delay'，Wn={wn} 固定，扫阶数）",
         "阶数 N", "bessel_scan_order.png", extra=extra,
         note="阶数越高，直流延迟不变（1.5388）、带宽越宽，但奈奎斯特附近的延迟峰越高")
    return orders


def cheby2_scan_cutoff(n=6, rs=40.0, wn_list=(0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6)):
    cases = []
    for wn in wn_list:
        z, p, k = cheby2_zpk(n, rs, wn)
        cases.append((f"N={n} Wn={wn:g}", z, p, k, wn))
    draw(cases, f"数字逆切比雪夫（cheby2，N={n}、rs={rs:g} dB 固定，扫阻带边沿 Wn）",
         "Wn", "cheby2_scan_cutoff.png", clip=(-4.0, 40.0),
         note="传输零点处群延迟 ±∞（垂直线，被坐标轴截断）；有限最大峰在通带边沿附近 ≈0.7·Wn")
    return wn_list


def cheby2_scan_order(wn=0.3, rs=40.0, orders=(2, 3, 4, 5, 6, 8, 10)):
    cases = []
    for n in orders:
        z, p, k = cheby2_zpk(n, rs, wn)
        cases.append((f"n={n}", z, p, k, n))
    draw(cases, f"数字逆切比雪夫（Wn={wn}、rs={rs:g} dB 固定，扫阶数）",
         "阶数 N", "cheby2_scan_order.png", clip=(-4.0, 40.0),
         note="阶数越高：直流延迟越低、尖峰越靠向阻带边沿、有限峰越高；零点数 = N/2（N 偶数）")
    return orders


# ------------------------------------------------------------
# 自检
# ------------------------------------------------------------

def report():
    print("--- 群延迟公式自检 ---")
    ok_band = W <= 0.9 * np.pi  # 近 Nyquist 处系数式求值会因抵消失真（scipy 会发 singularity 警告）
    for name, (z, p, k) in (("bessel n=4 Wn=0.5", bessel_zpk(4, 0.5)),
                            ("cheby2 n=4 rs=40 Wn=0.4", cheby2_zpk(4, 40.0, 0.4))):
        b, a = ss.zpk2tf(z, p, k)
        ref = ss.group_delay((b, a), w=W)[1]
        diff = np.abs(zpk_group_delay(z, p, k, W) - ref)
        print(f"  {name}: vs scipy.signal.group_delay max|diff| = {diff[ok_band].max():.3e}（ω≤0.9π）"
              f"，全带 {diff.max():.3e}（ω→π 处系数式失真）")

    # 与相位数值微分比：不依赖解析导数的独立参照，阶数高时也有效
    for name, (z, p, k) in (("bessel n=10 Wn=0.2", bessel_zpk(10, 0.2)),
                            ("cheby2 n=10 rs=60 Wn=0.2", cheby2_zpk(10, 60.0, 0.2))):
        wc = np.linspace(1e-9, 0.5 * np.pi, 2001)
        h = 1e-6
        fd = -(np.unwrap(np.angle(zpk_response(z, p, k, wc + h)))
               - np.unwrap(np.angle(zpk_response(z, p, k, wc - h)))) / (2.0 * h)
        print(f"  {name}: vs 相位中心差分 max|diff| = "
              f"{np.max(np.abs(zpk_group_delay(z, p, k, wc) - fd)):.3e}（差分步长 {h:g}）")

    print("--- 贝塞尔：数字群延迟就是模拟原型经双线性频率拉伸 ---")
    for n, wn in ((4, 0.2), (10, 0.2)):
        wc = np.linspace(1e-9, 0.999 * np.pi, 2001)
        z, p, k = bessel_zpk(n, wn, norm="delay")
        got = zpk_group_delay(z, p, k, wc)
        pred = bilinear_group_delay(z, p, wn, wc)
        print(f"  n={n:2d} Wn={wn:g}: max|tau_d - fs*sec²(w/2)*tau_a(Ω)| = {np.max(np.abs(got - pred)):.3e}"
              f"，数字峰值 {got.max():.3f} 样本（模拟原型 τ_a(0) 仅 {analog_group_delay(*bessel_analog_zpk(n, wn)[:2], [1e-9])[0]:.4f} 秒）")

    print("--- 贝塞尔：两种归一化下的直流群延迟（Wn=0.2）---")
    for n in (2, 4, 6, 10):
        z, p, k = bessel_zpk(n, 0.2, norm="delay")
        tau_delay = zpk_group_delay(z, p, k, W)[0]
        z, p, k = bessel_zpk(n, 0.2, norm="phase")
        tau_phase = zpk_group_delay(z, p, k, W)[0]
        print(f"  n={n:2d}: norm='delay' τ_g(0)={tau_delay:.4f}（=1/(2tan(πWn/2))={1/(2*math.tan(math.pi*0.1)):.4f}）"
              f"   norm='phase' τ_g(0)={tau_phase:.4f}  max|p|={np.max(np.abs(p)):.4f}")

    print("--- 逆切比雪夫：阻带波纹与通带 ---")
    for n, rs, wn in ((4, 40.0, 0.3), (6, 40.0, 0.3), (10, 60.0, 0.2)):
        z, p, k = cheby2_zpk(n, rs, wn)
        mag = zpk_mag_db(z, p, k, W)
        tau = zpk_group_delay(z, p, k, W)
        stop = WX >= wn
        pb = WX < wn
        # 传输零点（单位圆上的零点对应的频率）
        zeros_angle = np.sort(np.angle(z[np.abs(np.abs(z) - 1.0) < 1e-9]) / np.pi)
        print(f"  n={n:2d} rs={rs:g} Wn={wn:g}: 阻带波纹峰值={mag[stop].max():+.4f} dB（应=-rs）, "
              f"阻带最深={mag[stop].min():+.4f} dB（传输零点）, 通带最大={mag[pb].max():+.4f} dB, "
              f"Wn 前一点={mag[np.searchsorted(WX, wn) - 1]:+.4f} dB, "
              f"直流 τ={tau[0]:.4f}, max|p|={np.max(np.abs(p)):.4f}")
        print(f"      传输零点 ω/π = {np.round(zeros_angle[zeros_angle >= 0], 4)}")


if __name__ == "__main__":
    OUT_DIR.mkdir(exist_ok=True)
    report()
    bessel_scan_cutoff()
    bessel_scan_order()
    cheby2_scan_cutoff()
    cheby2_scan_order()
    print(f"\n图已写入 {OUT_DIR}")
