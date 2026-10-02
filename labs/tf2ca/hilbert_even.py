# -*- coding: utf-8 -*-
"""
hilbert_even.py
==============

`qwqdsp/include/qwqdsp/filter/iir_hilbert.hpp` 是怎么**构造**出来的？——本脚本给出可复现的构造链：

    偶数阶椭圆半带低通（crossover 落在 0.5，极点正好落在 z 平面的虚轴上）
        -> 复系数全通链分解  H = alpha*A1 + conj(alpha)*A2      （alpha 落在 45 度）
        -> 对两条链做频率旋转（极点乘 j，等价 z -> -jz），系数变实数
        -> 解析滤波器：一条实系数多相全通链 + 一条带 1 拍延迟的实系数多相全通链

最后一步的结构与 ``iir_hilbert.hpp`` 里的 ``APF`` / ``Tick`` 完全一致：
``APF`` 有两个状态，传递函数是 ``(z^-2 - a)/(1 - a*z^-2)``（zeta = z^-2 的多相全通），
``a`` 就是半带滤波器的半速率极点 ``rho^2``；``imag`` 支路多走一拍（``latch_``）。

为什么要这样：
1. **要让数字极点落在虚轴上，椭圆滤波器必须同时满足三条**（推导与验证见 ``hilbert_even.md``）：
   阶数偶数、波纹互补 ``eps_p*eps_s = 1``（于是 ``rs`` 由 ``rp`` 定死）、边沿对称
   ``Wp + Ws = 1``（即 prewarp 后 ``tan(pi*Wp/2)*tan(pi*Ws/2) = 1``，crossover 落在 0.5）。
   第三条有**闭式解**（椭圆次数方程的 nome 形式 ``q = q1^(1/N)`` + theta 级数，
   见 :func:`halfband_Wp`），所以可以像要求的那样——用 scipy 的数字 ``ellip`` 直接出滤波器。
   三条不满足时（例如只把 ``Wn`` 设成 0.5、波纹随意），极点会明显偏离虚轴。
2. 这种极点分布下，**实数**两路全通做不到（奇数阶那套 ``1/2(A0+A1)`` 在这里退化），
   但**复系数**两路全通可以，且权重必为 ``alpha``、``conj(alpha)``，实测 ``alpha`` 恰为 45 度
   （`(1+j)/(2*sqrt(2))`），这正是两条链相差 90 度、从而能拼成解析滤波器的原因。
3. 频率旋转把虚轴极点 ``±j*rho`` 变成实轴极点 ``∓rho``，系数变实；每组 ``±rho`` 合成一个
   ``zeta`` 多相节（``a = rho^2``），两支路交替分配 —— 就是头文件里那两组 alpha。

用法
----
    python hilbert_even.py                          # 半带 rs=30dB（默认），出图
    python hilbert_even.py --order 16 --rs 40
    python hilbert_even.py --header                 # 直接分析 iir_hilbert.hpp 里的固定系数
"""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 无界面后端
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal as sg
from scipy.special import ellipk, ellipkm1

SCRIPT_DIR = Path(__file__).resolve().parent


# ------------------------------------------------------------
# 头文件里的两组固定系数（iir_hilbert.hpp / iir_cpx_hilbert.hpp）
# ------------------------------------------------------------

HEADER_REAL = [0.4021921162426, 0.8561710882420, 0.9722909545651, 0.9952884791278]
HEADER_IMAG = [0.6923878, 0.9360654322959, 0.9882295226860, 0.9987488452737]
HEADER_DEEPER_REAL = [0.0406273391966415, 0.2984386654059753, 0.5938455547890998,
                      0.7953345677003365, 0.9040699927853059, 0.9568366727621767,
                      0.9815966237057977, 0.9938718801312583]
HEADER_DEEPER_IMAG = [0.1500685240941415, 0.4538477444783975, 0.7081016258869689,
                      0.8589957406397113, 0.9353623391637175, 0.9715130669899118,
                      0.9886689766148302, 0.9980623781456869]


# ------------------------------------------------------------
# 1. 精确半带椭圆滤波器（crossover 在 0.5）
# ------------------------------------------------------------


def k_from_q(q: float) -> float:
    """椭圆模数 ``k = theta2^2/theta3^2``（q 的级数，q 不接近 1 时几项就够）。"""
    n = 1
    while q ** (n * n) > 1e-18 and n < 500:
        n += 1
    m = np.arange(0, n + 1)
    s1 = np.sum(q ** (m * (m + 1)))              # theta2 的求和
    s2 = np.sum(q ** (np.arange(1, n + 1) ** 2))  # theta3 的求和
    return float(4 * np.sqrt(q) * s1**2 / (1 + 2 * s2) ** 2)


def ripple_rp(rs_db: float) -> float:
    """互补波纹（``eps_p*eps_s = 1``）下由**阻带深度**反推通带波纹。

    ``eps_s = sqrt(10^(rs/10)-1)``、``eps_p = 1/eps_s`` ⇒ ``rp = 10*log10(1 + eps_p^2)``。
    """
    es = np.sqrt(10 ** (rs_db / 10) - 1)
    return float(10 * np.log10(1 + 1 / es**2))


def halfband_Wp(order: int, rs_db: float):
    """闭式解：让数字极点在虚轴上的**通带边沿** ``Wp``（crossover 落在 0.5）。

    参数取**阻带深度** ``rs``（dB）：半带条件 ``eps_p*eps_s = 1`` 把通带波纹一起定死，
    所以自由参数只有 ``(order, rs)``。

    闭式来自椭圆次数方程的 nome 形式：``q = q1^(1/N)``、``q1 = exp(-pi*K'(k1)/K(k1))``、
    ``k1 = eps_p^2 = 1/eps_s^2``，再由 theta 级数得 ``k = theta2^2/theta3^2``，
    最后 ``Wp = (2/pi)*atan(sqrt(k))``（``Ws = 1 - Wp``）。

    :return: ``(Wp, k)``
    """
    es = np.sqrt(10 ** (rs_db / 10) - 1)
    ep = 1.0 / es
    k1 = ep**2                      # = eps_p/eps_s
    m1 = k1**2
    q1 = np.exp(-np.pi * ellipkm1(m1) / ellipk(m1))   # ellipkm1(m) = K(1-m)，小 m 时不退化
    q = q1 ** (1.0 / order)
    k = k_from_q(q)
    return float((2 / np.pi) * np.arctan(np.sqrt(k))), k


def halfband_design(order: int, rs_db: float) -> dict:
    """偶数阶椭圆**半带**低通：scipy 数字 ``ellip`` 直出，边沿用闭式解，极点取 ``output='zpk'``。

    自由参数只有 ``order``（偶数）与 ``rs``（阻带深度 dB）。

    半速率极点 ``rho^2`` 直接取 scipy 给的**数字极点**（``p = ±j*rho`` ⇒ ``rho^2 = |p|^2``），
    比 ``np.roots(A)`` 稳定得多：高阶时后者的极点分离条件数会崩。

    :return: ``{"b","a","rp","rs","Wp","k","rho2"}``
    """
    rp = ripple_rp(rs_db)
    Wp, k = halfband_Wp(order, rs_db)
    if not (0.0 < Wp < 0.5):
        raise ValueError(f"闭式 Wp={Wp} 越界（rs={rs_db} dB / N={order} 超出双精度可表示范围）")
    b, a = sg.ellip(order, rp, rs_db, Wp)
    _, poles, _ = sg.ellip(order, rp, rs_db, Wp, output="zpk")
    rho2 = np.sort(np.abs(poles[poles.imag > 0]) ** 2)
    return {"b": b, "a": a, "rp": rp, "rs": float(rs_db), "Wp": Wp, "k": k,
            "rho2": rho2, "poles": np.asarray(poles)}


def halfband_report(d: dict) -> dict:
    """核对半带的关键不变量：``A(z)=A_e(z^-2)``、极点在虚轴、crossover、功率互补。"""
    b, a = d["b"], d["a"]
    w = np.linspace(0, np.pi, 40001)
    _, H = sg.freqz(b, a, worN=w)
    _, Hm = sg.freqz(b, a, worN=np.pi - w)
    poles_roots = np.roots(a)
    rho2_roots = np.sort(-1 / np.roots(a[0::2][::-1]))    # 另一条路：zeta 多项式求根
    return {
        "a_odd": float(np.max(np.abs(a[1::2]))),          # 0 => A(z)=A_e(z^-2)
        "pole_re": float(np.max(np.abs(poles_roots.real))),  # 由多项式求根（高阶不可信）
        "pole_re_zpk": float(np.max(np.abs(d["poles"].real))),  # 设计器直给的极点，可信
        "h_half": float(np.abs(H[20000])),                # |H(0.5)| = 1/sqrt(2)
        "power_comp": float(np.max(np.abs(np.abs(H) ** 2 + np.abs(Hm) ** 2 - 1))),
        "rho2_zpk_vs_roots": float(np.max(np.abs(np.asarray(d["rho2"]) - rho2_roots))),
        "freq": w,
    }


# ------------------------------------------------------------
# 2. 复系数两路全通分解
# ------------------------------------------------------------


def conj_rev(d: np.ndarray) -> np.ndarray:
    """复全通分子：``D^# = conj(降幂系数反序)``。"""
    return np.conj(np.asarray(d)[::-1])


def complex_split(b: np.ndarray, a: np.ndarray, poles=None):
    """穷举共轭对取法，求 ``H = alpha*A1 + conj(alpha)*A2`` 的最佳解（**偶阶**专用）。

    偶阶半带的极点是纯虚数 ``±j*rho``（成对，无实极点），所以"每个共轭对取一个"就是"每个 ± 对各取一个"。
    奇阶另有一个自共轭的原点极点，无法这样拆 —— 见 ``hilbert_odd.py``（奇阶走"两组 AP 链 + 延迟"）。

    :return: ``(bits, d1, d2, alpha, residual, (sel1, sel2))``
    """
    tol = 1e-9
    poles = np.asarray(poles) if poles is not None else np.roots(a)
    upper = sorted((pp for pp in poles if pp.imag > tol), key=abs)
    scale = float(np.max(np.abs(b)))
    best = None
    for bits in range(1 << len(upper)):
        s1 = [pp if (bits >> k) & 1 else np.conj(pp) for k, pp in enumerate(upper)]
        s2 = [np.conj(pp) for pp in s1]
        d1, d2 = np.poly(s1), np.poly(s2)
        p = np.convolve(conj_rev(d1), d2)
        u, v = p + np.conj(p), 1j * (p - np.conj(p))
        xy, *_ = np.linalg.lstsq(np.stack([u.real, v.real], axis=1), b, rcond=None)
        alpha = xy[0] + 1j * xy[1]
        recon = alpha * p + np.conj(alpha) * np.conj(p)
        res = float(np.max(np.abs(recon - b)) / scale)
        if best is None or res < best[0]:
            best = (res, bits, d1, d2, alpha, np.array(s1), np.array(s2))
    res, bits, d1, d2, alpha, s1, s2 = best
    return bits, d1, d2, alpha, res, (s1, s2)


def rotate_branch(d: np.ndarray) -> np.ndarray:
    """频率旋转 ``z -> -jz``：系数 ``d_k -> d_k * j^k``（极点在 z 平面旋转 90 度）。"""
    d = np.asarray(d, dtype=complex)
    return d * (1j ** np.arange(d.size))


# ------------------------------------------------------------
# 3. 头文件式的解析滤波器（实系数多相全通 + 延迟）
# ------------------------------------------------------------


def analytic_response(real_alphas, imag_alphas, w):
    """解析滤波器 ``H = (A_real(z) + j * z^-1 * A_imag(z)) / 2``（``iir_hilbert.hpp`` 的结构除以 2）。

    与 ``iir_hilbert.hpp`` 的 ``Tick`` 逐行对应（``latch_`` 即 imag 支路的那一拍延迟）；
    **注意头文件本身不除 2**，所以 ``Tick`` 的原始输出是这里的 2 倍（教材里解析滤波器取
    ``H = 1/2[LP + j*HP]``，正频增益 1）。归一化后：正频 ``|H| = 1``、负频 ``|H| = 0``。
    """
    z2 = np.exp(-2j * np.asarray(w, dtype=float))
    real = np.ones_like(z2)
    imag = np.ones_like(z2)
    for aa in real_alphas:
        real *= (z2 - aa) / (1 - aa * z2)
    for aa in imag_alphas:
        imag *= (z2 - aa) / (1 - aa * z2)
    return (real + 1j * np.exp(-1j * np.asarray(w, dtype=float)) * imag) / 2.0


def _zfreq(poly, w):
    """``Sum_k poly[k] * z^-k`` 在 ``z = exp(jw)`` 上取值（poly 为降幂系数）。"""
    return np.polyval(np.asarray(poly)[::-1], np.exp(-1j * np.asarray(w, dtype=float)))


def chain_response_poles(poles, w):
    """同一全通链，但按**一阶节连乘**求值：``prod (z^-1 - conj(p)) / (1 - p z^-1)``。

    与按多项式求值 ``D^#/D`` 等价（``D`` 首项为 1 时逐节恒等），但高阶时不经过高次
    多项式求值，数值稳定得多 —— 高阶（N>=32）必须走这条。
    """
    zinv = np.exp(-1j * np.asarray(w, dtype=float))
    h = np.ones_like(zinv)
    for pp in poles:
        h *= (zinv - np.conj(pp)) / (1 - pp * zinv)
    return h


def band_stats(h, w, guard=0.05):
    """通带/阻带统计。``guard`` 是排除 DC / Nyquist 附近过渡区的比例。

    负频最大是**在负频网格上直接取 ``max|h|``**（不是寻峰），所以可能漏掉极窄的峰；
    同时报告取得最大值的频率，便于判断它落在保护带边缘还是阻带内部。
    """
    band = (np.abs(w) > guard * np.pi) & (np.abs(w) < (1 - guard) * np.pi)
    pos, neg = band & (w > 0), band & (w < 0)
    k = int(np.argmax(np.abs(h[neg]))) if np.any(neg) else 0
    return {
        "pos_min": float(np.min(np.abs(h[pos]))),
        "pos_max": float(np.max(np.abs(h[pos]))),
        "neg_max": float(np.max(np.abs(h[neg]))),
        "neg_argmax_w": float(w[neg][k]) if np.any(neg) else 0.0,
        "guard": guard,
    }


def analytic_metrics(real_alphas, imag_alphas, w=None, guard=0.05) -> dict:
    """解析滤波器(已除 2)的质量：通带(正频)幅度偏差、负频抑制、90 度相位偏差。

    只统计 ``|w| in [0.05pi, 0.95pi]``（排除 DC / Nyquist 附近的过渡区 —— 那里的取值由
    半带的 crossover 决定，恒为 ``sqrt(2)``）。
    """
    w = np.linspace(-np.pi, np.pi, 16001) if w is None else w
    h = analytic_response(real_alphas, imag_alphas, w)
    st = band_stats(h, w, guard)
    pos = (np.abs(w) > guard * np.pi) & (np.abs(w) < (1 - guard) * np.pi) & (w > 0)
    z2 = np.exp(-2j * w)
    ar = np.ones_like(z2)
    for aa in real_alphas:
        ar *= (z2 - aa) / (1 - aa * z2)
    ai = np.ones_like(z2)
    for aa in imag_alphas:
        ai *= (z2 - aa) / (1 - aa * z2)
    dphi = np.unwrap(np.angle(ar) - np.angle(np.exp(-1j * w) * ai))
    dev = np.abs(np.abs(dphi[pos]) - np.pi / 2)
    return {
        "mag_band": (st["pos_min"], st["pos_max"]),
        "neg_max": st["neg_max"],
        "neg_argmax_w": st["neg_argmax_w"],
        "phase_dev_deg": float(np.max(dev) * 180 / np.pi),
    }


def plot_header_branch_phase(out: Path) -> None:
    """头文件里两组固定系数的 branch phase（含结构延迟、不含权重 j）。

    ``iir_hilbert.hpp`` 的结构是 ``A_real(z) + j*z^-1*A_imag(z)``，每个 ``APF`` 为
    ``(z^-2 - a)/(1 - a*z^-2)``。两条链（含延迟）的相位差应为 **-90 度（正频）/ +90 度（负频）**；
    在 DC 与 ±pi 处曲线连续穿过 ±180 分支割线，绘图在这几个样本处断开。
    """
    w = np.linspace(-np.pi, np.pi, 16001)
    zeta = np.exp(-2j * w)

    def chain(alphas):
        h = np.ones_like(zeta)
        for aa in alphas:
            h *= (zeta - aa) / (1 - aa * zeta)
        return h

    fig, ax = plt.subplots(figsize=(9, 5))
    for name, (ra, ia), style in (("IIRHilbertDeeper (8+8)", (HEADER_DEEPER_REAL, HEADER_DEEPER_IMAG), "-"),
                                  ("IIRHilbert (4+4)", (HEADER_REAL, HEADER_IMAG), "--")):
        d = np.degrees(np.angle(np.exp(-1j * w) * chain(ia) / chain(ra)))
        d[np.abs(w) < 1e-12] = np.nan
        d[np.abs(np.abs(w) - np.pi) < 1e-12] = np.nan
        band = (np.abs(w) > 0.05 * np.pi) & (np.abs(w) < 0.95 * np.pi)
        dev = np.max(np.abs(np.abs(d[band]) - 90.0))
        ax.plot(w / np.pi, d, style, lw=2, label=f"{name}: max dev = {dev:.4f} deg")
    ax.axhline(-90, color="k", ls=":", label="-90 deg (w > 0)")
    ax.axhline(90, color="gray", ls=":", label="+90 deg (w < 0)")
    ax.set(title="iir_hilbert.hpp: branch phase difference (incl. the structural delay)",
           xlabel="$\\omega/\\pi$", ylabel="deg", ylim=(-200, 200))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="center left")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"图已保存：{out}")


# ------------------------------------------------------------
# 绘图
# ------------------------------------------------------------


def plot(b, a, info, d, d1, alpha, p1r, p2r, out: Path) -> None:
    """六张图：半带极点分组 / 半带原始幅度 / 旋转后实极点 / 解析幅度(与目标对照) / 相差 / 结构对照。"""
    fig, ax = plt.subplots(2, 3, figsize=(17.5, 9.5))
    w = np.linspace(-np.pi, np.pi, 16001)

    # (0,0) 半带 z 平面：极点在虚轴上，按复全通分组着色
    th = np.linspace(0, 2 * np.pi, 400)
    ax[0, 0].plot(np.cos(th), np.sin(th), "k:", lw=0.8)
    p1, p2 = np.roots(d1), np.roots(np.conj(d1))
    ax[0, 0].plot(p1.real, p1.imag, "x", color="#1f77b4", ms=11, label="chain 1 (D1)")
    ax[0, 0].plot(p2.real, p2.imag, "x", color="#d62728", ms=11, label="chain 2 (conj D1)")
    for pp in p1:
        ax[0, 0].plot([pp.real, pp.real], [pp.imag, -pp.imag], color="gray", lw=0.7, alpha=0.7)
    ax[0, 0].plot(np.roots(b).real, np.roots(b).imag, "o", mfc="none", mec="g", ms=8,
                  label="B zeros")
    ax[0, 0].axhline(0, color="gray", lw=0.5); ax[0, 0].axvline(0, color="gray", lw=0.5)
    ax[0, 0].set(title="halfband poles: on the imaginary axis, split per chain",
                 xlabel="Re", ylabel="Im", aspect="equal")
    ax[0, 0].grid(alpha=0.3); ax[0, 0].legend(loc="upper left", fontsize=8)

    # (0,1) 半带自身的原始幅度响应（dB）—— 确认它确实是半带
    ww = np.linspace(1e-9, np.pi - 1e-9, 20001)
    Hhb = sg.freqz(b, a, worN=ww)[1]
    ax[0, 1].plot(ww / np.pi, 20 * np.log10(np.maximum(np.abs(Hhb), 1e-30)), lw=2)
    ax[0, 1].axhline(-info["rs"], color="r", ls=":", lw=1, label=f"design rs = {info['rs']:g} dB")
    ax[0, 1].axhline(-3.01, color="k", ls=":", lw=1, label="-3.01 dB (crossover)")
    ax[0, 1].axvline(0.5, color="gray", lw=0.6)
    ax[0, 1].axvline(d["Wp"], color="g", ls="--", lw=1, label=f"Wp = {d['Wp']:.6f}")
    ax[0, 1].axvline(1 - d["Wp"], color="g", ls="--", lw=1)
    ax[0, 1].plot(0.5, 20 * np.log10(np.abs(sg.freqz(b, a, worN=[np.pi / 2])[1][0])), "ko", ms=5)
    ax[0, 1].set(title="halfband LP: raw |H| (dB)", xlabel="$\\omega/\\pi$", ylabel="dB",
                 ylim=(-1.35 * info["rs"] - 6, 3))
    ax[0, 1].grid(alpha=0.3, which="both"); ax[0, 1].legend(fontsize=7, loc="lower left")

    # (0,2) 旋转后的实轴极点（精确链：一阶实极点，两链交错）
    for name, pp, col, y in (("chain 1'", p1r, "#1f77b4", 0.0), ("chain 2'", p2r, "#d62728", 0.05)):
        ax[0, 2].plot(np.real(pp), [y] * len(pp), "x", color=col, ms=11, label=name)
    ax[0, 2].plot([-1.05, 1.05], [0, 0], "k:", lw=0.8)
    ax[0, 2].set(title="after rotation: real 1st-order poles (both signs), interlaced",
                 xlabel="Re (z)", ylim=(-0.4, 0.4))
    ax[0, 2].grid(alpha=0.3); ax[0, 2].legend(fontsize=8)

    # (1,0) 解析滤波器 |H|（dB）：精确构造 vs 目标（旋转半带）
    h_exact = alpha * chain_response_poles(p1r, w) + np.conj(alpha) * chain_response_poles(p2r, w)
    h_target = sg.freqz(b, a, worN=w - np.pi / 2)[1]
    ax[1, 0].plot(w / np.pi, 20 * np.log10(np.maximum(np.abs(h_exact), 1e-30)), lw=2,
                  label="alpha*A1' + conj(alpha)*A2'")
    ax[1, 0].plot(w / np.pi, 20 * np.log10(np.maximum(np.abs(h_target), 1e-30)), "k--", lw=1,
                  label="target: halfband rotated by $-\\pi/2$")
    ax[1, 0].axhline(0, color="k", ls=":", lw=1)
    ax[1, 0].set(title="analytic filter (dB): stopband = halfband stopband (= rs)",
                 xlabel="$\\omega/\\pi$", ylabel="dB", ylim=(-1.35 * info["rs"] - 6, 3))
    ax[1, 0].grid(alpha=0.3, which="both"); ax[1, 0].legend(fontsize=8, loc="lower left")

    # (1,1) 两条链的相位差：正频 -90°、负频 +90°（标准解析滤波器的约定：第二支滞后 90°）。
    # 曲线在 DC 处 *连续* 穿过 0，不需要 unwrap；曲线本身是两带内的平直线。
    dphi = np.angle(chain_response_poles(p1r, w) / chain_response_poles(p2r, w))
    ax[1, 1].plot(w / np.pi, np.degrees(dphi), lw=2, label="arg A1' - arg A2'")
    ax[1, 1].axhline(-90, color="k", ls=":", label="-90 deg")
    ax[1, 1].axhline(90, color="gray", ls=":", label="+90 deg")
    ax[1, 1].set(title="branch phase difference: +/-90 deg (sign flips with the band)",
                 xlabel="$\\omega/\\pi$", ylabel="deg", ylim=(-135, 135))
    ax[1, 1].grid(alpha=0.3); ax[1, 1].legend(fontsize=8)

    ax[1, 2].axis("off")   # 第 6 格留空（与 hilbert_odd.py 同版式）

    fig.suptitle(f"Hilbert from halfband ellip(n={len(a)-1}, rs={info['rs']:g} dB) -> "
                 f"complex chains (arg alpha = {np.angle(alpha)*180/np.pi:.2f} deg) -> rotation",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"图已保存：{out}")


# ------------------------------------------------------------
# main
# ------------------------------------------------------------


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="construct iir_hilbert.hpp structure from a halfband")
    p.add_argument("--order", type=int, default=8, help="半带阶数（偶数）")
    p.add_argument("--rs", type=float, default=30.0, help="半带阻带深度 dB（决定 rp；越深则过渡带越宽）")
    p.add_argument("--header", action="store_true", help="分析 iir_hilbert.hpp 里的固定系数")
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args(argv)

    if args.header:
        w = np.linspace(-np.pi, np.pi, 16001)
        for name, (rr, ii) in (("IIRHilbert", (HEADER_REAL, HEADER_IMAG)),
                               ("IIRHilbertDeeper", (HEADER_DEEPER_REAL, HEADER_DEEPER_IMAG))):
            m = analytic_metrics(rr, ii, w)
            print(f"== {name}（{len(rr)}+{len(ii)} 节）==")
            print(f"  alpha 合并排序 = {np.round(np.sort(rr + ii), 6)}")
            print(f"  正频 |H| 范围 = [{m['mag_band'][0]:.4f}, {m['mag_band'][1]:.4f}]"
                  f"（除以 2 后的解析滤波器，理想 1.0；头文件 Tick 原始输出 ×2）")
            print(f"  负频最大 |H| = {m['neg_max']:.3e}"
                  f"（{20*np.log10(max(m['neg_max'], 1e-30)):.2f} dB）")
            print(f"  相位差相对 -90 度最大偏差 = {m['phase_dev_deg']:.3f} 度")
        out = args.out or SCRIPT_DIR / "output" / "hilbert_header_branch_phase.png"
        plot_header_branch_phase(out)
        return

    if args.order % 2:
        p.error("本脚本只处理偶数阶（奇阶见 hilbert_odd.py）")

    d = halfband_design(args.order, args.rs)
    b, a = d["b"], d["a"]
    info = halfband_report(d)
    print(f"== 半带椭圆 n={args.order}, rs={args.rs:g} dB -> rp={d['rp']:.5f} dB, "
          f"Wp={d['Wp']:.9f} (crossover 0.5) ==")
    print(f"  A 的奇次系数 max       = {info['a_odd']:.2e}   （0 => A(z)=A_e(z^-2)）")
    print(f"  极点实部 max|Re p|     = {info['pole_re_zpk']:.2e}（zpk 直给）"
          f" / {info['pole_re']:.2e}（多项式求根，高阶不可信）")
    print(f"  |H(0.5)|               = {info['h_half']:.9f}（1/sqrt(2) = {1/np.sqrt(2):.9f}）")
    print(f"  功率互补误差           = {info['power_comp']:.2e}")
    print(f"  半速率极点 rho^2        = {np.round(d['rho2'], 6)}"
          f"（zpk 与多项式求根之差 {info['rho2_zpk_vs_roots']:.1e}）")
    print()

    bits, d1, d2, alpha, res, (sel1, sel2) = complex_split(b, a, poles=d["poles"])
    print("== 复系数全通分解：H = alpha*A1 + conj(alpha)*A2 ==")
    print(f"  最优分组 bits = {bits:0{len(np.roots(a))//2}b}（|p| 升序，1 = 取上半平面根），残差 = {res:.2e}")
    n_origin = int(np.sum(np.abs(sel1) < 1e-9) + np.sum(np.abs(sel2) < 1e-9))
    conj_ok = (abs(abs(alpha) - 0.5) < 1e-3 and abs(abs(np.angle(alpha)) * 180 / np.pi - 45.0) < 0.5
               and res < 1e-6)
    print(f"  alpha = {alpha:.9g}   |alpha| = {abs(alpha):.6f}, arg = {np.angle(alpha) * 180 / np.pi:.4f} 度")
    print(f"  原点极点(=> z^-1 延迟单元) 个数 = {n_origin}"
          f"（奇阶半带为 1，偶阶为 0）")
    print(f"  共轭形式是否适用 = {conj_ok}（要求 arg alpha = 45 度、|alpha| = 0.5、残差 ~1e-14）")
    if not conj_ok:
        print("  -> 不适用：原点极点是自共轭的，无法按「每个共轭对取一个」拆分。"
              "奇阶半带对应的是\"实数两路 + 延迟\"那种结构（其精确性本目录尚未建立，见 README）。")
    print(f"  => 权重相差 {2 * np.angle(alpha) * 180 / np.pi:.2f} 度（解析滤波器要求 90 度）")

    d1r, d2r = rotate_branch(d1), rotate_branch(d2)
    p1r, p2r = 1j * sel1, 1j * sel2        # 旋转后的极点 = j*(原极点)
    print()
    print("== 频率旋转 z -> -jz（极点 ×j）==")
    print(f"  旋转后 D1 虚部 max = {np.max(np.abs(d1r.imag)):.2e}"
          f"（应为 0 => 实系数）；极点 = {np.round(np.sort_complex(np.roots(d1r)), 6)}")
    print(f"  旋转后 D2 虚部 max = {np.max(np.abs(d2r.imag)):.2e}")

    w = np.linspace(-np.pi, np.pi, 16001)
    guard = max(0.05, 0.5 - d["Wp"])          # 保护带按设计的过渡带宽度自适应
    h_exact = alpha * chain_response_poles(p1r, w) + np.conj(alpha) * chain_response_poles(p2r, w)
    st = band_stats(h_exact, w, guard)
    h_target = sg.freqz(b, a, worN=w - np.pi / 2)[1]          # 目标 = 旋转后的半带
    bd = 20 * np.log10(np.maximum(np.abs(h_exact), 1e-30))
    if not conj_ok:
        print()
        print("== 解析滤波器（共轭形式不适用，下列数字仅供参考）==")
    print()
    print("== 解析滤波器（精确构造：alpha*A1' + conj(alpha)*A2'，即 H_halfband(-jz)）==")
    print(f"  与目标(旋转半带) 最大偏差      = {np.max(np.abs(h_exact - h_target)):.2e}")
    print(f"  正频 |H| 范围 = [{st['pos_min']:.5f}, {st['pos_max']:.5f}]"
          f"（{20*np.log10(st['pos_min']):.4f} … {20*np.log10(st['pos_max']):.4f} dB）")
    print(f"  负频最大 |H|  = {st['neg_max']:.3e} = {20*np.log10(st['neg_max']):.2f} dB"
          f" @ w = {st['neg_argmax_w']/np.pi:+.3f}pi（带内 {st['guard']:.2f}..{1-st['guard']:.2f}pi 网格最大值）")
    ratio = np.angle(chain_response_poles(p1r, w) / chain_response_poles(p2r, w))
    sel_b = (np.abs(w) > guard * np.pi) & (np.abs(w) < (1 - guard) * np.pi)
    dev_pos = np.max(np.abs(ratio[sel_b & (w > 0)] + np.pi / 2)) * 180 / np.pi
    dev_neg = np.max(np.abs(ratio[sel_b & (w < 0)] - np.pi / 2)) * 180 / np.pi
    print(f"  两链相位差（标准约定：第二支滞后）：正频 -90°±{dev_pos:.3f}°，负频 +90°±{dev_neg:.3f}°")

    out = args.out or SCRIPT_DIR / "output" / f"hilbert_n{args.order}_rs{args.rs:g}.png"
    plot(b, a, {"rs": args.rs, "rp": d["rp"]}, d, d1, alpha, p1r, p2r, out)


if __name__ == "__main__":
    main()
