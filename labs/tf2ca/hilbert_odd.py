# -*- coding: utf-8 -*-
"""
hilbert_odd.py
==============

**奇阶**半带滤波器的「两组 AP 链 + 延迟单元」实现 —— 本文件**不做频率旋转**，
目标只有一个：用全通链把半带滤波器本身（幅度与相位）复原出来：

    H(z) = 1/2 * [ A0(zeta) + z^-1 * A1(zeta) ]            zeta = z^-2

- `A0`/`A1` 是 `zeta = z^-2` 上的**实系数全通链**；每个 `±rho` 对合成一个 `zeta` 节，
  节参数就是 `a = rho^2`（与 `iir_hilbert.hpp` 的 `APF` 同形）。
- `z^-1` 是**延迟单元**。奇阶半带的专属结构来源：`A` 的末位系数为 0（`a[-1] ~ 1e-15`），
  有效分母阶是 `N-1`（偶数、只含 `z^-2`），因此结构里必然带一个 `z^-1`。

设计条件与偶阶相同（`hilbert_even.py` / `hilbert_even.md`）：
阶数 `N` 与阻带深度 `rs`（互补波纹 `eps_p*eps_s = 1` 把 `rp` 定死）、
边沿对称 `Wp + Ws = 1`（`Wn` 由闭式解给出）。**本文件只处理奇阶。**

实测（`--check` / 单次运行都会打印）：

| N | rs | 链的 ζ 阶数 | 残差 | 权重 |
|---|---|---|---|---|
| 9 | 30 dB | 2 + 2 | 5.1e-12 | 0.5 / 0.5 |
| 15 | 60 dB | **4 + 3** | 1.3e-09 | 0.5 / 0.5 |
| 17 | 80 dB | 4 + 4 | 6.0e-09 | 0.5 / 0.5 |

注意两条链的阶数**可以不等**（N=15 是 4+3）——只枚举等分是找不到解的，这也是本目录
早期版本误判"该形式不成立"的原因。

用法
----
    python hilbert_odd.py                        # 默认 N=15, rs=60dB，出图
    python hilbert_odd.py --N 9 --rs 30
    python hilbert_odd.py --check                # 扫几组 (N, rs)
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

from hilbert_even import band_stats, halfband_design, halfband_report   # 设计/指标共用（见 hilbert_even.md）

SCRIPT_DIR = Path(__file__).resolve().parent


# ------------------------------------------------------------
# zeta 域全通链
# ------------------------------------------------------------


def zeta_poles(a: np.ndarray) -> np.ndarray:
    """半带有效分母 `A(z) = A_e(z^-2)` 的 zeta 极点（负实数，zeta = z^-2）。"""
    a = np.asarray(a).ravel()
    return np.sort_complex(np.roots(np.array([x for k, x in enumerate(a) if k % 2 == 0])[::-1]))


def zeta_chain_arg(poles, zeta):
    """`zeta` 域实系数全通链 `D^#(zeta)/D(zeta)`（`D = prod(zeta - p)`）在给定自变量上取值。

    显式传入自变量，是为了让"频率旋转 `zeta -> -zeta`"不产生歧义：
    直接算 `A(-zeta)`，而不是把极点取负再去评估（后者会多出 `(-1)^n` 因子）。
    """
    den = np.atleast_1d(np.poly(list(poles)))
    return np.polyval(den[::-1], zeta) / np.polyval(den, zeta)


def zeta_chain(poles, w):
    """`zeta = z^-2` 上的全通链复频响（`z = exp(jw)`）。"""
    return zeta_chain_arg(poles, np.exp(-2j * np.asarray(w, dtype=float)))


def rotate_chain_poles(poles):
    """链的频率旋转 ``z -> -jz`` 在 ``zeta = z^-2`` 上的作用：``zeta -> -zeta``，即极点取负。

    未旋转时 ``zeta`` 极点是 ``-1/rho^2``（负），旋转后变成 ``+1/rho^2``；
    换成库里的节参数写法 ``(zeta-a)/(1-a*zeta)``（分母根在 ``zeta = 1/a``）就是
    ``a = rho^2 > 0`` —— 与 ``iir_hilbert.hpp`` 的 ``APF`` 完全一致。
    """
    return -np.asarray(poles)


def analytic_from_chains(p0, p1, g0, g1, w):
    """旋转后的解析滤波器：``H = 1/2 [ A0(-zeta) + j*z^-1*A1(-zeta) ]``（即 ``H(z) = H_halfband(-j*z)``）。

    带延迟的那条链（``A1``，延迟就是它在 ``z=0`` 的极点）接**虚部**槽 —— 这是唯一的解析方向接法：
    把延迟换到实部槽会得到频率镜像（正频 0 / 负频 1），因为那等价于对系数取共轭，而系数共轭把
    ``w -> -w`` 镜像掉。也与 ``iir_hilbert.hpp`` 的 ``Tick`` 一致。
    """
    zeta = np.exp(-2j * np.asarray(w, dtype=float))
    zinv = np.exp(-1j * np.asarray(w, dtype=float))
    return g0 * zeta_chain_arg(p0, -zeta) + g1 * 1j * zinv * zeta_chain_arg(p1, -zeta)


def find_split(a: np.ndarray, b: np.ndarray, num: int = 4001):
    """穷举把 zeta 极点分给两条链的**所有**方式（含不等分），求 ``H = g0*A0 + g1*z^-1*A1`` 的最佳解。

    :return: ``(residual, poles0, poles1, g0, g1)``
    """
    zp = zeta_poles(a)
    m = len(zp)
    w = np.linspace(1e-6, np.pi - 1e-6, num)
    target = sg.freqz(b, a, worN=w)[1]
    delay = np.exp(-1j * w)
    scale = float(np.max(np.abs(target)))
    best = None
    for d0 in range(0, m + 1):
        for sub in itertools.combinations(range(m), d0):
            rest = [i for i in range(m) if i not in sub]
            a0, a1 = zeta_chain([zp[i] for i in sub], w), zeta_chain([zp[i] for i in rest], w)
            m_ = np.stack([a0, delay * a1], axis=1)
            g, *_ = np.linalg.lstsq(m_, target, rcond=None)
            res = float(np.max(np.abs(m_ @ g - target)) / scale)
            if best is None or res < best[0]:
                best = (res, np.array([zp[i] for i in sub]), np.array([zp[i] for i in rest]),
                        g[0], g[1])
    return best


def check(N: int, rs: float) -> tuple:
    """设计 + 分解 + 校验，返回一行统计。"""
    d = halfband_design(N, rs)
    b, a = d["b"], d["a"]
    res, p0, p1, g0, g1 = find_split(a, b)
    return d, res, p0, p1, g0, g1


# ------------------------------------------------------------
# 绘图
# ------------------------------------------------------------


def plot(d: dict, p0, p1, g0, g1, res: float, w, Han, tgt, guard: float, out: Path) -> None:
    """与 ``hilbert_even.py`` 相同的 2x3 版式（第 6 格留空）：

    (0,0) 半带极点按链分组  (0,1) 半带原始 |H|(dB)（叠未旋转的 AP 链复原）  (0,2) 旋转后的实轴极点
    (1,0) 解析滤波器 |H|(dB) 与目标  (1,1) 两路相位差  (1,2) 空
    """
    b, a = d["b"], d["a"]
    N, rs = len(a) - 1, d["rs"]
    fig, ax = plt.subplots(2, 3, figsize=(17.5, 9.5))

    # (0,0) 半带极点：虚轴上的 ±rho，按 ζ 链分组
    th = np.linspace(0, 2 * np.pi, 400)
    ax[0, 0].plot(np.cos(th), np.sin(th), "k:", lw=0.8)
    for name, poles, col in (("chain A0", p0, "#1f77b4"), ("chain A1", p1, "#d62728")):
        rho = np.sqrt(-1.0 / np.real(poles))          # zeta 极点 = -1/rho^2 => rho = sqrt(-1/zeta)
        pts = np.concatenate([1j * rho, -1j * rho])
        ax[0, 0].plot(pts.real, pts.imag, "x", color=col, ms=10,
                      label=f"{name} ({len(poles)} x (+/-rho))")
    ax[0, 0].plot(0, 0, "x", color="#d62728", ms=10, label="chain A1: delay pole (z=0)")
    z = np.roots(b)
    ax[0, 0].plot(z.real, z.imag, "o", mfc="none", mec="g", ms=7, label="B zeros")
    ax[0, 0].axhline(0, color="gray", lw=0.5); ax[0, 0].axvline(0, color="gray", lw=0.5)
    ax[0, 0].set(title="halfband poles: on the imaginary axis, split per chain",
                 xlabel="Re", ylabel="Im", aspect="equal")
    ax[0, 0].grid(alpha=0.3); ax[0, 0].legend(loc="upper left", fontsize=8)

    # (0,1) 半带原始 |H|(dB) + 未旋转的「AP 链 + 延迟」复原
    ww = np.linspace(1e-9, np.pi - 1e-9, 20001)
    Hhb = sg.freqz(b, a, worN=ww)[1]
    Hrec = g0 * zeta_chain(p0, ww) + g1 * np.exp(-1j * ww) * zeta_chain(p1, ww)
    ax[0, 1].plot(ww / np.pi, 20 * np.log10(np.maximum(np.abs(Hhb), 1e-30)), lw=2.5,
                  label="halfband (direct)")
    ax[0, 1].plot(ww / np.pi, 20 * np.log10(np.maximum(np.abs(Hrec), 1e-30)), "r--", lw=1.1,
                  label=f"1/2[A0(z$^{{-2}}$)+z$^{{-1}}$A1] (max|$\\Delta$|={res:.1e})")
    ax[0, 1].axhline(-rs, color="r", ls=":", lw=1, label=f"design rs = {rs:g} dB")
    ax[0, 1].axhline(-3.01, color="k", ls=":", lw=1, label="-3.01 dB (crossover)")
    ax[0, 1].axvline(0.5, color="gray", lw=0.6)
    ax[0, 1].axvline(d["Wp"], color="g", ls="--", lw=1, label=f"Wp = {d['Wp']:.6f}")
    ax[0, 1].axvline(1 - d["Wp"], color="g", ls="--", lw=1)
    ax[0, 1].set(title="halfband LP: raw |H| (dB)", xlabel="$\\omega/\\pi$", ylabel="dB",
                 ylim=(-1.35 * rs - 6, 3))
    ax[0, 1].grid(alpha=0.3, which="both"); ax[0, 1].legend(fontsize=7, loc="lower left")

    # (0,2) 旋转后的实轴极点：每个 ±rho 对给出一对 ±rho 实极点，按链着色；另标延迟单元
    for name, poles, col, y in (("chain A0'", p0, "#1f77b4", 0.0), ("chain A1'", p1, "#d62728", 0.05)):
        rho = np.sqrt(-1.0 / np.real(poles))          # 旋转后节参数 a = rho^2，z 极点 = ±rho
        pts = np.concatenate([rho, -rho])
        ax[0, 2].plot(pts, [y] * len(pts), "x", color=col, ms=10,
                      label=f"{name} ({len(poles)} x (+/-rho))")
    ax[0, 2].plot(0, 0.05, "x", color="#d62728", ms=10, label="chain A1': delay pole (z=0)")
    ax[0, 2].plot([-1.05, 1.05], [0, 0], "k:", lw=0.8)
    ax[0, 2].set(title="after rotation: real poles (+/-rho), grouped per chain",
                 xlabel="Re (z)", ylim=(-0.4, 0.4))
    ax[0, 2].grid(alpha=0.3); ax[0, 2].legend(fontsize=8)

    # (1,0) 解析滤波器 |H|(dB) 与目标（旋转后的半带）
    ax[1, 0].plot(w / np.pi, 20 * np.log10(np.maximum(np.abs(Han), 1e-30)), lw=2,
                  label="A_real + j z$^{-1}$A_imag (rotated)")
    ax[1, 0].plot(w / np.pi, 20 * np.log10(np.maximum(np.abs(tgt), 1e-30)), "k--", lw=1,
                  label="target: halfband rotated by $-\\pi/2$")
    ax[1, 0].axhline(0, color="k", ls=":", lw=1)
    ax[1, 0].set(title="analytic filter (dB): stopband = halfband stopband (= rs)",
                 xlabel="$\\omega/\\pi$", ylabel="dB", ylim=(-1.35 * rs - 6, 3))
    ax[1, 0].grid(alpha=0.3, which="both"); ax[1, 0].legend(fontsize=7, loc="lower left")

    # (1,1) 两链相位差（含结构延迟、不含权重 j）：两带内恒为 ±90°，符号随频带翻转
    zeta_ = np.exp(-2j * w)
    A0r_p, A1r_p = zeta_chain_arg(p0, -zeta_), zeta_chain_arg(p1, -zeta_)
    # 两带内平直于 -90°（正频）/ +90°（负频）。曲线在 DC 与 ±pi 处是 *连续* 穿过 ±180 分支割线的，
    # 折返成 ±180 后若直接连线会画出一条假竖直跳线 —— 所以只在这几个样本断开，不改变折返值本身。
    ratio = np.degrees(np.angle(np.exp(-1j * w) * A1r_p / A0r_p))
    ratio[np.abs(w) < 1e-12] = np.nan                    # DC
    ratio[np.abs(np.abs(w) - np.pi) < 1e-12] = np.nan    # ±pi
    ax[1, 1].plot(w / np.pi, ratio, lw=2, label="arg(z$^{-1}$A1') - arg A0'")
    ax[1, 1].axhline(-90, color="k", ls=":", label="-90 deg")
    ax[1, 1].axhline(90, color="gray", ls=":", label="+90 deg")
    ax[1, 1].set(title="branch phase difference: -90 deg (w>0) / +90 deg (w<0)",
                 xlabel="$\\omega/\\pi$", ylabel="deg", ylim=(-200, 200))
    ax[1, 1].grid(alpha=0.3); ax[1, 1].legend(fontsize=8)

    # (1,2) 留空（与 hilbert_even.py 同版式）
    ax[1, 2].axis("off")

    fig.suptitle(f"odd halfband: two zeta allpass chains + delay (N={N}, rs={rs:g} dB) "
                 f"-> rotate -> iir_hilbert form", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"图已保存：{out}")


# ------------------------------------------------------------
# main
# ------------------------------------------------------------


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="odd-order halfband = two zeta allpass chains + delay")
    p.add_argument("--N", type=int, default=15, help="半带阶数（必须为奇数）")
    p.add_argument("--rs", type=float, default=60.0, help="阻带深度 dB")
    p.add_argument("--check", action="store_true", help="扫几组 (N, rs)")
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args(argv)

    if args.check:
        print(f"{'N':>4}{'rs':>7}{'zeta阶(链0+链1)':>16}{'残差':>11}{'g0':>10}{'g1':>10}")
        for N, rs in ((9, 30), (11, 40), (15, 60), (17, 80), (21, 60), (25, 80)):
            try:
                d, res, p0, p1, g0, g1 = check(N, rs)
                print(f"{N:>4}{rs:>7}{f'{len(p0)} + {len(p1)}':>16}{res:>11.2e}"
                      f"{g0.real:>10.6f}{g1.real:>10.6f}")
            except Exception as e:
                print(f"{N:>4}{rs:>7}   失败: {type(e).__name__}: {str(e)[:45]}")
        return

    if args.N % 2 == 0:
        p.error("--N 必须为奇数（偶阶见 hilbert_even.py）")

    d, res, p0, p1, g0, g1 = check(args.N, args.rs)
    b, a = d["b"], d["a"]
    info = halfband_report(d)
    print(f"== 奇阶半带 N={args.N}, rs={args.rs:g} dB -> rp={d['rp']:.3e} dB, Wp={d['Wp']:.9f} ==")
    print(f"  A 末位系数 a[-1]      = {abs(a[-1]):.2e}（奇阶 => 0，有效阶 N-1）")
    print(f"  |H(0.5)|             = {info['h_half']:.9f}（1/sqrt(2)）")
    print(f"  功率互补误差          = {info['power_comp']:.2e}")
    print()
    print("== 两组 AP 链 + 延迟：H = g0*A0(z^-2) + g1*z^-1*A1(z^-2) ==")
    print(f"  链 0：zeta 阶 {len(p0)}，zeta 极点 = {np.round(np.real(p0), 6)}")
    print(f"  链 1：zeta 阶 {len(p1)}，zeta 极点 = {np.round(np.real(p1), 6)}")
    print(f"  权重 g0 = {g0:.9g}   g1 = {g1:.9g}（应恰为 0.5 / 0.5）")
    print(f"  z 极点计数：链0 {2*len(p0)} 个（{len(p0)} 组 ±rho）；链1 {2*len(p1)} + 1（z=0 延迟极点）"
          f" = {2*len(p1)+1} 个；合计 {2*(len(p0)+len(p1))+1} = N = {2*(len(p0)+len(p1))+1}")
    print(f"  复原残差（相对）      = {res:.2e}")
    print(f"  节点参数 a = rho^2    : 链0 {np.round(-1/np.real(p0), 6)} / 链1 {np.round(-1/np.real(p1), 6)}")
    print()

    # ---------- 旋转：zeta -> -zeta，延迟项变成 j*z^-1 ----------
    w = np.linspace(-np.pi, np.pi, 16001)
    guard = max(0.05, 0.5 - d["Wp"])
    Han = analytic_from_chains(p0, p1, g0, g1, w)
    st = band_stats(Han, w, guard)
    tgt = sg.freqz(b, a, worN=w - np.pi / 2)[1]            # 目标 = 旋转后的半带
    zeta_ = np.exp(-2j * w)
    A0r, A1r = zeta_chain_arg(p0, -zeta_), zeta_chain_arg(p1, -zeta_)
    sb = (np.abs(w) > guard * np.pi) & (np.abs(w) < (1 - guard) * np.pi)
    ratio_r = np.exp(-1j * w) * A1r / A0r                      # 两链（含结构延迟、不含 j）
    delta = np.angle(ratio_r)                                  # 统计用折返值（图里才 unwrap）
    core = (np.abs(w) > 0.2 * np.pi) & (np.abs(w) < 0.8 * np.pi)
    dev = np.max(np.abs(np.abs(delta[core & (w > 0)]) - np.pi / 2)) * 180 / np.pi
    print("== 旋转后（= iir_hilbert.hpp 的形式）：A_real(z) + j*z^-1*A_imag(z) ==")
    print(f"  旋转后节参数 a = rho^2 : 链0 {np.round(-1.0 / np.real(p0), 6)} / "
          f"链1 {np.round(-1.0 / np.real(p1), 6)}（全为正 => 与 iir_hilbert.hpp 的 APF 同号）")
    print(f"  与目标（旋转半带）最大偏差 = {np.max(np.abs(Han - tgt)):.2e}"
          f"（延迟链在虚部槽 = 解析方向）")
    print(f"  正频 |H| = [{st['pos_min']:.5f}, {st['pos_max']:.5f}]"
          f"（{20*np.log10(st['pos_min']):.4f} … {20*np.log10(st['pos_max']):.4f} dB，理想 1/0）")
    print(f"  负频最大 |H| = {st['neg_max']:.3e} = {20*np.log10(st['neg_max']):.2f} dB"
          f" @ w = {st['neg_argmax_w']/np.pi:+.3f}pi（= 设计 rs={args.rs:g} dB）")
    print(f"  两链相位差恒为 ±90 度：正频核心段偏离 {dev:.3f} 度")
    print()

    out = args.out or SCRIPT_DIR / "output" / f"hilbert_odd_N{args.N}_rs{args.rs:g}.png"
    plot(d, p0, p1, g0, g1, res, w, Han, tgt, guard, out)


if __name__ == "__main__":
    main()
