"""Vicanek `match_biquad.hpp` 与 Ivantsov 双点解的频率映射对照。

Vicanek 的极点是脉冲响应不变（角度 = wc），零点由幅度平方在 2~3 个频点上的
线性方程解出，所以**没有全局映射**：每种类型的等价映射曲线都不同，而且
低通 / 高通只匹配 DC 与 wc（Nyquist 不管），映射右端停在 0.84×Nyquist 左右。
Ivantsov 的双点解则是一条与类型无关的曲线，右端精确落在 Nyquist。

映射定义与 `../ivantsovy/exp_mapping.py` 一致：给定设计，找模拟频率 y 使 ``|Ha(y)| = |Hd(x)|``。

用法: python qwqdsp/labs/match_filter/vicanek/exp_vicanek.py
产物: output/vicanek_mapping.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "ivantsovy"))

import ivantsovy as I  # noqa: E402
import two_point as TP  # noqa: E402
import vicanek as V  # noqa: E402

FS = 2.0
OUT = Path(__file__).parent / "output"
ZETA = 1.0 / np.sqrt(2.0)
Q = 1.0 / (2.0 * ZETA)
B2 = 2.0 * ZETA ** 2 - 1.0
BOOST = 12.0

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False


def inv_lp(m):
    """二阶低通原型幅度反解（``|Ha|^2 = 1/(r^4+2B r^2+1)``）。"""
    v = 1.0 / np.clip(m, 1e-14, 1.0 - 1e-15) ** 2 - 1.0
    return np.sqrt(np.maximum(-B2 + np.sqrt(B2 * B2 + v), 0.0))


def inv_hp(m):
    m2 = np.clip(m, 1e-14, 1.0 - 1e-12) ** 2
    a = 1.0 / m2 - 1.0
    return np.sqrt(np.maximum((-B2 + np.sqrt(B2 * B2 + a)) / a, 0.0))


def mapping_vicanek_lp(wc: float, n: int = 200001):
    w = np.geomspace(1e-7, np.pi, n)
    b, a = V.lowpass(wc, Q)
    return w / np.pi, wc * inv_lp(V.magnitude(b, a, w)) / np.pi


def mapping_vicanek_hp(wc: float, n: int = 200001):
    w = np.geomspace(1e-7, np.pi, n)
    b, a = V.highpass(wc, Q)
    return w / np.pi, wc * inv_hp(V.magnitude(b, a, w)) / np.pi


def mapping_vicanek_shelf(wc: float, n: int = 200001):
    """shelf 的模拟幅度单调增（Q = 0.707），直接对幅值反查频率。"""
    w = np.geomspace(1e-7, np.pi, n)
    b, a = V.highshelf(wc, Q, BOOST)
    ha = np.sqrt(V.shelf_target(BOOST, wc, Q, "high")(w))
    hd = V.magnitude(b, a, w)
    # 只在该档数字幅度落在模拟原型值域内时才有解；越界处留空
    inside = (hd >= ha[0]) & (hd <= ha[-1] * 1.000001)
    y = np.full_like(w, np.nan)
    y[inside] = np.interp(hd[inside], ha, w)
    return w / np.pi, y / np.pi


def mapping_ivantsov(fc: float, n: int = 200001):
    """Ivantsov 双点解：与类型无关，取二阶低通即可。"""
    w = np.geomspace(1e-7, np.pi, n)
    x = w / np.pi
    hd = I.digital_magnitude(TP.lowpass2(fc, FS, ZETA), x, FS)
    return x, fc * inv_lp(hd)


def main() -> int:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.9))
    for ax, fcn in zip(axes, (0.2, 0.9)):
        wc = np.pi * fcn
        ax.plot([0.0, 1.0], [0.0, 1.0], color="0.4", lw=1.0, ls="--", label="恒等映射")
        for name, fn, color in (("Vicanek 低通", mapping_vicanek_lp, "C3"),
                                ("Vicanek 高通", mapping_vicanek_hp, "C1"),
                                ("Vicanek shelf +12 dB", mapping_vicanek_shelf, "C4")):
            x, y = fn(wc)
            ax.plot(x, y, color=color, lw=1.7, label=name)
        x, y = mapping_ivantsov(fcn)
        ax.plot(x, y, color="C0", lw=2.2, label="Ivantsov 双点解（与类型无关）")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel("数字频率 / (fs/2)")
        ax.set_ylabel("等效模拟频率 / (fs/2)")
        ax.set_title(f"等价映射（fc = {fcn:g}·fs/2，Q = 1/√2）")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    path = OUT / "vicanek_mapping.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"图：{path}")

    print("映射 y/(fs/2)（x = 0.2 / 0.5 / 0.8 / 0.95 / 1.0）")
    for fcn in (0.2, 0.9):
        wc = np.pi * fcn
        for name, fn in (("Vicanek-LP", mapping_vicanek_lp),
                         ("Vicanek-HP", mapping_vicanek_hp),
                         ("Vicanek-shelf", mapping_vicanek_shelf),
                         ("Ivantsov", mapping_ivantsov)):
            x, y = fn(fcn)
            vals = [np.interp(t, x, y) for t in (0.2, 0.5, 0.8, 0.95, 0.99999)]
            print(f"  fc={fcn:4.1f} {name:<14}" +
                  " ".join(f"{v:7.4f}" if np.isfinite(v) else "      —" for v in vals))

    print("\n全带 max|ΔdB|（fs=48k，Q=1/√2；shelf 用同一模拟原型）")
    print(f"{'fc':>8}{'V-LP':>8}{'I-LP':>8}{'V-HP':>8}{'I-HP':>8}{'V-shelf':>9}{'I-shelf':>9}")
    for fc in (20.0, 117.6, 691.4, 4065.0, 23900.0):
        fs = 48000.0
        wc = 2.0 * np.pi * fc / fs
        w = np.linspace(1e-9, np.pi, 40001)
        fhz = w / np.pi * fs / 2.0
        r = fhz / fc
        ha_lp = 1.0 / np.sqrt(r ** 4 + 2.0 * r ** 2 * B2 + 1.0)
        ha_hp = r ** 2 * ha_lp
        cells = []
        for fn, ha in ((V.lowpass, ha_lp), (V.highpass, ha_hp)):
            b, a = fn(wc, Q)
            dev = 20 * np.log10(np.maximum(V.magnitude(b, a, w), 1e-14) / ha)
            cells.append(float(np.max(np.abs(dev[ha >= 1e-6]))))
        for build, ha in ((TP.lowpass2, ha_lp), (TP.highpass2, ha_hp)):
            hd = I.digital_magnitude(build(fc, fs, ZETA), fhz, fs)
            dev = 20 * np.log10(np.maximum(hd, 1e-300) / ha)
            cells.append(float(np.max(np.abs(dev[ha >= 1e-6]))))

        ha_sh = np.sqrt(V.shelf_target(BOOST, wc, Q, "high")(w))
        b, a = V.highshelf(wc, Q, BOOST)
        dev = 20 * np.log10(np.maximum(V.magnitude(b, a, w), 1e-14) / ha_sh)
        cells.append(float(np.max(np.abs(dev))))
        zeros, poles = V.shelf_roots(BOOST, Q)
        omega, sigma = TP.omega_sigma(fc, fs)
        design = I.design(zeros, poles, fs, fc, sigma, omega=omega, gain=1.0)
        hd = I.digital_magnitude(design, fhz, fs)
        dev = 20 * np.log10(np.maximum(hd, 1e-300) / ha_sh)
        cells.append(float(np.max(np.abs(dev))))
        print(f"{fc:8.1f}" + "".join(f"{v:8.3f}" for v in cells[:4]) +
              "".join(f"{v:9.3f}" for v in cells[4:]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
