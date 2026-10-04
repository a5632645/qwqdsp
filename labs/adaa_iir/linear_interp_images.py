"""线性插值上采样的镜像频率：为什么输出端 AA 滤波器管不到它们。

三步构造（fs_low 取超低值，方便把镜像看得清清楚楚）：
  1. 低采样率纯音 x[n] = sin(2π f0 n / fs_low)   —— 频谱用**线性频率轴**；
  2. 线性插值上采样 8192 倍 → x_up             —— 频谱用**对数频率轴**；
  3. 对 x_up 施加失真（硬削波）→ y_up           —— 频谱用**对数频率轴**。

信号都是严格周期的，所以用「一个周期长度的矩形窗 FFT」得到的是**无泄漏的线谱**（不需要窗函数）。
第 3 张图额外叠上：输出端 AA 滤波器的幅频 |H|（虚线）与「过完 AA 滤波器之后」的谱 ——
可以看到 AA 把高频镜像压掉，但**落在输出 Nyquist 以内的互调产物一动不动**，这就是它修不好的原因。

用法: python qwqdsp/labs/adaa_iir/linear_interp_images.py
输出: qwqdsp/labs/adaa_iir/output/linear_interp_images.png
"""
from __future__ import annotations

import sys
from math import gcd
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np


FS_LOW = 1000.0        # 超低采样率 [Hz]
F0 = 300.0             # 纯音频率 [Hz]（fs/F0 = 10/3：镜像落在 100 Hz 细格点，但不在谐波 300k 上——这是能看到带内互调的关键）
RATIO = 8192           # 线性插值上采样倍数
DRIVE = 2.0            # 失真前的增益

FS_UP = FS_LOW * RATIO
PERIOD = FS_LOW / gcd(int(FS_LOW), int(F0))     # 低采样率下的周期（样本数）
FLOOR_DB = -120.0


def tone_one_period():
    """低采样率下一个周期的纯音。"""
    n = int(PERIOD)
    return np.sin(2.0 * np.pi * F0 * np.arange(n) / FS_LOW)


def linear_upsample(x, ratio):
    """线性插值上采样（含周期性收尾：最后一段用 x[0] 闭合）。"""
    xs = np.concatenate([x, x[:1]])
    idx = np.arange(len(x) * ratio) / ratio
    return np.interp(idx, np.arange(len(x) + 1), xs)


def line_spectrum(sig, fs):
    """一个周期长度矩形窗 FFT → 无泄漏线谱（幅度，单边归一化到正弦峰值）。"""
    n = len(sig)
    sp = np.abs(np.fft.rfft(sig)) * 2.0 / n
    return np.fft.rfftfreq(n, 1.0 / fs), sp


def to_db(sp, ref):
    return 20.0 * np.log10(np.maximum(sp, 1e-30) / ref)


def bin_columns(freqs, db, f_lo, f_hi, num_cols):
    """按对数频率分列取最大值（线谱很密时画成包络）。"""
    edges = np.logspace(np.log10(f_lo), np.log10(f_hi), num_cols + 1)
    out_x = np.sqrt(edges[:-1] * edges[1:])
    out_y = np.full(num_cols, FLOOR_DB)
    idx = np.searchsorted(edges, freqs, side="right") - 1
    for k in range(num_cols):
        m = idx == k
        if m.any():
            out_y[k] = max(out_y[k], db[m].max())
    return out_x, out_y


def aa_filter_response(freqs_hz):
    """输出端 AA 滤波器的幅频响应（10 阶椭圆 rs=80dB、阻带边沿 0.49·fs_low）。

    这里用 scipy 的标准椭圆原型（偶数阶修正只改分子/阻带边沿，对本图结论无影响），
    自行把「阻带边沿」定标到 0.49·fs_low，与 demo 里那一档同规格。
    """
    import scipy.signal as ss
    z, p, k = ss.ellip(10, 0.1, 80.0, 1.0, btype="low", analog=True, output="zpk")
    prot = ss.freqs_zpk(z, p, k, np.logspace(-2, 1.2, 200000))
    w_stop = prot[0][np.argmax(20 * np.log10(np.abs(prot[1])) <= -80.0)]
    scale = 2.0 * np.pi * 0.49 / w_stop
    _, h = ss.freqs_zpk(z, p, k, 2.0 * np.pi * np.asarray(freqs_hz) / FS_LOW / scale)
    h0 = ss.freqs_zpk(z, p, k, np.array([1e-9]))[1][0]
    return np.abs(h) / abs(h0)


def main() -> int:
    x_low = tone_one_period()
    x_up = linear_upsample(x_low, RATIO)
    y_up = np.clip(DRIVE * x_up, -1.0, 1.0)

    f_low, sp_low = line_spectrum(x_low, FS_LOW)
    f_up, sp_up = line_spectrum(x_up, FS_UP)
    f_dis, sp_dis = line_spectrum(y_up, FS_UP)

    ref_low = sp_low.max()
    ref_dis = sp_dis.max()
    db_up = to_db(sp_up, sp_up.max())
    db_dis = to_db(sp_dis, ref_dis)

    # AA 滤波器（按 fs_low 归一化）作用在失真谱上（频域相乘 = 连续时间滤波的等效）
    h_aa = aa_filter_response(f_up)
    sp_filtered = sp_dis * h_aa
    db_filt = to_db(sp_filtered, ref_dis)

    print(f"低采样率 {FS_LOW:.0f} Hz，纯音 {F0:.0f} Hz，线性插值上采样 ×{RATIO} → {FS_UP/1e6:.3f} MHz")
    print(f"周期 {PERIOD:.0f} 个低采样率样本（= {len(x_up)} 个上采样样本，线谱间距 {F0:.0f} Hz）\n")

    print("线性插值镜像（相对基音 dB）—— 实测 vs 理论 sinc² 包络:")
    print(f"{'镜像频率[Hz]':>14}{'实测[dB]':>10}{'sinc² 包络[dB]':>16}")
    for m in (1, 2, 3, 4):
        for f_img, tag in ((m * FS_LOW - F0, f"{m}·fs−f0"), (m * FS_LOW + F0, f"{m}·fs+f0")):
            i = int(np.argmin(np.abs(f_up - f_img)))
            # 线性插值核 = 宽度 2 个低采样率样本的三角 → 幅度响应 sinc²(f/fs_low)
            env = 20 * np.log10(abs(np.sinc(f_img / FS_LOW)) ** 2)
            print(f"{tag:>10} {f_img:>8.0f}{db_up[i]:>10.1f}{env:>16.1f}")

    # 输出 Nyquist（折叠边界）以内的非谐波产物
    inband = f_up < FS_LOW / 2
    harm = np.zeros_like(f_up, dtype=bool)
    for k in range(1, int(FS_LOW / 2 / F0) + 1):
        harm |= np.abs(f_up - k * F0) < F0 * 0.4
    prod = inband & ~harm
    print(f"\n输出 Nyquist（{FS_LOW/2:.0f} Hz）以内的最大非谐波产物：")
    print(f"  过 AA 滤波器之前 {db_dis[prod].max():7.1f} dB   之后 {db_filt[prod].max():7.1f} dB"
          f"   ← 滤波器压不动它们")

    # ---- 画图 ----
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 12))

    ax = axes[0]
    db_low = to_db(sp_low, ref_low)
    loud = db_low > FLOOR_DB + 1
    ax.vlines(f_low[loud], FLOOR_DB, db_low[loud], color="tab:blue", lw=1.6)
    ax.set_xlim(0, FS_LOW / 2)
    ax.set_ylim(FLOOR_DB, 5)
    ax.set_xlabel("频率 [Hz]（线性轴）")
    ax.set_ylabel("幅度 [dB]")
    ax.set_title(f"① 低采样率纯音：{F0:.0f} Hz @ {FS_LOW:.0f} Hz（fs/2 = {FS_LOW/2:.0f} Hz）")
    ax.grid(alpha=0.3)

    for ax, (f, db, title, extra) in zip(axes[1:], [
        (f_up, db_up, f"② 线性插值上采样 ×{RATIO}（{FS_UP/1e6:.3f} MHz 采样）", None),
        (f_dis, db_dis, f"③ 上采样信号过失真（硬削波 ×{DRIVE:.0f}）", (h_aa, db_filt)),
    ]):
        fx, fy = bin_columns(f, db, F0 / 4, FS_UP / 2, 2000)
        ax.semilogx(fx, fy, color="tab:blue", lw=0.9, label="频谱")
        env_f = np.logspace(np.log10(F0 / 4), np.log10(FS_UP / 2), 2000)
        env_db = 20 * np.log10(np.maximum(np.abs(np.sinc(env_f / FS_LOW)) ** 2, 1e-12))
        ax.semilogx(env_f, np.maximum(env_db, FLOOR_DB), color="0.35", ls=":", lw=1.0,
                    label="插值核包络 sinc²(f/fs_low)")
        if extra is not None:
            h, dbf = extra
            fxm, fym = bin_columns(f, dbf, F0 / 4, FS_UP / 2, 2000)
            ax.semilogx(fxm, fym, color="tab:orange", lw=0.9, label="过 AA 滤波器之后")
            fh, hy = bin_columns(f, to_db(np.maximum(h, 1e-12), 1.0), F0 / 4, FS_UP / 2, 2000)
            ax.semilogx(fh, hy, color="0.55", ls="--", lw=1.0, label="AA 滤波器 |H|")
            ax.legend(fontsize=8, loc="lower left")
        ax.axvline(FS_LOW / 2, color="tab:red", ls=":", lw=1.4)
        ax.axvspan(F0 / 4, FS_LOW / 2, color="tab:red", alpha=0.06)
        ax.text(FS_LOW / 2 * 1.05, 0.0, "输出 Nyquist\n(折叠边界)", color="tab:red", fontsize=8, va="top")
        ax.set_xlim(F0 / 4, FS_UP / 2)
        ax.set_ylim(FLOOR_DB, 5)
        ax.set_xlabel("频率 [Hz]（对数轴）")
        ax.set_ylabel("幅度 [dB]")
        ax.set_title(title)
        ax.grid(True, which="both", alpha=0.25)

    fig.suptitle("线性插值的镜像频率，以及为什么输出端 AA 滤波器压掉它们也没用"
                 "（落在 Nyquist 以内的互调产物照样留下）", fontsize=11)
    fig.tight_layout()
    out = Path(__file__).parent / "output" / "linear_interp_images.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"\n已写出 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
