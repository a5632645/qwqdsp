"""复数整形器（Vicanek, *Complex Waveshapers*, 2025）能否清掉输出奈奎斯特带内的失真频率？

按用户给的链路：低采样率纯音 → 线性插值上采样 8192× → 转解析信号 → **复数**失真 → 取实部，
然后看输出 Nyquist 带内的「非谐波」成分。三条对照：

  (i)  实链路：实信号 + 实整形器 f(x)=arctan(x)                 —— 有镜像，差频产物落在带内
  (ii) 复链路：解析(线性插值上采样) + H(z)=arctan(z)            —— 解析化去掉了负频率镜像
  (iii)复链路：**理想正交输入** r·e^{iωt}（无镜像）+ 同一 H(z) —— 参照，应完全干净

文章里 H 取的解析函数，其 Taylor 系数就是谐波幅度：F = Re(H)/r、G = Im(H)/r 构成 Hilbert 对
（式 (2)(3)(4)）。arctan 的系数是 a_n = r^{n−1}/n（只奇次，式 (32)），所以 (iii) 的谱线电平
可以直接按公式核对。

用法: python qwqdsp/labs/adaa_iir/complex_waveshaper_check.py
输出: qwqdsp/labs/adaa_iir/output/complex_waveshaper_check.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as ss

FS_LOW = 1000.0
F0 = 300.0
RATIO = 8192
R = 0.9                      # 整形器的输入幅度参数（|r| < 1 才在 arctan 的收敛域内）
NREP = 8                     # 10 个低采样率样本 = 3 个周期 → 80 个样本 = 24 个周期


def analytic_via_hilbert(x):
    """FFT 法解析信号：只保留正频率（去掉共轭镜像）。"""
    return ss.hilbert(x)


def linear_upsample(x, ratio):
    xs = np.concatenate([x, x[:1]])
    idx = np.arange(len(x) * ratio) / ratio
    return np.interp(idx, np.arange(len(x) + 1), xs)


def inband_spurs(y, fs, f0, out_nyquist, n_show=4):
    """返回 (频率, dB) —— **输出** Nyquist 带内最强的非谐波谱线。"""
    n = len(y)
    sp = np.abs(np.fft.rfft(y)) * 2.0 / n
    fr = np.fft.rfftfreq(n, 1.0 / fs)
    db = 20.0 * np.log10(np.maximum(sp, 1e-30) / sp.max())
    inb = fr < out_nyquist
    harm = np.zeros_like(fr, dtype=bool)
    for k in range(1, int(out_nyquist / f0) + 1):
        harm |= np.abs(fr - k * f0) < f0 * 0.25
    cand = inb & ~harm & (db > -120)
    idx = np.argsort(-np.where(cand, db, -999))[:n_show]
    return [(float(fr[i]), float(db[i])) for i in idx if cand[i]], (fr, db, harm)


def main() -> int:
    n_low = 10 * NREP
    x_low = np.sin(2.0 * np.pi * F0 * np.arange(n_low) / FS_LOW)
    fs_up = FS_LOW * RATIO

    x_up = linear_upsample(x_low, RATIO)
    z_an_up = analytic_via_hilbert(R * x_up)          # (ii) 解析化（去负频率）
    t_up = np.arange(len(x_up)) / fs_up
    z_ideal = R * np.exp(2j * np.pi * F0 * t_up)      # (iii) 理想正交输入（无镜像）

    # 实整形器 / 复整形器：同一个 arctan 形状
    y_real = np.arctan(R * x_up) / R                                   # (i)
    y_cplx_an = np.arctan(z_an_up).real / R                            # (ii)
    y_cplx_id = np.arctan(z_ideal).real / R                            # (iii)

    print(f"fs_low = {FS_LOW:.0f} Hz，f0 = {F0:.0f} Hz，线性插值上采样 ×{RATIO}，r = {R}")
    print(f"H(z) = arctan(z)（解析），理论谐波幅度 a_n = r^(n-1)/n（只奇次）\n")
    print(f"{'方案':<44}{'带内最强非谐波产物':>22}")
    for tag, y in (("(i)  实链路：实信号 + arctan(x)", y_real),
                   ("(ii) 复链路：解析(线性插值) + arctan(z)", y_cplx_an),
                   ("(iii)复链路：理想正交输入 + arctan(z)", y_cplx_id)):
        spurs, _ = inband_spurs(y, fs_up, F0, FS_LOW / 2)
        s = "，".join(f"{f:.0f} Hz {v:.1f} dB" for f, v in spurs) if spurs else "（无）"
        print(f"{tag:<44}{s:>22}")

    print("\n(iii) 实测谐波电平 vs 理论 a_n = r^(n-1)/n（只奇次）:")
    n_fft = len(y_cplx_id)
    sp_id = np.abs(np.fft.rfft(y_cplx_id)) * 2.0 / n_fft
    fr_id = np.fft.rfftfreq(n_fft, 1.0 / fs_up)
    for k in (1, 3, 5, 7):
        i = int(np.argmin(np.abs(fr_id - k * F0)))
        th = 20 * np.log10(R ** (k - 1) / k)
        print(f"   {k}·f0 = {k*F0:6.0f} Hz   实测 {20*np.log10(sp_id[i]/sp_id[int(np.argmin(np.abs(fr_id-F0)))]):7.1f} dB"
              f"   理论 {th:7.1f} dB")

    # ---- 图 ----
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 9.6), sharex=True)
    for ax, (y, title) in zip(axes, [
        (y_real, "(i)  实链路：实信号 + 实 arctan —— 镜像经差频落回带内（红线以左）"),
        (y_cplx_an, "(ii) 复链路：解析化后的上采样信号 + 复 arctan —— 只剩求和项，带内干净"),
        (y_cplx_id, "(iii) 复链路：理想正交输入（无镜像）—— 参照"),
    ]):
        _, (fr, db, _) = inband_spurs(y, fs_up, F0, FS_LOW / 2)
        loud = db > -100
        ax.vlines(fr[loud], -120, db[loud], color="tab:blue", lw=0.8)
        ax.axvline(FS_LOW / 2, color="tab:red", ls=":", lw=1.4)
        ax.axvspan(0, FS_LOW / 2, color="tab:red", alpha=0.06)
        ax.set_xlim(0, 6 * FS_LOW)
        ax.set_ylim(-120, 5)
        ax.set_ylabel("幅度 [dB]")
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("频率 [Hz]（0..6·fs_low）")
    axes[0].text(FS_LOW / 2 + 30, -10, "输出 Nyquist", color="tab:red", fontsize=8)
    fig.suptitle(f"复数整形器 vs 实链路：fs_low={FS_LOW:.0f} Hz，f0={F0:.0f} Hz，"
                 f"H(z)=arctan(z)（Vicanek 2025）", fontsize=10.5)
    fig.tight_layout()
    out = Path(__file__).parent / "output" / "complex_waveshaper_check.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"\n已写出 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
