"""「消除负频率能不能解决互调」——把两条命题变成可验证的东西。

命题 A：线性插值上采样后，信号里既有正频率镜像、也有**负频率**镜像。
命题 B：在上采样**之前**把负频率去掉（即把低采样率信号变成解析信号，只留正频率），
        线性插值**仍会造出负频率**。

验证方式（fs = 1000 Hz、纯音 300 Hz）：
  ① 实信号线性插值上采样 → 双边谱（线在 ±300、±700、±1300 …，对称）；
  ② 解析信号（低采样率下只有 +300）线性插值上采样 → 双边谱：线在 300+1000m，
     其中 m = −1 落在 **−700 Hz（真正的负频率）**，且 ±300 的对称性没了；
  ③ 关键一条数值事实：**Ref(解析上采样) 与实信号上采样逐点相同** —— 对实数非线性
     来说"去掉负频率"什么也没改变（负频率本来就是正频率的共轭镜像，不含独立信息）。

结论：要消的不是"负频率"，而是**镜像（m ≠ 0 的整套复制品）**；那要靠带限重建
（sinc 类插值器）而不是解析化。另外基波自己的镜像 −f0 必须保留 —— 偶数次谐波就是
f0 − (−f0) = 2f0 这样来的。

用法: python qwqdsp/labs/adaa_iir/negative_freq_question.py
输出: qwqdsp/labs/adaa_iir/output/negative_freq_question.png
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
RATIO = 64                      # 图上看镜像够了，不必 8192
NREPS = 8                       # 取整数个「最小周期」避免泄漏：fs/f0 = 10/3 → 10 个样本 = 3 个周期


def main() -> int:
    n_low = 10 * NREPS          # 80 个低采样率样本 = 24 个周期
    x = np.sin(2.0 * np.pi * F0 * np.arange(n_low) / FS_LOW)

    # 解析信号（低采样率下：只保留正频率）
    x_an = ss.hilbert(x)

    def interp(sig):
        s = np.concatenate([sig, sig[:1]])
        idx = np.arange(len(sig) * RATIO) / RATIO
        return np.interp(idx, np.arange(len(sig) + 1), s)

    x_up = interp(x)
    x_an_up = interp(x_an.real) + 1j * interp(x_an.imag)

    print(f"fs_low = {FS_LOW:.0f} Hz，f0 = {F0:.0f} Hz，上采样 ×{RATIO}（{n_low} 个低采样率样本）")
    print("\n① 实信号上采样后的双边线谱（前 6 条）:")
    f = np.fft.fftfreq(len(x_up), 1.0 / (FS_LOW * RATIO))
    sp = np.abs(np.fft.fft(x_up))
    ok = sp > sp.max() * 1e-3
    for ff in sorted(f[ok], key=abs)[:6]:
        print(f"   {ff:+9.1f} Hz   {20*np.log10(sp[np.argmin(abs(f-ff))]/sp.max()):6.1f} dB")

    print("\n② 解析信号上采样后的双边线谱（前 6 条）—— 注意负频率仍存在:")
    sp2 = np.abs(np.fft.fft(x_an_up))
    ok2 = sp2 > sp2.max() * 1e-3
    for ff in sorted(f[ok2], key=abs)[:6]:
        print(f"   {ff:+9.1f} Hz   {20*np.log10(sp2[np.argmin(abs(f-ff))]/sp2.max()):6.1f} dB")

    same = np.max(np.abs(x_an_up.real - x_up))
    print(f"\n③ max|Ref(解析上采样) − 实信号上采样| = {same:.3e}"
          f"   → {'逐点相同，实数非线性看到的信号没变' if same < 1e-12 else '不一致（预期应相同）'}")

    # ---- 图：双边谱 ----
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 6.4), sharex=True)
    for ax, (sig, title) in zip(axes, [
        (x_up, f"① 实信号：线性插值 ×{RATIO}（线在 ±300、±700、±1300 …，对称）"),
        (x_an_up, f"② 解析信号（低采样率下只剩 +300）：上采样后仍有负频率（−700、−1700 …）"),
    ]):
        sp_ = np.abs(np.fft.fft(sig)) * 2.0 / len(sig)
        db = 20 * np.log10(np.maximum(sp_, 1e-6) / sp_.max())
        loud = db > -60
        ax.vlines(f[loud], -120, db[loud], color="tab:blue", lw=1.6)
        ax.set_ylim(-120, 5)
        ax.set_xlim(-2500, 2500)
        ax.set_xlabel("频率 [Hz]（双边，线性轴）")
        ax.set_ylabel("幅度 [dB]")
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
        ax.axhline(-60, color="0.7", ls=":", lw=1)
    axes[0].text(-2400, -20, "基波 −f0\n（偶数次谐波靠它）", fontsize=8, color="tab:red")
    axes[0].annotate("", xy=(-300, -5), xytext=(-1500, -18),
                     arrowprops=dict(arrowstyle="->", color="tab:red", lw=1))
    axes[1].text(-2400, -20, "−700 Hz：解析化之后\n新生的负频率镜像", fontsize=8, color="tab:red")
    axes[1].annotate("", xy=(-700, -5), xytext=(-1500, -18),
                     arrowprops=dict(arrowstyle="->", color="tab:red", lw=1))
    fig.suptitle("「消除负频率」为什么不是解药：实数信号的负频率是共轭镜像（去不掉也不该去），"
                 "而解析化之后线性插值照样造出负频率镜像", fontsize=10.5)
    fig.tight_layout()
    out = Path(__file__).parent / "output" / "negative_freq_question.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"\n已写出 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
