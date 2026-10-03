"""直接项（d ≠ 0 的原型）在三种比率下的行为。

复刻实现取输出的方式（`ResampleIIR::Process` 的相位累积），对比两种**模型**：

  dc_model : HP 原文的模型「输入 = 狄拉克串」—— 直接项是 d × 冲激串，只在 t_n 恰好落在
             输入栅格上时非零（值 d·x[k]），其余时刻为 0。这是 algorithm.md §9.4 的规则。
  dc_ideal : 另一种解释「输入 = 带限连续信号」—— 直接项是 d·x_buf(t_n)，处处非零且光滑。
             x_buf 用纯音的正弦（f < f_in/2 时采样正弦的带限插值就是它本身）。

两者只在栅格点之间有差别（基带完全一致，差别全在镜像上，见 algorithm.md §9.2）。
按 HP 模型实现重采样器时用 dc_model（§9.4）；本图用来量化两种解释差多少。


用法: python qwqdsp/labs/holters_parker/plot_direct_term.py
输出: qwqdsp/labs/holters_parker/output/direct_term_*.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = 1e-2          # 直接项的大小（= 偶数阶椭圆 rs=40 dB 的阻带波纹峰值）
F_TONE = 10000.0
N_IN = 4096
OUT_DIR = Path(__file__).parent / "output"


def run(f_in, f_out, n_in=N_IN, f=F_TONE, d=D):
    """按实现的方式累积相位，返回 (t_n, 命中掩码, dc_model, dc_ideal)。"""
    inc = f_in / f_out
    x = np.sin(2.0 * np.pi * f * np.arange(n_in) / f_in)

    times, rpos, phase = [], 0, 0.0
    while rpos < n_in - 1:
        times.append(rpos + phase)          # 输出在取样本之后才推进相位
        phase += inc
        new_rpos = min(rpos + int(np.floor(phase)), n_in - 1)
        phase -= np.floor(phase)
        rpos = new_rpos
    t = np.array(times)

    k = np.floor(t)
    hit = (t == k)
    dc_model = np.where(hit, d * x[np.minimum(k.astype(int), n_in - 1)], 0.0)
    dc_ideal = d * np.sin(2.0 * np.pi * f * t / f_in)
    return t, hit, dc_model, dc_ideal


def spectrum(sig, fs, nfft=None):
    n = len(sig)
    w = np.hanning(n)
    sp = np.abs(np.fft.rfft(sig * w, nfft)) * 2.0 / w.sum()
    return np.fft.rfftfreq(n, 1.0 / fs), sp


def figure(f_in, f_out):
    t, hit, dc_model, dc_ideal = run(f_in, f_out)
    err = np.abs(dc_model - dc_ideal).max()

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.2))

    # (a) 时间域放大
    nz = 48
    ax[0].plot(t[:nz], dc_ideal[:nz], "-", lw=1.0, color="tab:orange", label="带限模型  d·x_buf(t_n)")
    ml, sl, bl = ax[0].stem(t[:nz], dc_model[:nz], linefmt="C0-", markerfmt="C0o",
                            basefmt=" ", label="HP 模型（冲激串取样本）")
    plt.setp(ml, lw=1.2, ms=3.5)
    ax[0].plot(t[:nz][hit[:nz]], dc_model[:nz][hit[:nz]], "r.", ms=7, zorder=5)
    ax[0].axhline(0, color="0.7", lw=0.6)
    ax[0].set_xlabel("output time  (input samples)")
    ax[0].set_ylabel("direct-term contribution")
    ax[0].set_title(f"(a) time domain, first {nz} outputs")
    ax[0].legend(loc="lower right", fontsize=8)

    # (b) 频谱
    for sig, lab, c in ((dc_ideal, "ideal", "tab:orange"), (dc_model, "model", "C0")):
        fr, sp = spectrum(sig, f_out)
        ax[1].plot(fr / 1e3, 20 * np.log10(np.maximum(sp, 1e-12)), lw=1.2, color=c, label=lab)
    ax[1].set_xlim(0, f_out / 2 / 1e3)
    ax[1].set_ylim(-120, 0)
    ax[1].set_xlabel("frequency  (kHz)")
    ax[1].set_ylabel("dB (rel. full scale 1.0)")
    ax[1].set_title("(b) spectrum of the direct term")
    ax[1].legend(loc="upper right", fontsize=8)
    ax[1].grid(alpha=0.3)
    if f_out > f_in:                     # 插值：镜像落在 f_in − f
        ax[1].axvline((f_in - F_TONE) / 1e3, color="r", ls=":", lw=1.0)
        ax[1].text((f_in - F_TONE) / 1e3, -8, f" image\n {f_in/1e3:.0f}k−{F_TONE/1e3:.0f}k",
                   color="r", fontsize=8, va="top")

    # (c) 小数偏移（命中结构）
    m = min(300, len(t))
    ax[2].plot(t[:m] - np.floor(t[:m]), lw=0.9, color="C2")
    ax[2].plot(t[:m][hit[:m]] - np.floor(t[:m][hit[:m]]), "r.", ms=3)
    ax[2].set_ylim(-0.05, 1.05)
    ax[2].set_xlabel("output index n")
    ax[2].set_ylabel("frac(t_n) = d_n")
    ax[2].set_title("(c) fractional offset  (red = hits, frac == 0)")
    ax[2].grid(alpha=0.3)

    fig.suptitle(f"{f_in/1e3:.1f}k → {f_out/1e3:.1f}k   (inc = {f_in/f_out:.6f})   "
                 f"hit rate = {100*hit.mean():.3f}%   max|model − ideal| = {err:.3e}",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    OUT_DIR.mkdir(exist_ok=True)
    out = OUT_DIR / f"direct_term_{int(f_in/1000)}k_to_{int(f_out/1000)}k.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"{out}   hits={hit.sum()}/{len(hit)} ({100*hit.mean():.3f}%)  max|model-ideal|={err:.3e}")


if __name__ == "__main__":
    for pair in ((96000, 48000), (48000, 96000), (48000, 44100)):
        figure(*pair)
