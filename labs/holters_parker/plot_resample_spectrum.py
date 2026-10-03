"""重采样的输入 / 输出频谱（并显示 d ≠ 0 的原型会把什么带进输出）。

原型用仓库真实的 BestCoeffs（25 阶椭圆，`resample_coeffs.h`，脚本直接解析该头文件）。
输出用 §algorithm.md 式 (6)(7)(8) 的模型计算（已与 C++ `ResampleIIR` 逐位对齐到 1e-16）；
`d ≠ 0` 的变体按"直接项只在 t_n 落在输入栅格时贡献 d·x[k]"实现。

用法: python qwqdsp/labs/holters_parker/plot_resample_spectrum.py
输出: qwqdsp/labs/holters_parker/output/spectrum_*.png
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
COEFF_HDR = HERE.parents[1] / "include" / "qwqdsp" / "fx" / "resample_coeffs.h"
OUT_DIR = HERE / "output"

N_IN = 8192
SKIP = 512          # 跳过启动瞬态
D = 1e-2            # 直接项（偶数阶椭圆 rs=40 dB 的阻带波纹峰值）


# ------------------------------------------------------------
# 解析 resample_coeffs.h 里的 BestCoeffs
# ------------------------------------------------------------
def load_best():
    txt = COEFF_HDR.read_text(encoding="utf-8")
    body = re.search(r"struct BestCoeffs \{(.*?)\n\};", txt, re.S).group(1)

    def cf(key):
        blk = re.search(rf"{key}\{{\{{(.*?)\}}\}}", body, re.S).group(1)
        return np.array([complex(float(a), float(b)) for a, b in
                         re.findall(r"Sample\(([^)]+)\), Sample\(([^)]+)\)", blk)])

    def rl(key):
        blk = re.search(rf"{key}\{{\{{(.*?)\}}\}}", body, re.S).group(1)
        return np.array([float(x) for x in re.findall(r"Sample\(([^)]+)\)", blk)])

    fp = float(re.search(r"fpass = Sample\(([^)]+)\)", body).group(1))
    return dict(fpass=fp, cpx_p=cf("complexPoles"), real_p=rl("realPoles"),
                cpx_c=cf("complexCoeffsDirect"), real_c=rl("realCoeffsDirect"))


COEF = load_best()


def scaled(fs_in, fs_out):
    """SetCutoffByFpass(min(f_in,f_out)/2, fs_in) 的频率缩放。"""
    scale = (2.0 * np.pi * (min(fs_in, fs_out) / 2) / fs_in) / COEF["fpass"]
    lam = np.concatenate([COEF["cpx_p"], COEF["real_p"]]) * scale
    cc = np.concatenate([COEF["cpx_c"], COEF["real_c"]]) * scale
    return scale, lam, cc


def h_response(fs_in, fs_out, f_hz):
    """设计方程给出的连续原型响应（式 (2)），用于叠加显示。"""
    scale, _, _ = scaled(fs_in, fs_out)
    s = (1j * (2 * np.pi * np.asarray(f_hz, dtype=float) / fs_in) / scale)[:, None]
    cp, cc = COEF["cpx_p"][None, :], COEF["cpx_c"][None, :]
    h = np.sum(0.5 * (cc / (s - cp) + np.conj(cc) / (s - np.conj(cp))), axis=1)
    h = h + np.sum(COEF["real_c"][None, :] / (s - COEF["real_p"][None, :]), axis=1)
    return h


def resample(x, inc, lam, cc, d):
    """复刻 ResampleIIR::Process 的相位累积；同时给出 d≠0 直接项变体。

    状态 X_m 累积 c_m·x[k]·p̄_m^{n−k}，输出取 Σ_m Re(X_m·p̄_m^d)（系数只在注入时乘一次，
    与 elliptic_blep 的 `Add`/`Get` 一致；这里用精确 exp 代替 partial LUT）。
    """
    pbar = np.exp(lam)
    state = cc * x[0]                       # blep_.Add(x[0])
    t_all, y, y_d, hit = [], [], [], []
    rpos, phase = 0, 0.0
    while rpos < len(x) - 1:
        t = rpos + phase
        val = float(np.sum((state * np.exp(phase * lam)).real))
        h = (t == np.floor(t))
        t_all.append(t); y.append(val)
        y_d.append(val + (d * x[min(int(t), len(x) - 1)] if h else 0.0))
        hit.append(h)

        phase += inc
        new_rpos = min(rpos + int(np.floor(phase)), len(x) - 1)
        phase -= np.floor(phase)
        for k in range(rpos + 1, new_rpos + 1):
            state = state * pbar + cc * x[k]
        rpos = new_rpos
    return (np.array(t_all), np.array(y), np.array(y_d), np.array(hit), lam)


def spec(sig, fs, skip=0):
    s = sig[skip:]
    w = np.hanning(len(s))
    sp = np.abs(np.fft.rfft(s * w)) * 2.0 / w.sum()
    return np.fft.rfftfreq(len(s), 1.0 / fs), sp


def figure(fs_in, fs_out, tones):
    inc = fs_in / fs_out
    n = np.arange(N_IN)
    x = sum(np.sin(2 * np.pi * f * n / fs_in) for f, _ in tones)

    scale, lam, cc = scaled(fs_in, fs_out)
    t, y, y_d, hit, _ = resample(x, inc, lam, cc, D)

    fsig = np.array([f for f, _ in tones])
    fi, si = spec(x, fs_in)
    fo, so = spec(y, fs_out, SKIP)
    _, sd = spec(y_d, fs_out, SKIP)

    fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.3))

    # (a) 输入频谱 + 原型响应
    ax[0].plot(fi / 1e3, 20 * np.log10(np.maximum(si, 1e-12)), lw=1.0, color="C0", label="input")
    fg = np.linspace(0, fs_in / 2, 2000)
    ax[0].plot(fg / 1e3, 20 * np.log10(np.maximum(np.abs(h_response(fs_in, fs_out, fg)), 1e-12)),
               ls="--", lw=1.6, color="tab:green", label="|H(jω)| (prototype)")
    ax[0].axvline(min(fs_in, fs_out) / 2 / 1e3, color="r", ls=":", lw=1.0)
    ax[0].text(min(fs_in, fs_out) / 2 / 1e3, 2, f" min(f)/2\n {min(fs_in,fs_out)/2/1e3:.2f}k",
               color="r", fontsize=8)
    ax[0].set_xlim(0, fs_in / 2 / 1e3)
    ax[0].set_ylim(-120, 6)
    ax[0].set_xlabel("frequency (kHz) @ input rate")
    ax[0].set_ylabel("dB (rel. unit tone)")
    ax[0].set_title("(a) input spectrum  +  filter to be applied")
    ax[0].legend(fontsize=8, loc="lower left")
    ax[0].grid(alpha=0.3)

    # (b) 输出频谱
    ax[1].plot(fo / 1e3, 20 * np.log10(np.maximum(so, 1e-12)), lw=1.0, color="C0",
               label="output, d = 0 (strictly proper)")
    ax[1].plot(fo / 1e3, 20 * np.log10(np.maximum(sd, 1e-12)), lw=1.0, color="C3", alpha=0.85,
               label=f"output, d = {D} (direct term added)")
    if fs_out > fs_in:
        for f, _ in tones:
            ax[1].axvline((fs_in - f) / 1e3, color="0.5", ls=":", lw=0.8)
    ax[1].set_xlim(0, fs_out / 2 / 1e3)
    ax[1].set_ylim(-120, 6)
    ax[1].set_xlabel("frequency (kHz) @ output rate")
    ax[1].set_ylabel("dB (rel. unit tone)")
    ax[1].set_title("(b) output spectrum")
    ax[1].legend(fontsize=8, loc="lower left")
    ax[1].grid(alpha=0.3)

    # (c) 时域（截取约 1 ms；若有命中则对准第一个命中）
    span = int(round(1e-3 * fs_out))
    n0 = max(int(np.argmax(hit)) - span // 4, 0) if hit.any() else 0
    idx = np.arange(n0, min(n0 + span, len(t)))
    tt = t[idx] / fs_in * 1e3
    # 期望稳态输出：每个音乘 |H(f_i)| 与相移 ∠H(f_i)
    Ht = h_response(fs_in, fs_out, np.array([f for f, _ in tones]))
    ref = sum(np.abs(Ht[i]) * np.sin(2 * np.pi * f * t[idx] / fs_in + np.angle(Ht[i]))
              for i, (f, _) in enumerate(tones))
    ax[2].plot(tt, ref, lw=1.4, color="0.6", label="expected output  Σ|H(f_i)|·sin(…+∠H)")
    ax[2].plot(tt, y[idx], ".-", ms=3, lw=0.8, color="C0", label="output, d = 0")
    ax[2].plot(tt[hit[idx]], y_d[idx][hit[idx]], "x", ms=5, mew=1.1, color="C3",
               label="output, d ≠ 0 (only differs at hits)")
    ax[2].set_xlabel("time (ms)")
    ax[2].set_ylabel("amplitude")
    ax[2].set_title(f"(c) time domain:  hits in window = {int(hit[idx].sum())}")
    ax[2].legend(fontsize=8, loc="lower right")
    ax[2].grid(alpha=0.3)

    fig.suptitle(f"resampling {fs_in/1e3:.1f}k → {fs_out/1e3:.1f}k   (inc = {inc:.6f})   "
                 f"tones = {[f/1e3 for f, _ in tones]} kHz   d = {D}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    OUT_DIR.mkdir(exist_ok=True)
    out = OUT_DIR / f"spectrum_{int(fs_in/1000)}k_to_{int(fs_out/1000)}k.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"{out}  hits={int(hit.sum())}/{len(hit)}  "
          f"peak(d=0)={20*np.log10(np.maximum(si, 1e-12)).max():.1f} dB in")


if __name__ == "__main__":
    figure(96000, 48000, [(10000.0, 1.0), (30000.0, 1.0)])
    figure(48000, 96000, [(10000.0, 1.0), (20000.0, 1.0)])
    figure(48000, 44100, [(10000.0, 1.0), (23000.0, 1.0)])
