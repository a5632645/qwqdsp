# -*- coding: utf-8 -*-
"""
run.py
======

驱动脚本：读取 drumloop.wav，生成与论文 Fig.4 风格一致的双面板频谱图。

- 上面板：Blackman-Harris 加窗 STFT（常规 FFT 参照）
- 下面板：本文提出的 NC 无窗 DFT
- 纵轴：对数频率 27.5 ~ 3520 Hz；幅度为对数值(dB)，归一化到最大响应
- 色标：plasma，0 dB(最强) 深紫 -> -60 dB(最弱) 淡黄

用法
----
    python run.py
输出图像默认保存到 ./output/drumloop_combined.png
"""
from __future__ import annotations

import os
import sys

import numpy as np
import soundfile as sf
import matplotlib
matplotlib.use("Agg")  # 无界面后端，直接出图
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.signal import stft, get_window

from nc_dft import NCDftConfig, analyze

# ------------------------------------------------------------
# 配置
# ------------------------------------------------------------
F_MIN = 27.5      # A0
F_MAX = 3520.0    # A7 (参考图纵轴上限)
DB_FLOOR = -60.0  # 幅度下限(dB)
BINS_PER_OCT = 24

# 输入音频与输出图片路径。可用命令行参数覆盖: python run.py <name> [--out x.png]
# <name> 在 ../work_dir/input/ 下查找 (可带 .wav)。
INPUT_DIR = os.path.join(os.path.dirname(__file__),
                         "..", "work_dir", "input")
DEFAULT_WAV = os.path.join(INPUT_DIR, "drumloop.wav")
OUT = os.path.join(os.path.dirname(__file__),
                   "output", "drumloop_combined.png")


def to_db(amp: np.ndarray, ref: float | None = None) -> np.ndarray:
    """幅度转 dB，归一化到最大值(或 ref)，并夹到 [DB_FLOOR, 0]。"""
    ref = float(np.nanmax(amp)) if ref is None else float(ref)
    with np.errstate(divide="ignore"):
        db = 20.0 * np.log10(np.maximum(amp, 1e-30) / ref)
    return np.clip(db, DB_FLOOR, 0.0)


def stft_spectrogram(x: np.ndarray, fs: float):
    """Blackman-Harris 加窗 STFT，返回 (times, freqs, db)。

    只保留 [F_MIN, F_MAX] 范围内的频率；幅度为 |X|。
    """
    nperseg = 4096
    noverlap = nperseg // 2
    win = get_window("blackmanharris", nperseg)
    f, t, Z = stft(x, fs=fs, window=win, nperseg=nperseg,
                   noverlap=noverlap, boundary=None, padded=False)
    mag = np.abs(Z)
    db = to_db(mag)
    # 裁剪到目标频率范围
    keep = (f >= F_MIN) & (f <= F_MAX)
    return t, f[keep], db[keep]


def main(argv: list[str] | None = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    import argparse
    p = argparse.ArgumentParser(description="STFT vs NC-windowless-DFT spectrogram")
    p.add_argument("wav", nargs="?", default=DEFAULT_WAV,
                   help="输入 wav 路径，或 ../work_dir/input/ 下的文件名(可省略扩展名)")
    p.add_argument("--out", default=OUT, help="输出图片路径")
    args = p.parse_args(argv)

    wav = args.wav
    if not os.path.isabs(wav):
        name = wav if wav.lower().endswith(".wav") else wav + ".wav"
        cand = os.path.join(INPUT_DIR, name)
        if os.path.exists(cand):
            wav = cand
    out = args.out
    os.makedirs(os.path.dirname(out), exist_ok=True)

    data, fs = sf.read(wav)
    if data.ndim > 1:                 # 立体声取均值
        data = data.mean(axis=1)
    x = np.asarray(data, dtype=np.float64)

    # NC 算法
    cfg = NCDftConfig(sample_rate=fs, note_start=21.0, note_end=117.0,
                      bins_per_octave=BINS_PER_OCT)
    bin_freqs, grid_times, grid = analyze(x, fs, cfg)
    nc_db = to_db(grid)

    # STFT 参照
    t, fft_freqs, fft_db = stft_spectrogram(x, fs)

    # ------------------------------------------------------------
    # 绘图：上下双面板 + 右侧色标
    # ------------------------------------------------------------
    # 正向 plasma 映射:  0 dB(强) = 亮黄(顶部),  -60 dB(弱) = 深紫(底部)
    cmap = "plasma"
    norm = Normalize(vmin=DB_FLOOR, vmax=0.0)

    fig, axes = plt.subplots(2, 1, figsize=(13, 9.1), sharex=False,
                             gridspec_kw={"hspace": 0.12})
    # 上图：STFT
    ext = [t[0], t[-1], F_MIN, F_MAX]
    im0 = axes[0].imshow(fft_db, aspect="auto", origin="lower", cmap=cmap,
                         norm=norm, extent=ext, interpolation="nearest")
    axes[0].set_yscale("log")
    axes[0].set_ylim(F_MIN, F_MAX)
    axes[0].set_title("Blackman-Harris windowed FFT (top) vs. Proposed NC method (bottom)",
                      fontsize=13)
    axes[0].set_ylabel("Frequency (Hz)", fontsize=11)

    # 下图：NC
    im1 = axes[1].pcolormesh(grid_times, bin_freqs, nc_db,
                             shading="auto", cmap=cmap, norm=norm)
    axes[1].set_yscale("log")
    axes[1].set_ylim(F_MIN, F_MAX)
    axes[1].set_ylabel("Frequency (Hz)", fontsize=11)
    axes[1].set_xlabel("Time (s)", fontsize=11)

    for ax in axes:
        ax.set_xlim(0, max(t[-1], grid_times[-1]))

    # 右侧色标
    cbar = fig.colorbar(im0, ax=axes, orientation="vertical",
                        fraction=0.035, pad=0.02)
    cbar.set_label("Response (dB, normalized)", fontsize=11)

    fig.savefig(out, dpi=130, bbox_inches="tight")
    print("saved:", out)
    print(f"bins={len(bin_freqs)}  freq[{bin_freqs[0]:.1f},{bin_freqs[-1]:.1f}]"
          f"  time[{grid_times[0]:.2f},{grid_times[-1]:.2f}]")


if __name__ == "__main__":
    main()
