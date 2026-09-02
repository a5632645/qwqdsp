# -*- coding: utf-8 -*-
"""
run_sliding_dft.py
==================

用 **递归滑动 DFT**(recursive sliding DFT) 复现论文《Window Function-less DFT
with Reduced Noise and Latency for Real-Time Music Analysis》
(arXiv:2410.07982v3) 的 NC 方法，但采用改进版实现：

改进点
------
原论文(IV-B/IV-C 节)用"全局累积相位"的参考正弦波，每样本把相位 θ_L/θ_R
累加，因此输出前必须用旋转矩阵 R(θ) 校正(式8)。本文件改用递归滑动 DFT：
每个 bin 的频谱 $X(n)$ 通过恒定旋转因子 $W=e^{-j2\pi f/F_S}$ 一步递归：

    X(n) = W·X(n-1) + x(n) − x(n−N)·W^N

因为 X(n) 恒等于"以当前窗口起点为相位零点的窗口 DFT"，相位不随时间漂移，
所以**任意时刻读出的 X 已经正确，无需任何相位校正**。

与 nc_dft.py 的关系
-------------------
nc_dft.py 用 FFT 卷积等价实现累加器(逐样本等价，同样能复现 NC 行为)；
本文件用递归滑动 DFT，逐样本流式处理，避免卷积开销，且省去式(8)相位校正。

NC 部分(式3/8、归一化、时间抽取、IIR 平滑)与 nc_dft.py 完全一致。
"""
from __future__ import annotations

import os
import sys

import numpy as np
import soundfile as sf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from nc_dft import NCDftConfig, build_bins, downsample_and_smooth  # 复用 bin/配置/抽取平滑
from run import stft_spectrogram, to_db, F_MIN, F_MAX, DB_FLOOR  # 复用 STFT 参照与 dB 工具

# 默认输入/输出
INPUT_DIR = os.path.join(os.path.dirname(__file__),
                         "..", "work_dir", "input")
OUT = os.path.join(os.path.dirname(__file__),
                   "output", "drumloop_sliding.png")


# ------------------------------------------------------------
# 递归滑动 DFT：逐样本流式计算 X(n) = W·X(n-1) + x[n] - x[n-N]·W^N
# ------------------------------------------------------------
def sliding_spectrum(x: np.ndarray, f: float, fs: float,
                     N: int) -> np.ndarray:
    """返回复频谱 X(n)(长度 len(x))，用于 NC 相位点积。

    X(n) = W·X(n-1) + x[n] − x[n−N]·W^N
        W  = exp(-j·2π·f/F_S) 为每样本恒定旋转因子；
        W^N = exp(-j·2π·f·N/F_S) 为跨一个窗长的相位因子。
    因此 X(n) 恒等于对最近 N 个样本 x[n−N+1..n] 的窗口 DFT
    (相位零点在窗口起点)，**不随时间累积相位**，无需相位校正。
    """
    W = np.exp(-2j * np.pi * f / fs)
    WN = np.exp(-2j * np.pi * f * N / fs)   # = W**N
    xpad = np.concatenate((np.zeros(N), x))  # 让 x[n-N] 在 n<N 时为 0
    X = np.empty(len(x), dtype=np.complex128)
    acc = 0j
    for n in range(len(x)):
        acc = W * acc + x[n] - xpad[n] * WN
        X[n] = acc
    return X


# ------------------------------------------------------------
# 全部 bin 的 NC 幅度矩阵
# ------------------------------------------------------------
def _compute_from_bins(x: np.ndarray, bins, fs: float,
                       normalize: bool) -> np.ndarray:
    """对整段信号计算所有 bin 的 NC 幅度矩阵。

    NC 定义需要左右分量的"复相位"做点积(式3/8)，所以用复频谱:
        NC = max(0, -(Re_L·Re_R + Im_L·Im_R))
            = max(0, -Re(conj(X_L) · X_R))
    返回 shape (num_bins, num_samples)。
    """
    n_bins = len(bins)
    ns = x.size
    out = np.empty((n_bins, ns), dtype=np.float64)
    for i, b in enumerate(bins):
        xl = sliding_spectrum(x, b.f_left, fs, b.N)
        xr = sliding_spectrum(x, b.f_right, fs, b.N)
        nc = np.maximum(0.0, -(xl.real * xr.real + xl.imag * xr.imag))
        if normalize:
            nc = nc / b.N
        nc[: b.N - 1] = np.nan   # 预卷区(窗口尚未填满)
        out[i] = nc
    return out


# ------------------------------------------------------------
# 抽取 + IIR 平滑：复用 nc_dft.downsample_and_smooth
# (与论文 IV-C 节一致，此处无需重复定义)
# ------------------------------------------------------------
def main(argv: list[str] | None = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    import argparse
    p = argparse.ArgumentParser(
        description="NC 无窗 DFT(递归滑动 DFT 版，无相位校正) 频谱对比图")
    p.add_argument("wav", nargs="?", default=os.path.join(INPUT_DIR, "drumloop.wav"))
    p.add_argument("--out", default=OUT)
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
    if data.ndim > 1:
        data = data.mean(axis=1)
    x = np.asarray(data, dtype=np.float64)

    cfg = NCDftConfig(sample_rate=fs, note_start=21.0, note_end=117.0,
                      bins_per_octave=24)
    bins = build_bins(cfg)

    # 用递归滑动 DFT 计算 NC 矩阵
    matrix = _compute_from_bins(x, bins, fs, normalize=True)
    grid_times, grid = downsample_and_smooth(matrix, fs)
    nc_db = to_db(grid)

    # STFT 参照
    t, fft_freqs, fft_db = stft_spectrogram(x, fs)

    cmap = "plasma"
    norm = Normalize(vmin=DB_FLOOR, vmax=0.0)
    fig, axes = plt.subplots(2, 1, figsize=(13, 9.1),
                             gridspec_kw={"hspace": 0.12})
    ext = [t[0], t[-1], F_MIN, F_MAX]
    im0 = axes[0].imshow(fft_db, aspect="auto", origin="lower", cmap=cmap,
                         norm=norm, extent=ext, interpolation="nearest")
    axes[0].set_yscale("log")
    axes[0].set_ylim(F_MIN, F_MAX)
    axes[0].set_title("Blackman-Harris FFT (top) vs Recursive sliding NC DFT (bottom)")
    axes[0].set_ylabel("Frequency (Hz)")

    bin_freqs = np.array([b.f_center for b in bins])
    im1 = axes[1].pcolormesh(grid_times, bin_freqs, nc_db,
                             shading="auto", cmap=cmap, norm=norm)
    axes[1].set_yscale("log")
    axes[1].set_ylim(F_MIN, F_MAX)
    axes[1].set_ylabel("Frequency (Hz)")
    axes[1].set_xlabel("Time (s)")
    for ax in axes:
        ax.set_xlim(0, max(t[-1], grid_times[-1]))

    cbar = fig.colorbar(im0, ax=axes, orientation="vertical",
                        fraction=0.035, pad=0.02)
    cbar.set_label("Response (dB, normalized)")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print("saved:", out)
    print(f"bins={len(bins)}  freq[{bin_freqs[0]:.1f},{bin_freqs[-1]:.1f}]"
          f"  time[{grid_times[0]:.2f},{grid_times[-1]:.2f}]")


if __name__ == "__main__":
    main()
