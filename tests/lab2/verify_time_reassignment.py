# -*- coding: utf-8 -*-
"""
verify_time_reassignment.py
===========================

【历史实验】验证 windowless NC 的「时间重分配」对低频瞬态是否真的聚焦。

背景
----
曾为 windowless NC 设计了时间重分配(TimeReassignment)版本，用左右两 DFT 分量
交叉相位算群延迟：
    cross = X_R · conj(X_L),  gd = 0.5 - arg(cross)/(2π)
然后把能量按 gd 折算成列偏移写入图像缓冲。

本脚本做**可读图**的对照实验（人 + AI 都能看图判断）：
  用一颗清晰短促的低频鼓击(60Hz, 阻尼包络)作为输入，
  分别算：
    (a) plain      —— 无时间重分配(直接输出每帧 NC gain)
    (b) reassign   —— 时间重分配(群延迟 → 整数列偏移，写入持久缓冲)
  输出对比图，判断低频瞬态是否被聚焦。

结论
----
实验显示 reassign 不聚焦，反而把连续脉冲块撕成竖条、低频出现空隙、高频冒伪影，
因此该方案**已放弃**——C++ 侧的 WindowlessNcTimeFrame 已被移除，本脚本仅作为
证据保留。

用法
----
    python verify_time_reassignment.py
输出:
    output/tr_plain_vs_reassign.png   (上=plain, 下=reassign)
"""
from __future__ import annotations

import os

import numpy as np
import soundfile as sf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# 参数(与 C++ WindowlessNcFrame/TimeFrame 一致)
# ------------------------------------------------------------
SAMPLE_RATE = 48000
FFT_SIZE = 4096
HOP = 256
F_MIN = 20.0
F_MAX = 20000.0
DB_FLOOR = -72.0
MAX_WINDOW_S = 0.125
OUT_DIR = os.path.join(os.path.dirname(__file__), "output")

# 群延迟列偏移量(与 C++ colOffset_ = fftSize/hopSize 一致)
COL_OFFSET = FFT_SIZE // HOP


# ------------------------------------------------------------
# 网格基 bin 布局(与 C++ WindowlessNcFrame::Init 一致)
# ------------------------------------------------------------
def build_grid_bins(n_bins: int, bw_scale: float = 1.0) -> list[tuple[float, float, float, int]]:
    """返回 [(f_center, f_left, f_right, N), ...]。

    约定 y=0 为最低频(20Hz)，y 增大频率升高，与 imshow(origin='lower') 的
    底部=低频轴一致。n_bins = 每 y 像素一个 bin。
    """
    log_min = np.log10(F_MIN)
    log_max = np.log10(F_MAX)
    log_step = (log_max - log_min) / n_bins
    centers = 10.0 ** (log_min + (np.arange(n_bins) + 0.5) * log_step)
    max_window = int(MAX_WINDOW_S * SAMPLE_RATE)

    bins = []
    for y in range(n_bins):
        fc = float(centers[y])
        f_hi = float(centers[y + 1]) if y + 1 < n_bins else fc
        f_lo = float(centers[y - 1]) if y > 0 else fc
        # 带宽 = 相邻 bin 中心频率差(y 增大频率升高, f_hi=上方=更高频)
        w_nc = max(abs(f_hi - f_lo) * bw_scale, 1e-3)
        q = round(2.0 * fc / w_nc)
        n_f = round(q * SAMPLE_RATE / (2.0 * fc))
        N = int(max(8, min(n_f, max_window)))
        f_left = fc - SAMPLE_RATE / (2.0 * N)
        f_right = fc + SAMPLE_RATE / (2.0 * N)
        bins.append((fc, f_left, f_right, N))
    return bins


# ------------------------------------------------------------
# 单 bin 滑动 DFT + NC 增益 + (可选)群延迟
# ------------------------------------------------------------
def _sliding_nc(x: np.ndarray, f_left: float, f_right: float, N: int) -> tuple[np.ndarray, np.ndarray]:
    """逐样本滑动 DFT，返回 (NC gain 序列, 群延迟序列)。"""
    PI = np.pi
    Wl = complex(np.cos(-2 * PI * f_left / SAMPLE_RATE), np.sin(-2 * PI * f_left / SAMPLE_RATE))
    Wr = complex(np.cos(-2 * PI * f_right / SAMPLE_RATE), np.sin(-2 * PI * f_right / SAMPLE_RATE))
    WNl = complex(np.cos(-2 * PI * f_left * N / SAMPLE_RATE), np.sin(-2 * PI * f_left * N / SAMPLE_RATE))
    WNr = complex(np.cos(-2 * PI * f_right * N / SAMPLE_RATE), np.sin(-2 * PI * f_right * N / SAMPLE_RATE))
    acl = acr = 0j
    ring = np.zeros(N)
    n = 0
    gains = np.zeros(len(x))
    gds = np.zeros(len(x))
    for i, s in enumerate(x):
        old = ring[n % N] if n >= N else 0.0
        acl = Wl * acl + s - old * WNl
        acr = Wr * acr + s - old * WNr
        ring[n % N] = s
        n += 1
        nc = -(acl.real * acr.real + acl.imag * acr.imag)
        gains[i] = np.sqrt(max(nc, 0.0) + 1e-18) / N if nc > 0 else 0.0
        cross = acr * np.conj(acl)
        arg = np.angle(cross) / (2 * PI)
        arg %= 1.0
        gds[i] = 0.5 - arg
    return gains, gds


# ------------------------------------------------------------
# 按帧抽取(NC gain / 群延迟 在 hop 边界采样)
# ------------------------------------------------------------
def frame_sample(x: np.ndarray, fc: float, f_left: float, f_right: float, N: int):
    gains, gds = _sliding_nc(x, f_left, f_right, N)
    # 只取 hop 边界(每 FFT_SIZE 前进 HOP, 与 SpectrogramColumn 一致)
    # 简化: 每 hop 采一个, 从 FFT_SIZE 起(跳过预卷)
    starts = np.arange(FFT_SIZE - 1, len(x), HOP)
    return gains[starts], gds[starts]


# ------------------------------------------------------------
# 实际实现(逐 bin 预计算, 更稳)
# ------------------------------------------------------------
def compute_matrices(x: np.ndarray, n_bins: int, use_reassign: bool, imageW: int):
    """返回 (times, freqs_or_rows, matrix_2d_sampled)。

    matrix_2d_sampled: (n_bins, n_frames) gain; use_reassign 时已做偏移聚焦。
    """
    bins = build_grid_bins(n_bins)
    frames = np.arange(FFT_SIZE - 1, len(x), HOP)
    # 缓存每 bin 的 gain/gd 在 hop 边界
    bin_gains = []
    bin_gds = []
    for fc, fl, fr, N in bins:
        g, d = frame_sample(x, fc, fl, fr, N)
        bin_gains.append(g)
        bin_gds.append(d)

    n_frames = len(frames)
    if not use_reassign:
        return frames, np.array(bin_gains)

    # 时间重分配: 持久列缓冲(宽 imageW), 头指针前进, 写入 head+off
    buf = np.zeros((n_bins, imageW))
    out = np.zeros((n_bins, n_frames))
    head = 0
    for fi in range(n_frames):
        for yi in range(n_bins):
            g = bin_gains[yi][fi]
            gd = bin_gds[yi][fi]
            off = int(round(gd * COL_OFFSET))
            off = max(-COL_OFFSET, min(COL_OFFSET, off))
            pos = (head + off) % imageW
            buf[yi][pos] = max(buf[yi][pos], g)
        out[:, fi] = buf[:, head % imageW].copy()
        buf[:, head % imageW] = 0.0
        head += 1
    return frames, out


# ------------------------------------------------------------
# 主流程
# ------------------------------------------------------------
def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── 构造低频瞬态: 60Hz 阻尼鼓击, 100ms, 出现在 0.3s ──
    dur = 0.8
    t = np.arange(int(SAMPLE_RATE * dur)) / SAMPLE_RATE
    x = np.zeros_like(t)
    on = int(0.3 * SAMPLE_RATE)
    L = int(0.1 * SAMPLE_RATE)
    env = np.exp(-np.arange(L) / (SAMPLE_RATE * 0.02))
    x[on:on + L] = 0.5 * env * np.sin(2 * np.pi * 60 * t[on:on + L])

    n_bins = 90  # 用较小的网格高度便于读图
    imageW = 100

    frames, m_plain = compute_matrices(x, n_bins, False, imageW)
    frames, m_re = compute_matrices(x, n_bins, True, imageW)

    # 归一化到 0dB
    def to_db(m):
        mx = m.max() if m.max() > 0 else 1.0
        return np.clip(20 * np.log10(np.maximum(m, 1e-9) / mx), DB_FLOOR, 0.0)

    p_db = to_db(m_plain)
    r_db = to_db(m_re)

    # ── 输出对比图(叠加时间、频率轴) ──
    fig, axes = plt.subplots(2, 1, figsize=(13, 8))
    for ax, data, title in [
        (axes[0], p_db, "plain (no time reassignment)"),
        (axes[1], r_db, "time reassignment (group-delay offset)"),
    ]:
        im = ax.imshow(data, aspect="auto", origin="lower", cmap="plasma",
                       extent=[frames[0] / SAMPLE_RATE, frames[-1] / SAMPLE_RATE, F_MIN, F_MAX])
        ax.set_yscale("log")
        ax.set_ylim(F_MIN, F_MAX)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel("Frequency (Hz)")
        ax.axvline(0.3, color="cyan", ls="--", lw=1, alpha=0.7, label="onset")
    axes[0].legend()
    axes[1].set_xlabel("Time (s)")
    fig.colorbar(im, ax=axes, label="dB (normalized)")
    fig.suptitle("Low-freq transient (60Hz, 100ms): plain vs time-reassignment", fontsize=13)
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "tr_plain_vs_reassign.png")
    fig.savefig(out, dpi=130)
    print("saved:", out)

    # 聚焦度指标
    def spread(m):
        col_peak = m.max(axis=0)          # 每帧该行最大 gain
        col_peak = col_peak / col_peak.max() if col_peak.max() > 0 else col_peak
        return int(np.sum(col_peak > 0.1))   # 超过 10% 峰值的帧数
    print(f"plain   >10% fraction frames: {spread(m_plain)}")
    print(f"reassign >10% fraction frames: {spread(m_re)}")


if __name__ == "__main__":
    main()
