# -*- coding: utf-8 -*-
"""
nc_dft.py
=========

论文《Window Function-less DFT with Reduced Noise and Latency for Real-Time
Music Analysis》(arXiv:2410.07982v3) 的 Python 实现。

核心思想
--------
传统 FFT/DFT 输出线性间隔的频率区间，不适合音乐分析。本算法直接把每个输出
区间(NC bin)中心频率按十二平均律的指数公式铺设，使其对应音高；再对每个 bin
用 "邻域频谱分量合成(NC)" 方法——取左右两个相邻 DFT 分量的实/虚部做点积，
得到无窗函数、低旁瓣噪声、更高频分辨率的响应。

关键公式(与论文编号一致)
------------------------
- (1)  bin 中心频率:  f(n) = 440 Hz * 2^((n-69)/12)
- (7)  窗长:           N_k = round( round(2*f_c/W_NC) * F_S/(2*f_c) )
- (5)  左右 bin 频率:  f_left = f_c - F_S/(2N), f_right = f_c + F_S/(2N)
- (6)  NC bin 带宽:    W_NC = f_right - f_left = F_S/N
- (3)  NC 幅度:        max(0, -(Re_L*Re_R + Im_L*Im_R))
- (8)  滑动窗相位校正 + NC 幅度(输出前对累加器做旋转矩阵补正)

实现说明
--------
论文用逐样本累加器做成 FIR 滤波器，并在输出前旋转相位抵消滑动窗带来的相位
累积。数学上，每个 bin 在时刻 n 的"已旋转"累加值等价于：
    X[n] = sum_{m=0}^{N-1} x[n-m] * exp(-j*2*pi*f*m/F_S)
即信号 x 与核 g[m]=exp(-j*2*pi*f*m/F_S) 的线性卷积。这里用 FFT 卷积代替逐
样本累加，结果与论文算法等价(离线处理)。随后按 (3)/(8) 计算 NC 幅度并归一化。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.signal import fftconvolve


# ------------------------------------------------------------
# (1) 十二平均律中心频率
# ------------------------------------------------------------
def note_freq(note: float, a4: float = 440.0) -> float:
    """返回 MIDI 音号 note(n 可为小数)对应的频率(Hz)。n=69 -> A4=440Hz。"""
    return a4 * 2.0 ** ((note - 69.0) / 12.0)


# ------------------------------------------------------------
# 参数
# ------------------------------------------------------------
@dataclass
class NCDftConfig:
    """控制 NC-DFT 的 bin 布局。

    默认使用论文实验参数：8 个八度(27.5Hz~7040Hz)，每八度 24 个 NC bin
    (每音 2 个 bin)，共 192 个 bin。
    """
    sample_rate: float = 48000.0
    note_start: float = 21.0       # A0 = 27.5 Hz
    note_end: float = 117.0        # 8 个八度之上 => 21 + 96
    bins_per_octave: int = 24     # 每音高 2 个 NC bin
    max_window_s: float = 0.125    # 低音窗长上限(论文 IV-A 节)

    @property
    def num_bins(self) -> int:
        n0, n1 = self.note_start, self.note_end
        step = 12.0 / self.bins_per_octave
        return int(round((n1 - n0) / step)) + 1

    @property
    def spans_octaves(self) -> float:
        return self.note_end - self.note_start


# ------------------------------------------------------------
# bin 描述
# ------------------------------------------------------------
@dataclass
class NcBin:
    """单个 NC bin 的几何与滤波器定义。"""
    index: int
    f_center: float   # 中心频率(Hz)
    f_left: float     # 左分量 bin 频率
    f_right: float    # 右分量 bin 频率
    N: int            # 窗长(样本数)

    # 卷积核(延迟 m 处的复系数 g[m]=exp(-j*2*pi*f*m/F_S))
    g_left: np.ndarray = field(repr=False, default=None)
    g_right: np.ndarray = field(repr=False, default=None)


# ------------------------------------------------------------
# 构造 bin
# ------------------------------------------------------------
def build_bins(cfg: NCDftConfig) -> list[NcBin]:
    """按论文 公式(1)(5)(7)(6) 构造全部 NC bin。"""
    fs = cfg.sample_rate
    step = 12.0 / cfg.bins_per_octave
    notes = np.arange(cfg.note_start, cfg.note_end + 1e-9, step)
    centers = note_freq(notes)                     # (1)

    bins: list[NcBin] = []
    for i, f_c in enumerate(centers):
        # 带宽 = 相邻两 bin 中心频率之差(论文 IV 节: W_NC = f(i+1) - f(i-1))
        lo = centers[i - 1] if i > 0 else f_c
        hi = centers[i + 1] if i + 1 < len(centers) else f_c
        w_nc = hi - lo                             # (6) 的左半——带宽

        # (7) 窗长; 低音窗长受 max_window_s 限制(IV-A 节)
        q = round(2.0 * f_c / w_nc)
        n_candidate = round(q * fs / (2.0 * f_c))
        n = max(8, min(n_candidate, int(cfg.max_window_s * fs)))

        # (5) 左右 bin 频率
        f_left = f_c - fs / (2.0 * n)
        f_right = f_c + fs / (2.0 * n)

        # 构造核 g[m] = exp(-j*2*pi*f*m/F_S)，m = 0..N-1
        m = np.arange(n)
        ang = -2.0 * np.pi * m / fs
        g_l = np.exp(1j * ang * f_left)
        g_r = np.exp(1j * ang * f_right)

        bins.append(NcBin(index=i, f_center=f_c, f_left=f_left,
                          f_right=f_right, N=n, g_left=g_l, g_right=g_r))
    return bins


# ------------------------------------------------------------
# 逐 bin 卷积 -> NC 幅度
# ------------------------------------------------------------
def _bin_magnitude(x: np.ndarray, bin_: NcBin) -> np.ndarray:
    """计算单个 bin 在全部样本时刻的 NC 幅度(未归一化)。

    X_L[n] = (x * g_left)[n], X_R[n] = (x * g_right)[n]
    NC = max(0, -(Re_L*Re_R + Im_L*Im_R)) = max(0, -Re(conj(X_L) * X_R))
    返回长度 len(x) 的数组；前 N-1 个样本为卷积预卷区，用 NaN 屏蔽。
    """
    xl = fftconvolve(x, bin_.g_left)[: x.size]
    xr = fftconvolve(x, bin_.g_right)[: x.size]
    nc = np.maximum(0.0, -(xl.real * xr.real + xl.imag * xr.imag))
    nc[: bin_.N - 1] = np.nan   # 卷积边界的预卷区
    return nc


def compute_nc_matrix(x: np.ndarray, bins: list[NcBin],
                      normalize: bool = True) -> np.ndarray:
    """对整段信号计算所有 bin 的 NC 幅度矩阵。

    返回 shape = (num_bins, num_samples) 的 float 数组，行=bin(低->高音)，
    列=样本时刻。已故按各自窗长归一化(NC / N)。
    """
    n_bins = len(bins)
    ns = x.size
    out = np.empty((n_bins, ns), dtype=np.float64)
    for i, b in enumerate(bins):
        mag = _bin_magnitude(x, b)
        if normalize:
            mag = mag / b.N
        out[i] = mag
    return out


# ------------------------------------------------------------
# 时间轴抽取 + IIR 平滑
# ------------------------------------------------------------
def downsample_and_smooth(matrix: np.ndarray, sample_rate: float,
                          hop_ms: float = 5.0,
                          iir_alpha: float = 0.5) -> tuple[np.ndarray, np.ndarray]:
    """把逐样本 NC 矩阵抽到显示网格，并对每 bin 做一阶低通(IIR)平滑。

    论文 IV-C 节: 高频 bin 窗长极短(毫秒级)，若输出帧率不够会漏数据；用
    IIR 滤波器拉平各 bin 的响应时间。

    参数
    ----
    matrix    : (num_bins, num_samples)
    hop_ms    : 输出时间栅格的间隔(ms)。论文以"最小窗长一半"为推送间隔。
    iir_alpha : 一阶 IIR 系数(0,1]，越大越平滑(反应越慢)。

    返回
    ----
    (grid_times, grid_matrix)
        grid_times   : 长度 H 的数组(秒)
        grid_matrix  : (num_bins, H)
    """
    hop = max(1, int(round(hop_ms * sample_rate / 1000.0)))
    # 只取有效区(跳过最靠前的部分预卷样本)
    n_bins, ns = matrix.shape
    windowed = matrix[:, :ns]
    valid = ~np.isnan(windowed)
    windowed = np.where(valid, windowed, 0.0)

    # 抽取到栅格
    grid_cols = np.arange(ns, step=hop)
    grid = windowed[:, grid_cols]
    grid_times = grid_cols / sample_rate

    # 每 bin 一阶 IIR: y[k] = (1-a)*y[k-1] + a*x[k]
    grid = np.apply_along_axis(
        lambda row: _iir(row, iir_alpha), 1, grid)
    return grid_times, grid


def _iir(row: np.ndarray, alpha: float) -> np.ndarray:
    y = np.empty_like(row, dtype=np.float64)
    acc = 0.0
    for k, v in enumerate(row):
        acc = (1.0 - alpha) * acc + alpha * v
        y[k] = acc
    return y


# ------------------------------------------------------------
# 上层便利函数
# ------------------------------------------------------------
def analyze(x: np.ndarray, sample_rate: float, cfg: NCDftConfig | None = None):
    """一次性: 构造 bin -> 计算 NC 矩阵 -> 抽取平滑。

    返回 (bin_frequencies, grid_times, grid_matrix)。
    bin_frequencies: 各 bin 中心频率(Hz)，做纵轴刻度。
    """
    cfg = cfg or NCDftConfig(sample_rate=sample_rate)
    bins = build_bins(cfg)
    matrix = compute_nc_matrix(x, bins)
    grid_times, grid = downsample_and_smooth(matrix, sample_rate)
    freqs = np.array([b.f_center for b in bins])
    return freqs, grid_times, grid
