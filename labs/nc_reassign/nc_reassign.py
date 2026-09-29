# -*- coding: utf-8 -*-
"""
nc_reassign.py
==============

无窗 NC(邻域频谱分量合成, 即论文《Window Function-less DFT with Reduced
Noise and Latency for Real-Time Music Analysis》的 NC 方法)谱图的
**时间-频率重分配**核心库(离线、逐帧、对数频率轴)。

三个目标(与 C++ GUI 侧三个标准重分配帧一一对应)
------------------------------------------------
================  ==========================  =============================
变体              频率轴                      时间轴
================  ==========================  =============================
``freq``          瞬时频率                    分析窗中心(= 不重分配时间)
``time``          bin 中心频率                 群延迟(窗口内重新定位)
``tf``            瞬时频率                    群延迟
================  ==========================  =============================

标准参照实现
------------
加 Blackman-Harris 3 项窗的 STFT, 用与
``qwqdsp/example/gui/spectral/reassignment/tf_reassignment_frame.hpp`` 相同的
三个算子(``zeroPad = 1``, 与 GUI 的 ``kZeroPadFull`` 一致):

    X_h  = FFT(x·w)                    -> 幅度
    X_t  = FFT(x[n-1]·w[n])            -> 瞬时频率 = arg(X_h·conj(X_t))/2π·Fs
    X_pf = roll(X_h, 1), X_pf[0] = 0   -> 群延迟   = 0.5 - frac(arg(X_h·conj(X_pf))/2π)

NC 侧算子(推导与验证见 README)
------------------------------
滑动 DFT(锚定窗口起点相位, 无需相位校正):

    X(n) = W·X(n-1) + x[n] - x[n-N]·W^N,        W = exp(-j·2π·f/Fs)

与 ``nc_dft.py`` 的 FFT 卷积形式逐点等价, 即 x 与核 ``g[m] = exp(-j·2π·f·m/Fs)``
(m = 0..N-1)的线性卷积。bin 取左右分量 ``f_left = f_c ∓ Fs/(2N)``。

- NC 幅度:    ``gain = sqrt(max(0, -(Re_L·Re_R + Im_L·Im_R))) / N``
- 瞬时频率:    ``f_inst = -Fs/(2π)·arg(X[n]·conj(X[n-1]))``
  (本实现的滑动 DFT 相位随时间**递减**(对正频率), 与标准 DFT 反号, 故乘 -1;
  数值自检见 ``verify_log_chirp.py --selftest``。该估计等价于窗口**中心**处的瞬时频率。)
- 群延迟:      ``ψ = arg(X_R·conj(X_L))``;  ``δ = wrap(ψ + π)``
  静止纯音的 ψ 恒为 -π(与频率无关), 故以 π 为参考零点:

      ``<m> = N/2 - N·δ/(2π)``   (能量相对**当前样本**的延迟, 样本数)

  即能量相对**窗口中心** ``(N-1)/2`` 的偏移为 ``+N·δ/(2π)``(样本)。
  该算子是 KFN 谱导数法的两点(相邻 bin, 间距 Fs/N)有限差分形式, 与标准方法
  用 roll 求群延迟的截断阶数相同(两者都有 ``Δω·N = 2π`` 的模糊范围)。

时间参考约定(全库统一, 见 README「时间参考」一节)
--------------------------------------------------
图像列 c 对应绝对时刻 ``t_c = (c·hop + (fft_size-1)/2)/Fs``, 即第 0 帧的
**分析窗中心**。所有能量都落在它被估计出的绝对时刻上:

- 标准 STFT: 窗中心 = 帧中心 -> 偏移 ``gd·fft_size`` 样本(``gd ∈ (-0.5, 0.5]``)
- NC:        窗中心 = ``r_j - (N-1)/2`` -> 偏移 ``(fft_size-N)/2 + N·δ/(2π)`` 样本

这样"不重分配"就是"落在自己窗中心", 静止音对所有 bin 都落在同一列(不会因
各 bin 窗长不同而在图上倾斜), 而对数 chirp 的能量会落到真实穿越时刻上。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import fftconvolve

# ------------------------------------------------------------
# 基础工具
# ------------------------------------------------------------


def wrap_pi(x: np.ndarray | float) -> np.ndarray | float:
    """把相位卷绕到 (-π, π]。"""
    return np.angle(np.exp(1j * np.asarray(x)))


def blackman_harris_3term(n: int) -> np.ndarray:
    """BH 3 项周期窗(与 qwqdsp_window::BlackmanHarrisThreeTerm + Helper::Normalize 一致)。

    归一化系数 2/Σw, 使 bin 中心的单位实正弦的谱峰 = 1。
    """
    a0, a1, a2 = 0.4243801, 0.4973406, 0.0782793
    t = np.arange(n) / n
    w = a0 - a1 * np.cos(2 * np.pi * t) + a2 * np.cos(4 * np.pi * t)
    return w * (2.0 / w.sum())


# ------------------------------------------------------------
# 参数
# ------------------------------------------------------------
@dataclass
class Config:
    """谱图与重分配网格参数(默认值与 C++ GUI 的 reassignment.cpp 一致)。"""
    fs: float = 48000.0
    fft_size: int = 4096
    hop: int = 256
    n_rows: int = 310
    f_min: float = 20.0
    f_max: float = 20000.0
    db_floor: float = -72.0
    max_window_s: float = 0.075          # NC 低音窗长上限(秒, 论文 IV-A 节)
    max_periods: float | None = None     # 另一条窗长上限: N <= max_periods·Fs/f_c
                                         # (即 bin 的 Q 上限; None = 不启用)
    min_periods: float | None = None     # 窗长下限: N >= min_periods·Fs/f_c(至少 k 个周期
                                         # 的窗, 超低频用它把秒数上限顶开; None = 不启用)
    bandwidth_scale: float = 1.0         # NC bin 带宽缩放
    # 瞬时频率估计器(见 README「频率估计器」):
    #   "taper"(默认) 用 Y = X_R − X_L 的相位差求窗中心处的瞬时频率。Y 等价于对同一段
    #              样本做**正弦锥**(sin(πm/N))加窗的滑动 DFT(因为 e^{-jω_L m} − e^{-jω_R m}
    #              = −2j·sin(πm/N)·e^{-jω_c m})，旁瓣比矩形窗低约 10 dB、滚降快一倍,
    #              低频镜像/旁瓣污染因此显著变小; 再把 IF 按序列斜率一阶**外推到落点时刻**
    #              (重分配把能量搬到窗中心+群延迟处, 频率必须跟着走, 否则落点偏离真值线)
    #   "taper_center" 同上但不外推(诊断用)
    #   "weighted"     原实现: 左右分量各自求相位差后按幅度加权平均(隐含地跟着落点跑,
    #                  但估计本身脏: 低频 25 Hz 处误差是 taper 的 3 倍)
    if_estimator: str = "taper"
    subcell_scale: float = 0.25          # 子格宽 = FFT bin 宽 × 此系数

    @property
    def n_cols_sub(self) -> int:
        """一帧 = 多少个子列(= fft_size/hop)。"""
        return self.fft_size // self.hop

    @property
    def bin_hz(self) -> float:
        return self.fs / self.fft_size

    def n_frames(self, n_samples: int) -> int:
        return max(0, (n_samples - self.fft_size) // self.hop + 1)

    def frame_starts(self, n_samples: int) -> np.ndarray:
        return np.arange(self.n_frames(n_samples)) * self.hop

    def col_times(self, n_cols: int, col_shift: int = 0) -> np.ndarray:
        """列 c 的绝对时刻(秒); ``col_shift`` 为图左侧预垫列数(与 LogReassignGrid 一致)。"""
        return ((np.arange(n_cols) - col_shift) * self.hop
                + (self.fft_size - 1) / 2.0) / self.fs


# ------------------------------------------------------------
# 对数频率行
# ------------------------------------------------------------
def log_row_centers(cfg: Config) -> np.ndarray:
    """每行中心频率(Hz), 行 0 = 最低频, 行 ``n_rows-1`` = 最高频。

    像素中心铺在 log 轴的等分点上: ``f_y = 10^(log_min + y/(n_rows-1)·span)``。
    这是全库唯一的行-频率约定(作图、指标、网格内部几何都用它)。
    """
    span = np.log10(cfg.f_max) - np.log10(cfg.f_min)
    y = np.arange(cfg.n_rows)
    return 10.0 ** (np.log10(cfg.f_min) + span * y / max(cfg.n_rows - 1, 1))


def log_row_bands(cfg: Config) -> tuple[np.ndarray, np.ndarray]:
    """每行覆盖的频率区间 [下界, 上界): 相邻行中心的几何中点, 两端只取半格。"""
    centers = log_row_centers(cfg)
    logc = np.log10(centers)
    step = (logc[-1] - logc[0]) / max(cfg.n_rows - 1, 1)
    lo = 10.0 ** (logc - 0.5 * step)
    hi = 10.0 ** (logc + 0.5 * step)
    lo[0] = centers[0]
    hi[-1] = centers[-1]
    return lo, hi


# ------------------------------------------------------------
# NC bin
# ------------------------------------------------------------
@dataclass
class NcBin:
    """一个 NC bin: 中心频率、窗长、左右分量频率。"""
    index: int
    f_center: float
    f_left: float
    f_right: float
    N: int


def build_nc_bins(cfg: Config) -> list[NcBin]:
    """在 log 行网格上铺 NC bin(每行一个, 中心 = 行中心), 与 C++ WindowlessNcFrame 一致。

    公式(论文 arXiv:2410.07982v3):
      (6) 带宽   W_NC = f(i+1) - f(i-1)(这里用相邻行频率差)
      (7) 窗长   N    = round( round(2·f_c/W_NC)·Fs/(2·f_c) ), 夹到 [8, 上限]
      (5) 分量   f_left/right = f_c ∓ Fs/(2N)

    两条上限(见 README「窗长 N 的下限」/「低频窗长」):
      ``max_window_s``: N ≤ 固定秒数(时间分辨率/延迟上限)
      ``max_periods`` : N ≤ k·Fs/f_c(即 bin 的 Q 上限, k 个周期)
    自然值 N ≈ Fs/W_NC 相当于每个 bin 恒为 44 个周期(Q = 44, 与 log 行带宽匹配)。
    """
    centers = log_row_centers(cfg)
    max_n = max(8, int(cfg.max_window_s * cfg.fs))
    bins: list[NcBin] = []
    for i, fc in enumerate(centers):
        f_hi = centers[i - 1] if i > 0 else fc
        f_lo = centers[i + 1] if i + 1 < len(centers) else fc
        w_nc = max(abs(f_hi - f_lo) * cfg.bandwidth_scale, 1e-3)
        q = round(2.0 * fc / w_nc)
        n = int(round(q * cfg.fs / (2.0 * fc)))
        cap = max_n
        if cfg.max_periods is not None:
            cap = min(cap, int(cfg.max_periods * cfg.fs / fc))
        n = int(max(8, min(n, cap)))
        if cfg.min_periods is not None:                  # 周期数下限可顶开秒数上限(见 README)
            n = min(round(q * cfg.fs / (2.0 * fc)), max(n, round(cfg.min_periods * cfg.fs / fc)))
            n = int(max(8, n))
        bins.append(NcBin(index=i, f_center=float(fc), f_left=fc - cfg.fs / (2.0 * n),
                          f_right=fc + cfg.fs / (2.0 * n), N=n))
    return bins


def sliding_dft(x: np.ndarray, f: float, N: int, fs: float) -> np.ndarray:
    """递归滑动 DFT 的 FFT 卷积等价实现。

    ``X(n) = Σ_{m=0}^{N-1} g[m]·x[n-m]``, ``g[m] = exp(-j·2π·f·m/fs)``,
    与 ``X(n) = W·X(n-1) + x[n] - x[n-N]·W^N`` 逐点等价(nc_dft.py 已验证),
    且相位锚定窗口起点、不随时间漂移。
    返回长度 len(x); 前 ``N-1`` 个样本窗口未填满(预卷区)。
    """
    g = np.exp(-2j * np.pi * f * np.arange(N) / fs)
    return fftconvolve(x, g)[: x.size]


# ------------------------------------------------------------
# 逐 (bin, frame) 的重分配输入表
# ------------------------------------------------------------
@dataclass
class FrameTable:
    """重分配输入: 每个 (bin, frame) 的 目标频率 / 目标列 / 幅度。"""
    freq_hz: np.ndarray      # (n_bins, n_frames)
    col: np.ndarray          # (n_bins, n_frames) 列(相对第 0 帧参考, 可为小数)
    gain: np.ndarray         # (n_bins, n_frames) 线性幅度
    kind: str                # "nc" / "stft"
    variant: str             # "plain" / "freq" / "time" / "tf"
    labels: np.ndarray | None = None   # (n_bins,) 纵轴刻度用频率


def _check_variant(variant: str) -> None:
    if variant not in ("plain", "freq", "time", "tf"):
        raise ValueError(f"未知变体: {variant}")


@dataclass
class NcSeries:
    """无窗 NC 的逐 (bin, frame) 观测量(与变体无关, 只算一次)。"""
    gain: np.ndarray        # (n_bins, n_frames) NC 幅度
    if_hz: np.ndarray       # (n_bins, n_frames) 瞬时频率(窗口中心处)
    rate: np.ndarray        # (n_bins, n_frames) 瞬时频率对样本的斜率(Hz/样本)
    dt: np.ndarray          # (n_bins, n_frames) 相对窗口中心的时间偏移(样本)
    extrapolate: bool       # 是否把频率外推到落点时刻(见 if_estimator)
    ref_off: np.ndarray     # (n_bins,) 窗口中心相对帧中心的偏移(样本) = (fft_size-N)/2
    labels: np.ndarray      # (n_bins,) 中心频率
    n_frames: int


@dataclass
class StftSeries:
    """标准加窗 STFT 的逐 (bin, frame) 观测量。"""
    mag: np.ndarray         # (n_bins, n_frames) |X_h|
    if_hz: np.ndarray       # 瞬时频率(Hz, 折叠到 [0, Fs))
    gd: np.ndarray          # 群延迟 ∈ (-0.5, 0.5], 单位 = fft_size 样本
    labels: np.ndarray      # FFT bin 中心频率
    n_frames: int


def nc_series(x: np.ndarray, cfg: Config, bins: list[NcBin]) -> NcSeries:
    """算无窗 NC 的全部观测量(每个 bin 两次 FFT 卷积)。

    - 幅度      ``gain = sqrt(max(0, -(Re_L·Re_R + Im_L·Im_R)))/N``
    - 瞬时频率  ``Config.if_estimator`` 选: 默认 ``"taper"`` 用
      ``Y = X_R − X_L`` 的相位差(等价于正弦锥加窗, 旁瓣污染小), 另有 ``"weighted"``
      (左右分量相位差按幅度加权) 与 ``"circ"``(幅度² 复数加权)
    - 群延迟    ``dt = N·δ/(2π)`` 样本, ``δ = wrap(arg(X_R·conj(X_L)) + π)``;
      dt 是**能量相对窗口中心**的偏移(静止纯音 dt = 0, ψ = -π)
    """
    fs, hop = cfg.fs, cfg.hop
    n_frames = cfg.n_frames(x.size)
    idx = cfg.frame_starts(x.size) + cfg.fft_size - 1
    idx_prev = np.maximum(idx - 1, 0)
    n_bins = len(bins)
    gain = np.zeros((n_bins, n_frames))
    if_hz = np.zeros((n_bins, n_frames))
    rate = np.zeros((n_bins, n_frames))
    dt = np.zeros((n_bins, n_frames))
    ref_off = np.zeros(n_bins)

    for i, b in enumerate(bins):
        xl = sliding_dft(x, b.f_left, b.N, fs)
        xr = sliding_dft(x, b.f_right, b.N, fs)
        s = -(xl.real * xr.real + xl.imag * xr.imag)
        gain[i] = np.sqrt(np.maximum(s[idx], 0.0)) / b.N
        if cfg.if_estimator in ("taper", "taper_center"):
            # Y = X_R − X_L: 正弦锥加窗的滑动 DFT(见 Config.if_estimator 注释)
            y = xr - xl
            if_hz[i] = -fs / (2 * np.pi) * np.angle(y[idx] * np.conj(y[idx_prev]))
        else:  # "weighted"
            wl, wr = np.abs(xl[idx]), np.abs(xr[idx])
            il = -fs / (2 * np.pi) * np.angle(xl[idx] * np.conj(xl[idx_prev]))
            ir = -fs / (2 * np.pi) * np.angle(xr[idx] * np.conj(xr[idx_prev]))
            if_hz[i] = (wl * il + wr * ir) / np.maximum(wl + wr, 1e-300)
        dt[i] = b.N * wrap_pi(np.angle(xr[idx] * np.conj(xl[idx])) + np.pi) / (2 * np.pi)
        ref_off[i] = (cfg.fft_size - b.N) / 2.0
        bad = idx < b.N - 1                       # 滑动 DFT 未填满(fs 内 fft_size >= N 时不出现)
        gain[i][bad] = 0.0
        if_hz[i][bad] = 0.0
        dt[i][bad] = 0.0
        rate[i] = np.gradient(if_hz[i]) / hop     # IF 序列斜率(Hz/样本), 供外推用

    return NcSeries(gain=gain, if_hz=if_hz, rate=rate, dt=dt, ref_off=ref_off,
                    extrapolate=(cfg.if_estimator == "taper"),
                    labels=np.array([b.f_center for b in bins]), n_frames=n_frames)


def stft_series(x: np.ndarray, cfg: Config) -> StftSeries:
    """算标准加窗 STFT 的全部观测量(zeroPad = 1, 与 GUI 的 kZeroPadFull 一致)。

      X_h  = FFT(x·w) / X_t = FFT(x[n-1]·w[n]) / X_pf = roll(X_h, 1), X_pf[0] = 0
    """
    fs = cfg.fs
    w = blackman_harris_3term(cfg.fft_size)
    starts = cfg.frame_starts(x.size)
    n_frames = len(starts)
    fft_len = cfg.fft_size
    n_bins = fft_len // 2 + 1
    freqs_bin = np.arange(n_bins) * fs / fft_len

    frames = np.stack([x[s: s + cfg.fft_size] for s in starts])
    shifted = np.concatenate([np.zeros((n_frames, 1)), frames[:, :-1]], axis=1)
    Xh = np.fft.rfft(frames * w, fft_len, axis=1)
    Xt = np.fft.rfft(shifted * w, fft_len, axis=1)
    Xpf = np.roll(Xh, 1, axis=1)
    Xpf[:, 0] = 0.0

    if_norm = np.angle(Xh * np.conj(Xt)).T / (2 * np.pi)
    if_norm -= np.floor(if_norm)
    arg_f = np.angle(Xh * np.conj(Xpf)).T / (2 * np.pi)
    gd = 0.5 - (arg_f - np.floor(arg_f))
    mag = np.abs(Xh).T
    keep = mag > 1e-8
    keep[0] = False                       # 与 C++ 一致: X_pf[0] = 0, 第 0 bin 无群延迟可用
    return StftSeries(mag=np.where(keep, mag, 0.0), if_hz=if_norm * fs, gd=gd,
                      labels=freqs_bin, n_frames=n_frames)


def table_from_series(series: NcSeries | StftSeries, cfg: Config, variant: str) -> FrameTable:
    """观测量 → 某个变体的重分配输入表。

    ``freq``/``tf`` 用瞬时频率, 否则用 bin 中心; ``time``/``tf`` 用群延迟, 否则落在窗中心。
    """
    _check_variant(variant)
    n_frames = series.n_frames
    j = np.arange(n_frames)[None, :]
    # 各 bin 的「窗中心相对帧中心」偏移(样本): NC 为 (fft_size-N)/2, STFT 为 0
    base = (series.ref_off[:, None] if isinstance(series, NcSeries)
            else np.zeros((len(series.labels), 1)))
    if isinstance(series, NcSeries):
        freq = series.if_hz if variant in ("freq", "tf") else np.repeat(series.labels[:, None], n_frames, 1)
        off = series.dt if variant in ("time", "tf") else np.zeros_like(series.dt)
        if series.extrapolate and variant in ("freq", "tf"):
            # 重分配把能量搬到"窗中心 + off"处, 频率必须同步外推, 否则落点偏离真值线:
            # freq 变体 off = 0(落在窗中心) → 自动退化为窗中心处的频率
            freq = freq + series.rate * off
        col = j + (base + off) / cfg.hop
    else:
        freq = series.if_hz if variant in ("freq", "tf") else np.repeat(series.labels[:, None], n_frames, 1)
        off = series.gd * cfg.fft_size if variant in ("time", "tf") else np.zeros_like(series.gd)
        col = j + (base + off) / cfg.hop
    keep = series.gain if isinstance(series, NcSeries) else series.mag
    valid = (freq >= cfg.f_min) & (freq <= cfg.f_max)
    return FrameTable(freq_hz=np.where(valid, freq, 0.0),
                      col=np.where(valid, col, 0.0),
                      gain=np.where(valid, keep, 0.0),
                      kind="nc" if isinstance(series, NcSeries) else "stft",
                      variant=variant, labels=series.labels)


def build_series(x: np.ndarray, cfg: Config, kind: str, bins: list[NcBin] | None = None):
    """按方法算观测量。"""
    if kind == "nc":
        return nc_series(x, cfg, build_nc_bins(cfg) if bins is None else bins)
    if kind == "stft":
        return stft_series(x, cfg)
    raise ValueError(f"未知方法: {kind}")


# ------------------------------------------------------------
# 对数重分配网格(C++ LogReassignGrid 的离线等价)
# ------------------------------------------------------------
class LogReassignGrid:
    """log 频率行 → 线性子格; 子格内求和、行内取 max; 时间按列双线性。

    每行把行带宽按 ``bin_hz·subcell_scale`` 细分成 k 个子格(低频行带宽不足一个
    子格时 k = 1), 落到 (行, 子格) 后**子格内求和、行内取 max**: 同一主瓣的
    bin 重分配后落进同一子格(能量集中), 而宽行内不同频率的 bin 分处不同子格,
    不会叠在一起抬高噪声底。全程不加权。
    """

    def __init__(self, cfg: Config, n_cols: int, col_shift: int = 0,
                 bin_hz: float | None = None):
        self.cfg = cfg
        self.col_shift = col_shift          # 左侧预垫列数(容纳"过去"的重分配落点)
        self.n_cols = n_cols + col_shift
        log_min, log_max = np.log10(cfg.f_min), np.log10(cfg.f_max)
        self.log_min, self.log_span = log_min, log_max - log_min
        subcell_hz = (cfg.bin_hz if bin_hz is None else bin_hz) * cfg.subcell_scale

        band_lo, band_hi = log_row_bands(cfg)
        self.row_f_lo = band_lo
        self.row_cell_w = np.zeros(cfg.n_rows)
        self.row_k = np.zeros(cfg.n_rows, dtype=int)
        self.buf: list[np.ndarray] = []
        for y in range(cfg.n_rows):
            row_w = band_hi[y] - band_lo[y]
            k = max(int(round(row_w / subcell_hz)), 1)
            self.row_k[y] = k
            self.row_cell_w[y] = row_w / k
            self.buf.append(np.zeros((k, self.n_cols)))

    def add(self, freq_hz: np.ndarray, col: np.ndarray, mag: np.ndarray) -> None:
        """把一批 (频率, 列, 幅度) 累加进网格(频率双线性到相邻行、时间双线性到相邻列)。"""
        keep = (mag > 0.0) & (freq_hz >= self.cfg.f_min) & (freq_hz <= self.cfg.f_max)
        if not np.any(keep):
            return
        f = freq_hz[keep]
        m = mag[keep]
        c_pos = np.clip(col[keep] + self.col_shift, 0.0, self.n_cols - 1.0)
        c_idx = np.minimum(np.floor(c_pos).astype(int), self.n_cols - 1)
        c_frac = c_pos - c_idx

        norm = (np.log10(f) - self.log_min) / self.log_span
        y_pos = np.clip((self.cfg.n_rows - 1) * norm, 0.0, self.cfg.n_rows - 1.0)
        y0 = np.floor(y_pos).astype(int)
        y_frac = y_pos - y0
        y0 = np.minimum(y0, self.cfg.n_rows - 1)
        y1 = np.minimum(y0 + 1, self.cfg.n_rows - 1)

        for row, wy in ((y0, 1.0 - y_frac), (y1, y_frac)):
            for r in np.unique(row[wy > 0.0]):
                s = (row == r) & (wy > 0.0)
                sub = np.clip(((f[s] - self.row_f_lo[r]) / self.row_cell_w[r]).astype(int), 0,
                              self.row_k[r] - 1)
                vals = m[s] * wy[s]
                np.add.at(self.buf[r], (sub, c_idx[s]), vals * (1.0 - c_frac[s]))
                nz = c_frac[s] > 0.0
                if np.any(nz):
                    np.add.at(self.buf[r], (sub[nz], c_idx[s][nz] + 1), vals[nz] * c_frac[s][nz])

    def emit(self) -> np.ndarray:
        """归约: 行内取 max → (n_rows, n_cols + col_shift) 线性幅度图。"""
        out = np.zeros((self.cfg.n_rows, self.n_cols))
        for y in range(self.cfg.n_rows):
            out[y] = self.buf[y].max(axis=0)
        return out


def to_db(img: np.ndarray, cal: float, db_floor: float) -> np.ndarray:
    """线性幅度图 → dB(以 cal 为 0 dB), 夹到 [db_floor, 0]。"""
    return np.clip(20.0 * np.log10(np.maximum(img, 1e-12) / max(cal, 1e-12)), db_floor, 0.0)


# ------------------------------------------------------------
# 上层: 一次算完一张图
# ------------------------------------------------------------
VARIANTS = ("plain", "freq", "time", "tf")


def build_table(x: np.ndarray, cfg: Config, kind: str, variant: str,
                bins: list[NcBin] | None = None) -> FrameTable:
    """按方法/变体算重分配输入表(单一变体的便捷入口; 多变体请复用 build_series)。"""
    return table_from_series(build_series(x, cfg, kind, bins), cfg, variant)


def render(table: FrameTable, cfg: Config, n_cols: int, col_shift: int = 0,
           bin_hz: float | None = None) -> np.ndarray:
    """把输入表铺进对数重分配网格, 返回线性幅度图 (n_rows, n_cols + col_shift)。"""
    grid = LogReassignGrid(cfg, n_cols, col_shift=col_shift, bin_hz=bin_hz)
    grid.add(table.freq_hz.ravel(), table.col.ravel(), table.gain.ravel())
    return grid.emit()


def calibrate(cfg: Config, kind: str, variant: str, f_cal: float = 1000.0,
              duration: float = 1.0, bins: list[NcBin] | None = None) -> float:
    """用单位幅度稳态正弦标定: 返回使该纯音显示为 0 dB 的**线性**幅度系数。

    标定频率取最靠近 ``f_cal`` 的 log 行中心(避免频率轴插值衰减)。
    """
    centers = log_row_centers(cfg)
    f0 = float(centers[np.argmin(np.abs(np.log(centers) - np.log(f_cal)))])
    n = int(duration * cfg.fs) + cfg.fft_size
    t = np.arange(n) / cfg.fs
    x = np.sin(2 * np.pi * f0 * t)
    n_cols = cfg.n_frames(n) + cfg.n_cols_sub
    table = build_table(x, cfg, kind, variant, bins)
    img_lin = render(table, cfg, n_cols, col_shift=cfg.n_cols_sub + 8)
    return float(img_lin.max())

