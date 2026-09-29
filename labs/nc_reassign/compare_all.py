# -*- coding: utf-8 -*-
"""
compare_all.py
==============

**只出图、不做任何数据分析**: 把六种方法画在同一张图里, 供人直接看图比较。

测试信号(三段叠加, 2 s @48 kHz):
  1. 稳态纯音 440 Hz
  2. 对数 chirp 22 Hz → 11.3 kHz(9 个八度)
  3. 单个 delta(冲激), 位于正中 t = 1.0 s

六种方法:
  1. ``stft-tf``     标准加窗 STFT 的时间+频率重分配(C++ ``tf_reassignment_frame`` 等价)
  2. ``nc-cap``      无窗 NC 时间+频率重分配(现 C++ 口径: ``max_window_s = 0.075 s``)
  3. ``nc-free``     无窗 NC 时间+频率重分配(理论窗长, 不截短)
  4. ``nc-floor4``   无窗 NC 时间+频率重分配(超低频 ``N ≥ 4 周期``)
  5. ``cqt``         **librosa** 的 CQT(恒定 Q 变换) —— 不做任何重分配
  6. ``nc-plain``    无窗 NC 原样显示 —— 不做任何重分配

方法来源: 重分配与无窗 NC 用本目录 ``nc_reassign.py``(与 C++ 帧逐算子对应);
CQT 用已安装的 ``librosa.cqt``(``bins_per_octave=31`` 使 Q≈44, 与显示对数网格/NC
自然设计的 Q 一致, 于是 CQT 的 310 个 bin 与显示 310 行一一对齐)。

幅度: 每种方法各自用**单位幅度稳态正弦**标定到 0 dB, 因此六个面板的 dB 直接可比。

用法
----
    python compare_all.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nc_reassign as nr          # noqa: E402

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
CMAP = "magma"
DB_MIN, DB_MAX = -60.0, 0.0

# ── 测试信号参数 ──
DURATION = 2.0
TONE_HZ = 440.0
CHIRP_F0 = 22.0
CHIRP_OCTAVES = 9.0
# delta 用单样本冲激: NC 对瞬态的增益 ≈1/N, 幅度要给足才看得见
AMP_TONE, AMP_CHIRP, AMP_DELTA = 0.5, 0.5, 6.0


def chirp_freq(t):
    """对数 chirp 的瞬时频率(Hz)。"""
    k = CHIRP_OCTAVES / DURATION
    return CHIRP_F0 * 2.0 ** (k * np.asarray(t, dtype=float))


def test_signal(cfg: nr.Config) -> tuple[np.ndarray, float]:
    """纯音 + 对数 chirp + 中间 delta 的叠加信号; 返回 (x, delta 时刻)。"""
    n = int(DURATION * cfg.fs)
    t = np.arange(n) / cfg.fs
    a = (CHIRP_OCTAVES / DURATION) * np.log(2.0)
    chirp = AMP_CHIRP * np.sin(2 * np.pi * CHIRP_F0 * (np.exp(a * t) - 1.0) / a)
    tone = AMP_TONE * np.sin(2 * np.pi * TONE_HZ * t)
    x = chirp + tone
    t_delta = DURATION / 2.0
    x[int(round(t_delta * cfg.fs))] += AMP_DELTA
    return x, t_delta


# ------------------------------------------------------------
# CQT(librosa, 不做重分配)
# ------------------------------------------------------------
def cqt_bins_per_octave(cfg: nr.Config) -> int:
    """取 bins_per_octave 使 Q = 1/(2^(1/bpo)-1) 与显示行网格的 Q 一致。

    显示 310 行覆盖 20 Hz–20 kHz = 9.966 octave → Q ≈ 44; bpo = 31 时 CQT 的 Q ≈ 44
    且 bin 数 = 310, 与显示行一一对应。
    """
    return int(round(cfg.n_rows / np.log2(cfg.f_max / cfg.f_min)))


def cqt_magnitude(x: np.ndarray, cfg: nr.Config, bpo: int, hop: int) -> np.ndarray:
    """librosa CQT 幅度, 返回 (n_bins, n_frames); 无相位、无重分配。"""
    import librosa
    n_bins = int(np.ceil(np.log2(cfg.f_max / cfg.f_min) * bpo)) + 1
    C = librosa.cqt(x.astype(np.float32), sr=int(cfg.fs), hop_length=hop,
                    fmin=float(cfg.f_min), n_bins=n_bins, bins_per_octave=bpo,
                    filter_scale=1.0)
    return np.abs(C)


def cqt_image(x: np.ndarray, cfg: nr.Config, n_cols: int, col_shift: int,
              bpo: int, cal: float) -> np.ndarray:
    """CQT → 与其它方法同一张显示网格。

    - 频率: CQT 的 bin k 与显示行 k 一一对齐(Q 相符, 都从 f_min 起)
    - 时间: librosa 第 j 帧中心在 j·hop 样本; 本目录第 c 列对应样本
      (c-col_shift)·hop + (fft_size-1)/2 → 整块右移 col_shift - (fft_size-1)/(2·hop) 列
    """
    mag = cqt_magnitude(x, cfg, bpo, cfg.hop)
    n_bins, n_frames = mag.shape
    img = np.zeros((cfg.n_rows, n_cols + col_shift))
    delay = int(round((cfg.fft_size - 1) / (2 * cfg.hop)))
    off = col_shift - delay
    rows = min(n_bins, cfg.n_rows)
    cols = max(0, min(n_frames, img.shape[1] - off))
    img[:rows, off:off + cols] = mag[:rows, :cols]
    return nr.to_db(img, cal, DB_MIN)


# ------------------------------------------------------------
# 六种方法的图像
# ------------------------------------------------------------
def build_images(x: np.ndarray) -> dict:
    """返回 {名称: (img_db, cfg, col_shift, 备注)}。"""
    cfg0 = nr.Config()                                   # 0.075 s 上限(现 C++ 口径)
    cfg_free = nr.Config(max_window_s=1e3)               # 理论窗长
    cfg_fl4 = nr.Config(min_periods=4)                   # 超低频 ≥4 周期
    out = {}

    def add_nc(name, cfg, variant, note):
        bins = nr.build_nc_bins(cfg)
        n_frames = cfg.n_frames(x.size)
        n_max = max(b.N for b in bins)
        pad = cfg.n_cols_sub + 8 + int(np.ceil(n_max / (2 * cfg.hop))) + 8
        pad_r = int(np.ceil(n_max / (2 * cfg.hop))) + 8
        cal = nr.calibrate(cfg, "nc", variant, duration=0.5, bins=bins)
        table = nr.table_from_series(nr.build_series(x, cfg, "nc", bins), cfg, variant)
        img = nr.render(table, cfg, n_frames + pad_r, pad)
        out[name] = (nr.to_db(img, cal, DB_MIN), cfg, pad, note)

    def add_stft(name, variant, note):
        cfg = cfg0
        n_frames = cfg.n_frames(x.size)
        pad = cfg.n_cols_sub + 8
        cal = nr.calibrate(cfg, "stft", variant)
        table = nr.table_from_series(nr.build_series(x, cfg, "stft"), cfg, variant)
        img = nr.render(table, cfg, n_frames, pad)
        out[name] = (nr.to_db(img, cal, DB_MIN), cfg, pad, note)

    add_stft("stft-tf", "tf", "STFT tf reassign (C++ tf_reassignment_frame equivalent)")
    add_nc("nc-cap", cfg0, "tf", "NC tf, max_window_s = 0.075 s (C++ cap)")
    add_nc("nc-free", cfg_free, "tf", "NC tf, theoretical N (no cap)")
    add_nc("nc-floor4", cfg_fl4, "tf", "NC tf, N >= 4 periods floor below 53 Hz")
    add_nc("nc-plain", cfg0, "plain", "NC plain (no reassignment)")

    # CQT: 标定 + 网格对齐
    cfg = cfg0
    n_frames = cfg.n_frames(x.size)
    pad = cfg.n_cols_sub + 8
    bpo = cqt_bins_per_octave(cfg)
    tone = np.sin(2 * np.pi * 1000.0 * np.arange(int(0.5 * cfg.fs)) / cfg.fs)
    cal = float(cqt_magnitude(tone, cfg, bpo, cfg.hop).max())
    out["cqt"] = (cqt_image(x, cfg, n_frames, pad, bpo, cal), cfg, pad,
                  f"librosa CQT (bpo={bpo}, no reassignment)")
    return out


# ------------------------------------------------------------
def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    cfg = nr.Config()
    x, t_delta = test_signal(cfg)
    print(f"test signal: {TONE_HZ:g} Hz tone + log chirp {CHIRP_F0:g} Hz -> "
          f"{chirp_freq(DURATION) / 1000:.1f} kHz + delta @ {t_delta:.2f} s "
          f"({DURATION:g} s @ {cfg.fs:g} Hz)")

    images = build_images(x)
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 8.2))
    norm = Normalize(vmin=DB_MIN, vmax=DB_MAX)
    t_all = np.linspace(0.05, DURATION - 0.05, 200)
    for ax, (name, (db, cfg_i, pad, note)) in zip(axes.ravel(), images.items()):
        times = cfg_i.col_times(db.shape[1], pad)
        m = ax.pcolormesh(times, nr.log_row_centers(cfg_i), db, shading="auto", cmap=CMAP,
                          norm=norm, rasterized=True)
        ax.set_yscale("log")
        ax.set_ylim(cfg_i.f_min, cfg_i.f_max)
        ax.set_xlim(0.0, DURATION)          # 各方法时间轴一致(切掉两侧垫列)
        ax.axhline(TONE_HZ, color="cyan", ls="--", lw=0.7, alpha=0.85)
        ax.plot(t_all, chirp_freq(t_all), color="cyan", ls="--", lw=0.7, alpha=0.85)
        ax.axvline(t_delta, color="cyan", ls=":", lw=0.7, alpha=0.85)
        ax.set_title(name, fontsize=11)
        ax.set_xlabel("Time (s)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Frequency (Hz)")
    fig.colorbar(m, ax=axes.ravel().tolist(), label="dB (unit tone = 0 dB)",
                 fraction=0.015, pad=0.01)
    fig.suptitle(f"Test: {TONE_HZ:g} Hz tone + log chirp {CHIRP_F0:g} Hz\u2192"
                 f"{chirp_freq(DURATION) / 1000:.1f} kHz + delta @ {t_delta:.2f} s "
                 f"(cyan dashed = reference)", fontsize=12)
    fig.subplots_adjust(hspace=0.28)
    out = os.path.join(OUT_DIR, "compare_all.png")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("saved:", out)


if __name__ == "__main__":
    main()
