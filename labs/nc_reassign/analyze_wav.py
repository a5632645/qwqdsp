# -*- coding: utf-8 -*-
"""
analyze_wav.py
==============

对**真实 wav 文件**跑三种方法并输出频谱图(**只出图, 不做任何数据分析**):

  1. ``stft-tf``    标准加窗 STFT 的时间+频率重分配(C++ ``tf_reassignment_frame`` 等价)
  2. ``nc-floor4``  无窗 NC 时间+频率重分配(超低频 ``N ≥ 4 周期`` 下限)
  3. ``nc-plain``   无窗 NC 原样显示 —— 不做任何重分配

三块面板用同一显示网格, 每种方法各自用**单位幅度稳态正弦**标定到 0 dB, 因此面板间 dB 直接
可比。时间轴按**文件内的绝对时刻**标注。

用法
----
    # 默认: wormhole.wav 第 25013 号样本起 14217 个采样
    python analyze_wav.py

    # 其它文件/片段
    python analyze_wav.py path/to/x.wav --start 0 --length 48000 --out output/x.png
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nc_reassign as nr          # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "output")
INPUT_DIR = os.path.join(HERE, "..", "..", "work_dir", "input")

# 默认片段(用户指定)
DEFAULT_WAV = os.path.join(INPUT_DIR, "wormhole.wav")
DEFAULT_START = 25013
DEFAULT_LENGTH = 14217

CMAP = "magma"


def methods(fs: float):
    """(名称, kind, variant, 配置, 说明) 三个面板。"""
    return (
        ("stft-tf", "stft", "tf", nr.Config(fs=fs), "STFT tf reassign"),
        ("nc-floor4", "nc", "tf", nr.Config(fs=fs, min_periods=4), "NC tf, N >= 4 periods floor"),
        ("nc-plain", "nc", "plain", nr.Config(fs=fs), "NC plain (no reassignment)"),
    )


def read_segment(path: str, start: int, length: int) -> tuple[np.ndarray, float]:
    """只读文件里 [start, start+length) 这段(单声道化)。"""
    import soundfile as sf
    data, fs = sf.read(path, start=start, frames=length, dtype="float64")
    if data.ndim > 1:
        data = data.mean(axis=1)
    return np.asarray(data, dtype=np.float64), float(fs)


def render_method(x: np.ndarray, cfg: nr.Config, kind: str, variant: str):
    """渲染一张 dB 图; 返回 (img_db, col_shift, cfg)。"""
    if kind == "stft":
        table = nr.table_from_series(nr.build_series(x, cfg, "stft"), cfg, variant)
        n_frames = cfg.n_frames(x.size)
        pad = cfg.n_cols_sub + 8
        pad_r = 0
        cal = nr.calibrate(cfg, "stft", variant)
    else:
        bins = nr.build_nc_bins(cfg)
        n_max = max(b.N for b in bins)
        table = nr.table_from_series(nr.build_series(x, cfg, "nc", bins), cfg, variant)
        n_frames = cfg.n_frames(x.size)
        pad = cfg.n_cols_sub + 8 + int(np.ceil(n_max / (2 * cfg.hop))) + 8
        pad_r = int(np.ceil(n_max / (2 * cfg.hop))) + 8
        cal = nr.calibrate(cfg, "nc", variant, duration=0.5, bins=bins)
    img = nr.render(table, cfg, n_frames + pad_r, pad)
    return nr.to_db(img, cal, cfg.db_floor), pad, cfg


def main() -> None:
    ap = argparse.ArgumentParser(description="真实 wav 的三方法频谱图(只出图)")
    ap.add_argument("wav", nargs="?", default=DEFAULT_WAV)
    ap.add_argument("--start", type=int, default=DEFAULT_START, help="起始样本")
    ap.add_argument("--length", type=int, default=DEFAULT_LENGTH, help="采样点数")
    ap.add_argument("--out", default=None, help="默认 output/wav_<文件名>_<起始样本>.png")
    args = ap.parse_args()
    if args.out is None:                      # 按文件/偏移自动命名, 便于连续跑多个片段
        stem = os.path.splitext(os.path.basename(args.wav))[0]
        args.out = os.path.join(OUT_DIR, f"wav_{stem}_{args.start}.png")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    x, fs = read_segment(args.wav, args.start, args.length)
    t0 = args.start / fs
    t1 = (args.start + x.size) / fs
    print(f"{os.path.basename(args.wav)}: 样本 [{args.start}, {args.start + x.size}) "
          f"= {x.size} 点 / {x.size / fs * 1000:.1f} ms @ {fs:g} Hz, t = [{t0:.3f}, {t1:.3f}] s")

    items = []
    for name, kind, variant, cfg, note in methods(fs):
        db, pad, cfg_i = render_method(x, cfg, kind, variant)
        items.append((name, note, db, pad, cfg_i))
        print(f"  {name:<10} {note:<34} 帧数={cfg_i.n_frames(x.size)}")

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8))
    norm = Normalize(vmin=cfg_i.db_floor, vmax=0.0)
    # 统一时间范围 = 各方法都覆盖到的"帧中心跨度"(否则垫列会使面板长度不一致)
    n_frames = items[0][4].n_frames(x.size)
    span_lo = t0 + (cfg_i.fft_size - 1) / 2.0 / fs
    span_hi = span_lo + (n_frames - 1) * cfg_i.hop / fs
    for ax, (name, note, db, pad, cfg_i) in zip(axes, items):
        times = cfg_i.col_times(db.shape[1], pad) + t0
        ax.pcolormesh(times, nr.log_row_centers(cfg_i), db, shading="auto", cmap=CMAP,
                      norm=norm, rasterized=True)
        ax.set_yscale("log")
        ax.set_ylim(cfg_i.f_min, cfg_i.f_max)
        ax.set_xlim(span_lo, span_hi)
        ax.set_title(f"{name} — {note}", fontsize=9.5)
        ax.set_xlabel("Time (s, in file)")
        ax.set_aspect("auto")
    axes[0].set_ylabel("Frequency (Hz)")
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), ax=list(axes),
                 label="dB (unit tone = 0 dB)", fraction=0.015, pad=0.01)
    fig.suptitle(f"{os.path.basename(args.wav)}  samples [{args.start}, "
                 f"{args.start + x.size})  ({x.size / fs * 1000:.0f} ms @ {fs / 1000:g} kHz)",
                 fontsize=12)
    fig.subplots_adjust(wspace=0.18)
    fig.savefig(args.out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("saved:", args.out)


if __name__ == "__main__":
    main()
