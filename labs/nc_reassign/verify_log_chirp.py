# -*- coding: utf-8 -*-
"""
verify_log_chirp.py
===================

用**对数 chirp** 验证无窗 NC 的时间-频率重分配, 并与标准加窗 STFT 重分配对照。

验收判据(用户需求)
------------------
1. 三个变体(freq / time / tf)都要能跑;
2. 对数 chirp 的频谱图应当是**一条锐利的直线**;
3. 无窗 NC 的时间+频率重分配**不能比**标准 STFT 的时间+频率重分配**更差**。

判据 2/3 的量化: 谱图当图像来分析(见 ``image_metrics.py``)
  - 脊线偏差(cent):  逐列能量重心 vs 解析真值
  - 线上能量占比:    真值线 ±50 cent 内的功率占比(越高越锐利)
  - Hough 支撑度:    二值化后落在同一直线上的像素占比(越高越直)
  - 虚假能量占比:    离真值线 ≥ 200 cent 的功率占比(低频镜像/伪影)

用法
----
    python verify_log_chirp.py                 # 对数 chirp, 出图 + 打印指标表
    python verify_log_chirp.py --case crossed  # 上行+下行两条 chirp(多分量)
    python verify_log_chirp.py --selftest      # 只跑约定自检(纯音/冲激), 不出图
"""
from __future__ import annotations

import argparse
import dataclasses
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nc_reassign as nr          # noqa: E402
import image_metrics as im        # noqa: E402

CFG = nr.Config()
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
COL_SHIFT = CFG.n_cols_sub + 8            # 左垫列数(容纳"过去"的落点)
CMAP = "magma"


# ------------------------------------------------------------
# 测试信号
# ------------------------------------------------------------
def log_chirp(duration: float = 2.0, f0: float = 25.0, octaves: float = 9.0,
              fs: float = 48000.0, amp: float = 1.0):
    """对数(指数)chirp: f(t) = f0·2^(k·t), k = octaves/duration。

    返回 (x, t, f_true(t) 的函数)。
    """
    k = octaves / duration
    a = k * np.log(2.0)
    t = np.arange(int(duration * fs)) / fs
    phase = 2 * np.pi * f0 * (np.exp(a * t) - 1.0) / a
    f_true = lambda tt: f0 * 2.0 ** (k * np.asarray(tt, dtype=float))   # noqa: E731
    return amp * np.sin(phase), t, f_true


def crossed_chirps(duration: float = 2.0, f0: float = 25.0, octaves: float = 9.0,
                   fs: float = 48000.0, amp: float = 0.7):
    """上行 + 下行两条对数 chirp(交叉), 用于多分量鲁棒性检查。"""
    x, t, f_up = log_chirp(duration, f0, octaves, fs, amp)
    f1 = f0 * 2.0 ** octaves
    k = octaves / duration
    a = k * np.log(2.0)
    phase = 2 * np.pi * f1 * (1.0 - np.exp(-a * t)) / a
    f_dn = lambda tt: f1 * 2.0 ** (-k * np.asarray(tt, dtype=float))    # noqa: E731
    return x + amp * np.sin(phase), t, f_up, f_dn


# ------------------------------------------------------------
# 每种方法/变体的图像 + 指标
# ------------------------------------------------------------
def images_for(x: np.ndarray, kind: str, variants=nr.VARIANTS):
    """算一个方法的所有变体: 返回 {variant: (img_lin, img_db, cal)}。"""
    series = nr.build_series(x, CFG, kind)
    n_cols = CFG.n_frames(x.size)
    out = {}
    for v in variants:
        cal = nr.calibrate(CFG, kind, v)
        table = nr.table_from_series(series, CFG, v)
        img_lin = nr.render(table, CFG, n_cols, COL_SHIFT)
        out[v] = (img_lin, nr.to_db(img_lin, cal, CFG.db_floor), cal)
    return out


def metrics_for(img_lin, img_db, row_logf, n_cols, f_true, t_span, tol_cents=50.0,
                f_band=(20.0, 20000.0)):
    return im.describe(img_lin, img_db, row_logf, CFG, n_cols, COL_SHIFT,
                       f_true, t_span, tol_cents=tol_cents, f_band=f_band)


# ------------------------------------------------------------
# 作图
# ------------------------------------------------------------
def plot_panels(panels, out_path, title, f_trues, t_span, n_cols, col_shift=COL_SHIFT):
    """panels: [(标题, img_db), ...]; 每格叠画真值线(``f_trues`` 为曲线列表, 可为 None)。"""
    freqs = nr.log_row_centers(CFG)
    times = CFG.col_times(n_cols + col_shift, col_shift)
    norm = Normalize(vmin=CFG.db_floor, vmax=0.0)
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(4.6 * n, 5.2), squeeze=False)
    tt = np.linspace(t_span[0], t_span[1], 400)
    for ax, (label, img_db) in zip(axes[0], panels):
        m = ax.pcolormesh(times, freqs, img_db, shading="auto", cmap=CMAP, norm=norm,
                          rasterized=True)
        ax.set_yscale("log")
        ax.set_ylim(CFG.f_min, CFG.f_max)
        ax.set_xlim(times[col_shift], times[-1])
        ax.set_title(label, fontsize=10)
        ax.set_xlabel("Time (s)")
        for fn in f_trues or ():
            ax.plot(tt, fn(tt), color="cyan", ls="--", lw=0.9, alpha=0.9)
    axes[0][0].set_ylabel("Frequency (Hz)")
    fig.colorbar(m, ax=list(axes[0]), label="dB", fraction=0.02)
    fig.suptitle(title, fontsize=13)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("saved:", out_path)


def plot_compare(img_by_kind, out_path, title, f_trues, t_span, n_cols):
    """2×2: 两种方法 × {plain, tf}。"""
    freqs = nr.log_row_centers(CFG)
    times = CFG.col_times(n_cols + COL_SHIFT, COL_SHIFT)
    norm = Normalize(vmin=CFG.db_floor, vmax=0.0)
    fig, axes = plt.subplots(2, 2, figsize=(13, 7.4))
    tt = np.linspace(t_span[0], t_span[1], 400)
    for r, kind in enumerate(("nc", "stft")):
        for c, v in enumerate(("plain", "tf")):
            ax = axes[r][c]
            img_db = img_by_kind[kind][v][1]
            m = ax.pcolormesh(times, freqs, img_db, shading="auto", cmap=CMAP, norm=norm,
                              rasterized=True)
            ax.set_yscale("log")
            ax.set_ylim(CFG.f_min, CFG.f_max)
            ax.set_xlim(times[COL_SHIFT], times[-1])
            ax.set_title(f"{'NC' if kind == 'nc' else 'STFT'} · {v}", fontsize=10)
            for fn in f_trues or ():
                ax.plot(tt, fn(tt), color="cyan", ls="--", lw=0.9, alpha=0.9)
            if r == 1:
                ax.set_xlabel("Time (s)")
            if c == 0:
                ax.set_ylabel("Frequency (Hz)")
    fig.colorbar(m, ax=axes.ravel().tolist(), label="dB", fraction=0.02)
    fig.suptitle(title, fontsize=13)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("saved:", out_path)


def plot_metrics(metrics, out_path, title):
    """指标柱状图: 脊线偏差 / 线上能量 / 虚假能量 / Hough 支撑度。"""
    keys = [("err_cents_med", "ridge |err| median (cent)", False),
            ("width_oct_med", "10-90 width (octave)", False),
            ("on_line", "energy on line (±50 cent)", True),
            ("spurious", "spurious energy (≥200 cent)", False),
            ("support", "Hough support", True)]
    names = list(metrics.keys())
    fig, axes = plt.subplots(1, len(keys), figsize=(4.0 * len(keys), 3.9))
    for ax, (key, label, higher_better) in zip(axes, keys):
        vals = [metrics[n].get(key, np.nan) for n in names]
        ax.bar(range(len(names)), vals, color=["tab:blue"] * 4 + ["tab:orange"] * 4)
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
        ax.set_title(label, fontsize=9)
    fig.suptitle(title, fontsize=12)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("saved:", out_path)


# ------------------------------------------------------------
# 自检: 时间参考与算子符号(用纯音与冲激)
# ------------------------------------------------------------
def selftest() -> None:
    """验证全库的时间参考约定与算子符号。

    1) 静止纯音(1 kHz): NC 的 dt 应 ≈ 0(能量落在窗口中心), 瞬时频率误差 ≈ 0 cent;
    2) 冲激: NC/STFT 的落点时刻应 ≈ 冲激真实时刻(不受帧边界影响)。
    """
    cfg = CFG
    print("== 自检 1: 静止纯音 (f0 在 log 行中心) ==")
    centers = nr.log_row_centers(cfg)
    f0 = float(centers[np.argmin(np.abs(np.log(centers) - np.log(1000.0)))])
    n = int(1.0 * cfg.fs) + cfg.fft_size
    t = np.arange(n) / cfg.fs
    tone = np.sin(2 * np.pi * f0 * t)
    bins = nr.build_nc_bins(cfg)
    ser = nr.nc_series(tone, cfg, bins)
    centers_bins = np.array([b.f_center for b in bins])
    i0 = int(np.argmin(np.abs(centers_bins - f0)))
    b = bins[i0]
    g = ser.gain[i0]
    sel = g > 0.5 * g.max()
    print(f"  NC bin f_c={b.f_center:8.1f} Hz N={b.N:5d}: "
          f"IF={np.median(ser.if_hz[i0][sel]):9.2f} Hz "
          f"(误差 {1200 * np.log2(np.median(ser.if_hz[i0][sel]) / f0):+.2f} cent, 应为 0); "
          f"dt={np.median(ser.dt[i0][sel]):+.2f} 样本(应为 0), "
          f"|dt| max={np.max(np.abs(ser.dt[i0][sel])):.1f} (N/2={b.N / 2:.0f})")
    sst = nr.stft_series(tone, cfg)
    k0 = int(round(f0 / (cfg.fs / cfg.fft_size)))
    gs = sst.mag[k0]
    sel2 = gs > 0.5 * gs.max()
    print(f"  STFT bin k={k0}: IF={np.median(sst.if_hz[k0][sel2]):9.2f} Hz "
          f"(误差 {1200 * np.log2(np.median(sst.if_hz[k0][sel2]) / f0):+.2f} cent, 应为 0); "
          f"gd={np.median(sst.gd[k0][sel2]):+.3f}(应为 ~0.0, 窗中心)")

    print("== 自检 2: 单冲激的落点时刻 ==")
    n0 = 20000
    imp = np.zeros(200000)
    imp[n0] = 1.0
    n_frames = cfg.n_frames(imp.size)
    j = np.arange(n_frames)
    ref = j * cfg.hop + (cfg.fft_size - 1) / 2.0        # 帧中心(绝对样本)
    ser = nr.nc_series(imp, cfg, bins)
    errs, skipped = [], 0
    for i, b in enumerate(bins):
        g = ser.gain[i]
        if g.max() < 1e-6:          # 该 bin 对冲激无响应(窗内 ⟨m⟩/N 不在 1/4..3/4)
            skipped += 1
            continue
        sel = g > 0.2 * g.max()
        dep = ref[sel] + ser.ref_off[i] + ser.dt[i][sel]
        errs.append(np.median(dep - n0))
    errs = np.array(errs)
    print(f"  NC   ({len(errs)} 个有响应 bin, {skipped} 个无响应): 落点误差 中位="
          f"{np.median(np.abs(errs)):.1f} 样本, 90 分位={np.percentile(np.abs(errs), 90):.1f} 样本 "
          f"(hop={cfg.hop} 样本)")
    sst = nr.stft_series(imp, cfg)
    errs2 = []
    for k in range(1, len(sst.labels)):
        g = sst.mag[k]
        if g.max() < 1e-6:
            continue
        sel = g > 0.2 * g.max()
        dep = ref[sel] + sst.gd[k][sel] * cfg.fft_size
        errs2.append(np.median(dep - n0))
    errs2 = np.array(errs2)
    print(f"  STFT ({len(errs2)} 个有响应 bin): 落点误差 中位={np.median(np.abs(errs2)):.1f} 样本, "
          f"90 分位={np.percentile(np.abs(errs2), 90):.1f} 样本")


# ------------------------------------------------------------
# 主流程
# ------------------------------------------------------------
def run_chirp() -> tuple[dict, dict]:
    x, t, f_true = log_chirp()
    t_span = (0.25, 1.75)          # 只统计真值落在显示范围内的列
    f_band = (100.0, 10000.0)      # 低频段(≤100Hz)NC 窗长被夹住, 单独看更公平
    n_cols = CFG.n_frames(x.size)
    row_logf = np.log10(nr.log_row_centers(CFG))
    print(f"信号: 对数 chirp {f_true(0):.1f} -> {f_true(t[-1]):.1f} Hz / {t[-1]:.2f} s, "
          f"{n_cols} 帧, hop={CFG.hop}, fft={CFG.fft_size}")

    img_by_kind = {}
    metrics, metrics_band = {}, {}
    for kind in ("nc", "stft"):
        img_by_kind[kind] = images_for(x, kind)
        for v in nr.VARIANTS:
            lin, db, cal = img_by_kind[kind][v]
            metrics[f"{kind}-{v}"] = metrics_for(lin, db, row_logf, n_cols, f_true, t_span)
            metrics_band[f"{kind}-{v}"] = metrics_for(lin, db, row_logf, n_cols, f_true,
                                                      t_span, f_band=f_band)

    plot_panels([(f"NC · {v}", img_by_kind["nc"][v][1]) for v in ("freq", "time", "tf")],
                os.path.join(OUT_DIR, "chirp_nc.png"),
                "Windowless NC: freq-only / time-only / time+freq reassignment",
                [f_true], t_span, n_cols)
    plot_panels([(f"STFT · {v}", img_by_kind["stft"][v][1]) for v in ("freq", "time", "tf")],
                os.path.join(OUT_DIR, "chirp_stft.png"),
                "Standard STFT: freq-only / time-only / time+freq reassignment",
                [f_true], t_span, n_cols)
    plot_compare(img_by_kind, os.path.join(OUT_DIR, "chirp_compare.png"),
                 "Log chirp: plain vs time+freq reassignment (NC top, STFT bottom)",
                 [f_true], t_span, n_cols)
    plot_metrics(metrics, os.path.join(OUT_DIR, "chirp_metrics.png"), "Log chirp metrics")
    return metrics, metrics_band


def print_table(metrics: dict, title: str) -> None:
    cols = [("err_cents_med", "脊线偏差med(cent)"), ("err_cents_p90", "脊线偏差p90"),
            ("width_oct_med", "10-90宽度(oct)"), ("width_perp_oct_med", "垂直宽度(oct)"),
            ("on_line", "线上能量(±50c)"), ("spurious", "虚假能量(≥200c)"),
            ("support", "Hough支撑"), ("perp_rms_rows", "垂直RMS(行)"),
            ("peaks_med", "每列峰数"), ("pixels", "阈值上像素数"),
            ("slope", "Hough斜率"), ("slope_theory", "斜率真值")]
    w = max(len(k) for k in metrics) + 2
    print(f"\n== 指标表: {title} ==")
    print(" " * w + "".join(f"{c[1]:>19}" for c in cols))
    for name, m in metrics.items():
        print(f"{name:<{w}}" + "".join(
            f"{m.get(k, float('nan')):>19.4g}" if isinstance(m.get(k), float) else
            f"{m.get(k, float('nan')):>19}" for k, _ in cols))


def run_crossed() -> None:
    x, t, f_up, f_dn = crossed_chirps()
    n_cols = CFG.n_frames(x.size)
    t_span = (0.25, 1.75)
    times = CFG.col_times(n_cols + COL_SHIFT, COL_SHIFT)
    cols = np.arange(n_cols + COL_SHIFT)
    sel = cols >= COL_SHIFT
    img_by_kind = {}
    for kind in ("nc", "stft"):
        img_by_kind[kind] = images_for(x, kind)
        for v in ("plain", "tf"):
            p = im.peak_count(img_by_kind[kind][v][1], cols[sel], thr_db=-25.0)
            print(f"  {kind}-{v}: 每列峰数 中位={p['peaks_med']:.0f} 均值={p['peaks_mean']:.2f} "
                  f"(期望 ≈ 2 条 chirp)")
    plot_compare(img_by_kind, os.path.join(OUT_DIR, "crossed_compare.png"),
                 "Crossed log chirps (up+down): plain vs time+freq",
                 [f_up, f_dn], t_span, n_cols)
    plot_panels([(f"NC · {v}", img_by_kind['nc'][v][1]) for v in ("freq", "time", "tf")],
                os.path.join(OUT_DIR, "crossed_nc.png"),
                "Windowless NC on crossed log chirps", [f_up, f_dn], t_span, n_cols)


def estimate_case_variants():
    """频率估计器对比用的三档配置(其余参数同默认)。"""
    return [("weighted", "old: weighted component IF (no extrapolation)"),
            ("taper_center", "sine-taper IF, center reference only"),
            ("taper", "sine-taper IF + 1st-order extrapolation to deposit time (default)")]


def image_for_cfg(cfg, x: np.ndarray, variant: str = "tf"):
    """按给定配置渲染一张 NC 图; 返回 (lin, db, pad, n_cols)。"""
    bins = nr.build_nc_bins(cfg)
    series = nr.build_series(x, cfg, "nc", bins)
    cal = nr.calibrate(cfg, "nc", variant, duration=0.5, bins=bins)
    table = nr.table_from_series(series, cfg, variant)
    n_frames = cfg.n_frames(x.size)
    n_max = max(b.N for b in bins)
    pad = cfg.n_cols_sub + 8 + int(np.ceil(n_max / (2 * cfg.hop))) + 8
    pad_r = int(np.ceil(n_max / (2 * cfg.hop))) + 8
    lin = nr.render(table, cfg, n_frames + pad_r, pad)
    return lin, nr.to_db(lin, cal, cfg.db_floor), pad, n_frames + pad_r


def run_if_estimator_case() -> None:
    """频率估计器改进的前后对比(默认 taper vs 旧 weighted), 出图 + 打印指标。

    上图整带、下图放大低频段(20–300 Hz) —— 差异主要在低频与旁瓣污染处。
    """
    x, _, f_true = log_chirp()
    t_span = (0.25, 1.75)
    row_logf = np.log10(nr.log_row_centers(CFG))
    items = []
    for est, label in estimate_case_variants():
        cfg = dataclasses.replace(CFG, if_estimator=est)
        lin, db, pad, nbuf = image_for_cfg(cfg, x)
        m_all = im.describe(lin, db, row_logf, cfg, nbuf, pad, f_true, t_span)
        m_low = im.describe(lin, db, row_logf, cfg, nbuf, pad, f_true, t_span, f_band=(30, 300))
        items.append((est, label, cfg, db, pad, m_all, m_low))
        print(f"  {est:<13} 全带 线上能量={m_all['on_line']:.4f} 脊线偏差={m_all['err_cents_med']:5.2f}c "
              f"支撑={m_all['support']:.3f} 宽度={m_all['width_oct_med']:.5f}oct | "
              f"30-300Hz 线上能量={m_low['on_line']:.4f} 偏差={m_low['err_cents_med']:5.2f}c   {label}")

    fig, axes = plt.subplots(2, 3, figsize=(16.5, 8.4))
    tt = np.linspace(t_span[0], t_span[1], 400)
    for col, (est, label, cfg, db, pad, m_all, m_low) in enumerate(items):
        times = cfg.col_times(db.shape[1], pad)
        freqs = nr.log_row_centers(cfg)
        for row_i, (ylim, span, met, tag) in enumerate(
                [((cfg.f_min, cfg.f_max), t_span, m_all, "full band"),
                 ((cfg.f_min, 300.0), (0.0, 1.2), m_low, "20-300 Hz zoom")]):
            ax = axes[row_i][col]
            m = ax.pcolormesh(times, freqs, db, shading="auto", cmap=CMAP,
                              vmin=cfg.db_floor, vmax=0.0, rasterized=True)
            ax.set_yscale("log"); ax.set_ylim(*ylim)
            ax.set_xlim(span[0] if span[0] > 0 else times[pad], span[1])
            ax.set_title(f"{est} — {tag}: on-line {met['on_line']:.3f}, "
                         f"ridge {met['err_cents_med']:.1f}c", fontsize=9)
            ax.plot(tt, f_true(tt), color="cyan", ls="--", lw=0.8)
            if col == 0:
                ax.set_ylabel("Frequency (Hz)")
            if row_i == 1:
                ax.set_xlabel("Time (s)")
    fig.colorbar(m, ax=axes.ravel().tolist(), label="dB", fraction=0.015, pad=0.01)
    fig.suptitle("Windowless NC freq-reassignment estimator: old vs sine-taper "
                 "(+ extrapolation to deposit time)", fontsize=12)
    fig.subplots_adjust(hspace=0.3)
    out = os.path.join(OUT_DIR, "if_estimator_compare.png")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("saved:", out)


def window_case_configs() -> list[tuple[str, nr.Config, str]]:
    """低频 NC bin 的两种窗长策略(都做 NC 时频重分配)。"""
    return [
        ("clamped", nr.Config(), "max_window_s=0.075s: N<=3600smp (below ~590 Hz truncated)"),
        ("free", nr.Config(max_window_s=1e3), "N = theoretical eq.(7), no cap (25 Hz -> ~85k smp)"),
    ]


def run_window_case(dur: float = 2.0, pad_s: float = 2.3) -> None:
    """低频 NC bin **不截短**(用理论 DFT 窗长 N) vs 截短, 都做 NC 时频重分配。

    理论窗长在低频极大(25 Hz ≈ 85k 样本 = 1.77 s), 所以测试信号前面垫
    ``pad_s`` 秒静音, 让这些长窗在 chirp 到来前就填满(等价于一个连续运行的
    实时系统); 静音垫段不参与指标统计。
    """
    x0, _, f_true0 = log_chirp(duration=dur)
    pad = int(pad_s * CFG.fs)
    x = np.concatenate([np.zeros(pad), x0, np.zeros(int(0.2 * CFG.fs))])
    f_true = lambda tt: f_true0(np.asarray(tt) - pad_s)          # noqa: E731
    t_span = (pad_s + 0.25, pad_s + dur - 0.25)
    row_logf = np.log10(nr.log_row_centers(CFG))
    low_band = (20.0, 120.0)
    print(f"\n== 低频窗长对比: 对数 chirp {dur:.0f}s (前置 {pad_s:.1f}s 静音), "
          f"分析窗 {t_span[0]:.2f}-{t_span[1]:.2f}s ==")

    imgs, metrics, metrics_low = {}, {}, {}
    for tag, cfg, note in window_case_configs():
        bins = nr.build_nc_bins(cfg)
        n_frames = cfg.n_frames(x.size)
        n_max = max(b.N for b in bins)
        pad_cols = cfg.n_cols_sub + 8 + int(np.ceil(n_max / (2 * cfg.hop))) + 8
        pad_right = int(np.ceil(n_max / (2 * cfg.hop))) + 8
        cal_tf = nr.calibrate(cfg, "nc", "tf", duration=0.5, bins=bins)
        cal_plain = nr.calibrate(cfg, "nc", "plain", duration=0.5, bins=bins)
        print(f"[{tag}] {note}")
        print(f"      N: 最低行 f={bins[0].f_center:.1f}Hz -> N={bins[0].N} "
              f"({bins[0].N / cfg.fs * 1000:.0f} ms = {bins[0].N * bins[0].f_center / cfg.fs:.0f} 周期), "
              f"最高行 N={bins[-1].N}; 最长窗延迟={n_max / cfg.fs * 1000:.0f} ms; "
              f"标定 tf/plain={cal_tf:.4f}/{cal_plain:.4f}")
        imgs[tag] = {}
        series = nr.build_series(x, cfg, "nc", bins)      # 观测量只算一次, 4 个变体共用
        for v in nr.VARIANTS:
            table = nr.table_from_series(series, cfg, v)
            lin = nr.render(table, cfg, n_frames + pad_right, pad_cols)
            db = nr.to_db(lin, cal_tf, cfg.db_floor)
            imgs[tag][v] = (lin, db, cfg, pad_cols, n_frames + pad_right)
            metrics[f"{tag}-{v}"] = im.describe(lin, db, row_logf, cfg, n_frames + pad_right,
                                                pad_cols, f_true, t_span)
            metrics_low[f"{tag}-{v}"] = im.describe(lin, db, row_logf, cfg, n_frames + pad_right,
                                                    pad_cols, f_true, t_span, f_band=low_band)
    print_table(metrics, "全带 20 Hz – 20 kHz (低频窗长对比)")
    print_table(metrics_low, "低频段 20 – 120 Hz (低频窗长对比)")

    # 图: N(f) 曲线 + 两张低频放大图
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.0))
    ax = axes[0]
    for tag, cfg, _ in window_case_configs():
        bins = nr.build_nc_bins(cfg)
        ax.semilogy([b.f_center for b in bins], [b.N for b in bins], lw=1.4, label=tag)
    ax.axhline(CFG.fs * 0.075, color="0.6", ls=":", lw=1.0)
    ax.text(30, CFG.fs * 0.075 * 1.15, "C++ cap 0.075 s = 3600 smp", fontsize=7, color="0.4")
    ax.set_xscale("log"); ax.set_xlabel("bin center (Hz)"); ax.set_ylabel("window length N (samples)")
    ax.set_title("NC window length: theoretical vs clamped", fontsize=10)
    ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=8)
    for ax, tag in zip(axes[1:], ("clamped", "free")):
        lin, db, cfg, pad_cols, n_buf = imgs[tag]["tf"]
        times = cfg.col_times(db.shape[1], pad_cols)
        vmin = max(cfg.db_floor, float(db.max()) - 55.0)     # 各自峰值相对显示(两法低频增益差很大)
        m = ax.pcolormesh(times, nr.log_row_centers(cfg), db, shading="auto",
                          cmap=CMAP, vmin=vmin, vmax=db.max(), rasterized=True)
        ax.set_yscale("log"); ax.set_ylim(cfg.f_min, 400.0)
        ax.set_xlim(t_span[0] - 0.15, t_span[1] + 0.15)
        ax.set_xlabel("Time (s)")
        ax.set_title(f"NC-tf, low band zoom - {tag} (dB rel. own peak)", fontsize=10)
        tt = np.linspace(*t_span, 200)
        ax.plot(tt, f_true(tt), color="cyan", ls="--", lw=0.9)
    axes[1].set_ylabel("Frequency (Hz)")
    fig.colorbar(m, ax=list(axes[1:]), label="dB (rel. own peak)", fraction=0.02)
    fig.suptitle("Low-frequency NC bins: clamped window vs theoretical (uncapped) window",
                 fontsize=12)
    out = os.path.join(OUT_DIR, "window_lowfreq_compare.png")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("saved:", out)


def main() -> None:
    ap = argparse.ArgumentParser(description="无窗 NC 时间-频率重分配的 log chirp 验证")
    ap.add_argument("--case", choices=("chirp", "crossed", "window", "ifest", "all"),
                    default="chirp")
    ap.add_argument("--selftest", action="store_true", help="只跑时间参考自检")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    if args.selftest:
        selftest()
        return
    if args.case in ("chirp", "all"):
        metrics, metrics_band = run_chirp()
        print_table(metrics, "全带 20 Hz – 20 kHz")
        print_table(metrics_band, "限带 100 Hz – 10 kHz")
        print("\n== 验收: 无窗 NC 时间+频率重分配 vs 标准 STFT 时间+频率重分配 ==")
        for band, m in (("全带", metrics), ("限带", metrics_band)):
            nctf, sttf = m["nc-tf"], m["stft-tf"]
            print(f"  [{band}] 线上能量(越高越锐): NC {nctf['on_line']:.4f} vs STFT {sttf['on_line']:.4f}; "
                  f"脊线偏差med(越低越准): NC {nctf['err_cents_med']:.2f} vs STFT {sttf['err_cents_med']:.2f} cent; "
                  f"虚假能量: NC {nctf['spurious']:.2e} vs STFT {sttf['spurious']:.2e}; "
                  f"Hough支撑: NC {nctf['support']:.3f} vs STFT {sttf['support']:.3f}")
    if args.case in ("crossed", "all"):
        run_crossed()
    if args.case in ("window", "all"):
        run_window_case()
    if args.case in ("ifest", "all"):
        print("\n== 频率估计器对比 (--case ifest) ==")
        run_if_estimator_case()


if __name__ == "__main__":
    main()
