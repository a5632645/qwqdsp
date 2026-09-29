# -*- coding: utf-8 -*-
"""
probe_n_limit.py
================

复现并解释: **无窗 NC 的 bin 窗长 N 不能为了时间分辨率随便截短**——截短到某个
程度后整个 NC bin 会失效(输出抖动、增益错误)。

结论(机理)
-----------
一个 NC bin 由**同一个长度 N 的滑动 DFT 的两个相邻 bin** 构成:

    f_left = f_c − Fs/(2N),   f_right = f_c + Fs/(2N),   W_NC = Fs/N

`X_R·conj(X_L)` 的相位差恒为 `−π(N−1)/N ≈ −π`(与信号频率无关), NC 幅度
`max(0, −Re(·))` 正是靠这个"负相关"才成为 bin 的带通检测器。但 N 变短会同时破坏
两件事:

1. **半周期硬边界**: `f_left ≤ 0 ⟺ N ≤ Fs/(2 f_c)`。左分量越过 DC 后不再是
   "更低频的邻箱", 而与 f_c 的**负频镜像**重合(此时二者幅度相等, 比值 = 1)。
   绝对下限就是 **N > Fs/(2 f_c)**, 即 **bin 带宽 Fs/N < 2 f_c**, 也就是**窗内
   至少要装半个周期的中心频率**。
2. **镜像泄漏(渐进)**: `N ≫ Fs/(2f_c)` 时负频镜像的幅度
   `|D(−f_c−f)|/|D(f_c−f)|` 虽小但随 N 减小迅速变大。该镜像与正频主瓣在叉积里
   形成 `e^{±2jφ}` 交叉项 → NC 输出叠加一个 **2·f_c 的拍**, 幅度检测器开始
   "抖"。实测: 镜像/主瓣 ≤3%(≈10 周期)时抖动 <0.1%, 33%(1 周期)时 5%,
   ≥100%(≤半周期)时 50%~80% 且均值明显偏离 A/π。

因此规则是:

    N >  Fs/(2 f_c)   硬性(半周期, 否则 f_left ≤ 0 / 镜像与主瓣重合)
    N ≳  Fs/f_c       可用质量(1 个周期, 抖动约 5%)
    N ≳ 2Fs/f_c       干净(2 个周期, 抖动约 1%)
    自然设计 ≈ 44 周期(每 octave 44 行的 log 网格 + 空窗)总是安全

第二条独立约束(与 f_c 无关, 只与帧步进有关): **N ≥ hop**。窗短于每帧新进样本数时,
逐帧采样的 NC 输出会漏掉瞬态(`verify_log_chirp.py --selftest` 里 310 个 bin 只 233 个有响应),
C++ `WindowlessNcFrame` 为此准备了逐样本 EMA(当前代码里 `alpha = 1.0`, 即**被关闭**)。

用法
----
    python probe_n_limit.py            # 扫描表 + 边界图 + 显示级复现, 图存 output/n_limit.png
    python probe_n_limit.py --quick    # 只打印两张表, 不出图
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nc_reassign as nr          # noqa: E402
import image_metrics as im        # noqa: E402

FS = 48000.0
HOP = 256
CMAP = "magma"
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")


# ------------------------------------------------------------
# 解析量: 单分量 DFT 的幅度 |Σ_{m<N} e^{j2π·Hz·m/Fs}|
# ------------------------------------------------------------
def dft_mag(hz, N, fs: float = FS) -> np.ndarray:
    """长度 N 的 DFT 在频率 hz 处的幅度(矩形窗, 解析式)。"""
    hz = np.asarray(hz, dtype=float)
    num = np.sin(np.pi * hz * N / fs)
    den = np.sin(np.pi * hz / fs)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.abs(num / den)
    out = np.where(np.isclose(den, 0.0), float(N), out)
    return np.where(np.isclose(num, 0.0), 0.0, out)


def pair(fc: float, N: float, fs: float = FS) -> tuple[float, float]:
    return fc - fs / (2 * N), fc + fs / (2 * N)


def contamination(fc: float, N: float) -> float:
    """左右分量里「负频镜像 / 正频主瓣」的较大者(纯音居中在 f_c 时)。"""
    out = []
    for f in pair(fc, N):
        main = dft_mag(fc - f, N)
        imag = dft_mag(-fc - f, N)
        out.append(imag / max(main, 1e-300))
    return float(max(out))


# ------------------------------------------------------------
# 数值仿真: 稳态纯音下该 bin 的 NC 输出健康度
# ------------------------------------------------------------
def simulate(fc: float, N: int, dur: float = 1.0, amp: float = 1.0) -> dict:
    """返回 (增益均值, 抖动 CV, 抖动主频, f_left, 是否有定义)。"""
    fs = FS
    t = np.arange(int(dur * fs)) / fs
    x = amp * np.sin(2 * np.pi * fc * t)
    fl, fr = pair(fc, N)
    XL = nr.sliding_dft(x, fl, N, fs)
    XR = nr.sliding_dft(x, fr, N, fs)
    s = slice(4 * N, x.size - 4 * N if x.size > 8 * N else x.size - 1)
    g = np.sqrt(np.maximum(-(XL.real * XR.real + XL.imag * XR.imag), 0.0))[s] / N
    gg = g - g.mean()
    sp = np.abs(np.fft.rfft(gg * np.hanning(gg.size)))
    fdom = np.fft.rfftfreq(gg.size, 1 / fs)[np.argmax(sp)]
    return {"mean": float(g.mean()), "cv": float(g.std() / max(g.mean(), 1e-30)),
            "f_chatter": float(fdom), "f_left": fl, "expect": amp / np.pi * np.cos(np.pi / N)}


# ------------------------------------------------------------
# 表 1: 固定 f_c, 扫 N
# ------------------------------------------------------------
def table_sweep_n(fc: float = 200.0) -> None:
    n_list = (4800, 2400, 1200, 600, 480, 240, 160, 120, 100, 60, 30, 16, 8)
    print(f"\n== 扫 N (f_c = {fc:g} Hz, 半周期 Fs/2f_c = {FS / (2 * fc):.0f} 样本, "
          f"理论增益 A/π = {1 / np.pi:.4f}) ==")
    print(f"{'N':>6} {'周期Nf/Fs':>10} {'f_left(Hz)':>11} {'镜像/主瓣':>9} {'增益均值':>9} "
          f"{'CV':>7} {'增益/(A/π)':>10} {'抖动主频':>9}")
    for N in n_list:
        r = simulate(fc, N)
        c = contamination(fc, N)
        flag = "  <- 破(过 DC)" if r["f_left"] <= 0 else ""
        print(f"{N:6d} {N * fc / FS:10.2f} {r['f_left']:11.1f} {c:9.3f} {r['mean']:9.4f} "
              f"{r['cv']:7.3f} {r['mean'] / (1 / np.pi):10.3f} {r['f_chatter']:7.1f}Hz{flag}")


# ------------------------------------------------------------
# 表 2: 固定 N, 扫 f_c (等价的失效边界)
# ------------------------------------------------------------
def table_sweep_fc(N: int = 256) -> None:
    print(f"\n== 扫 f_c (N = {N} 固定): 边界 = Fs/(2N) = {FS / (2 * N):.1f} Hz ==")
    print(f"{'f_c(Hz)':>9} {'f_left':>9} {'镜像/主瓣':>9} {'增益均值':>9} {'CV':>7} {'增益/(A/π)':>10}")
    for fc in (25.0, 50.0, 75.0, 93.75, 120.0, 200.0, 500.0, 1000.0, 4000.0):
        r = simulate(fc, N)
        print(f"{fc:9.2f} {r['f_left']:9.1f} {contamination(fc, N):9.3f} {r['mean']:9.4f} "
              f"{r['cv']:7.3f} {r['mean'] / (1 / np.pi):10.3f}")


# ------------------------------------------------------------
# 边界图: (f_c, N) 平面上的镜像污染 / 与仿真的吻合
# ------------------------------------------------------------
def boundary_map() -> plt.Figure:
    fcs = np.geomspace(20.0, 4000.0, 90)
    ns = np.geomspace(16.0, 8192.0, 90)
    FC, NN = np.meshgrid(fcs, ns, indexing="ij")
    cont = np.empty_like(FC)
    for i in range(FC.shape[0]):
        for j in range(FC.shape[1]):
            cont[i, j] = contamination(FC[i, j], NN[i, j])
    half_period = FS / (2 * fcs)                   # f_left = 0 的 N (1-D, 随 f_c)

    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.6))
    ax = axes[0]
    m = ax.pcolormesh(FC, NN, cont, shading="auto", norm=LogNorm(vmin=1e-3, vmax=3.0), cmap="inferno")
    ax.plot(fcs, half_period, color="cyan", lw=1.8, ls="--", label="N = Fs/(2·f_c): f_left = 0 (hard)")
    ax.plot(fcs, FS / fcs, color="w", lw=1.2, ls=":", label="N = Fs/f_c: 1 period (CV ~5%)")
    ax.plot(fcs, 2 * FS / fcs, color="0.7", lw=1.2, ls=":", label="N = 2Fs/f_c: 2 periods (CV ~1%)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("bin center f_c (Hz)"); ax.set_ylabel("window length N (samples)")
    ax.set_title("negative-frequency image / main lobe (left+right comp.)")
    ax.legend(loc="lower left", fontsize=8)
    fig.colorbar(m, ax=ax, label="image / main")

    ax = axes[1]
    pts = [(200.0, N) for N in (8, 16, 30, 60, 100, 120, 160, 240, 480, 1200, 4800)]
    pts += [(20.2, 3600), (30.0, 3600), (100.0, 3600), (500.0, 512), (2000.0, 128)]
    xs = [contamination(fc, n) for fc, n in pts]
    ys = [max(simulate(fc, int(n))["cv"], 1e-5) for fc, n in pts]
    ax.loglog(xs, ys, "o", color="tab:red", ms=6)
    ax.set_ylim(1e-5, 2.0)
    for (fc, n), x, y in zip(pts, xs, ys):
        ax.annotate(f"f_c={fc:g},N={int(n)}", (x, y), fontsize=7, alpha=0.8,
                    xytext=(3, 3), textcoords="offset points")
    ax.axvline(1.0, color="cyan", ls="--", lw=1.2)
    ax.set_xlabel("image / main (analytic)"); ax.set_ylabel("NC output CV (measured)")
    ax.set_title("analytic contamination predicts the measured chatter")
    ax.grid(True, which="both", alpha=0.25)
    return fig


# ------------------------------------------------------------
# 显示级复现: 为了时间分辨率把窗长上限调短 → 低频 bin 整片失效
# ------------------------------------------------------------
def display_demo() -> plt.Figure:
    import verify_log_chirp as V
    x, t, f_true = V.log_chirp()
    row_logf = np.log10(nr.log_row_centers(nr.Config()))
    t_span = (0.25, 1.75)
    cfgs = [("max_window_s = 0.075 s (N<=3600, default)", 0.075),
            ("max_window_s = 0.020 s (N<=960)", 0.020),
            ("max_window_s = 0.008 s (N<=384)", 0.008)]
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 5.0))
    print("\n== 显示级: 缩短窗长上限后, [20,120] Hz 段的线上能量 ==")
    for ax, (label, mws) in zip(axes, cfgs):
        cfg = nr.Config(max_window_s=mws)
        bins = nr.build_nc_bins(cfg)
        n_cols = cfg.n_frames(x.size)
        series = nr.build_series(x, cfg, "nc", bins)
        cal = nr.calibrate(cfg, "nc", "tf", bins=bins)
        table = nr.table_from_series(series, cfg, "tf")
        img_lin = nr.render(table, cfg, n_cols, V.COL_SHIFT)
        img_db = nr.to_db(img_lin, cal, cfg.db_floor)
        times = cfg.col_times(n_cols + V.COL_SHIFT, V.COL_SHIFT)
        m = ax.pcolormesh(times, nr.log_row_centers(cfg), img_db, shading="auto",
                          cmap="magma", vmin=cfg.db_floor, vmax=0.0, rasterized=True)
        ax.set_yscale("log"); ax.set_ylim(cfg.f_min, 400.0)
        ax.set_xlim(times[V.COL_SHIFT], times[-1])
        ax.set_title(label, fontsize=9); ax.set_xlabel("Time (s)")
        ax.plot(np.linspace(*t_span, 200), f_true(np.linspace(*t_span, 200)),
                color="cyan", ls="--", lw=0.9)
        met = im.describe(img_lin, img_db, row_logf, cfg, n_cols, V.COL_SHIFT, f_true,
                          t_span, f_band=(20.0, 120.0))
        print(f"  {label}: [20,120] Hz 线上能量={met['on_line']:.3f} "
              f"脊线偏差med={met['err_cents_med']:.0f} cent; 该上限下 f_left=0 的边界 "
              f"f_c = {1.0 / (2.0 * mws):.1f} Hz")
    axes[0].set_ylabel("Frequency (Hz)")
    fig.suptitle("shortening the NC window cap kills the low-frequency bins "
                 "(NC-tf, log chirp)", fontsize=12)
    return fig




# ------------------------------------------------------------
# N 与 hop: 短窗逐帧采样的漏检
# ------------------------------------------------------------
def hop_demo() -> None:
    """第二条约束: N 与帧步进 hop 的关系(与 f_c 无关的采样问题)。

    逐帧读出的 NC 输出只每隔 hop 个样本采一次; 对**瞬态**而言, 检测器在
    ⟨m⟩ ∈ (N/4, 3N/4) 内才给正值(即窗内位置落在中间一半), 于是单帧栅格上的
    检出概率 ≈ min(1, N/(2·hop))。N < 2·hop 时瞬态会被整帧漏掉(显示闪断)。
    """
    rng = np.random.default_rng(0)
    n = 60000
    frame_idx = (np.arange(1, (n - 4096) // HOP) * HOP + 4095).astype(int)
    positions = rng.integers(2000, 56000, 60)
    print("\n== N vs hop(帧步进 256): 冲激的检出概率(P = 至少一帧采到检测器开窗) ==")
    print(f"{'N':>6} {'N/hop':>7} {'实测 P(检出)':>13} {'理论 min(1,N/2hop)':>19}")
    for N in (3600, 1024, 512, 256, 128, 64, 32, 16, 8):
        hit = []
        for n0 in positions:
            imp = np.zeros(n)
            imp[n0] = 1.0
            XL = nr.sliding_dft(imp, 1000.0 - FS / (2 * N), N, FS)
            XR = nr.sliding_dft(imp, 1000.0 + FS / (2 * N), N, FS)
            g = np.sqrt(np.maximum(-(XL.real * XR.real + XL.imag * XR.imag), 0.0))
            hit.append(bool(np.any(g[frame_idx] > 1e-12)))     # 任何一帧采到即算检出
        print(f"{N:6d} {N / HOP:7.2f} {np.mean(hit):13.3f} {min(1.0, N / (2 * HOP)):19.3f}")


# ------------------------------------------------------------
# ------------------------------------------------------------
# 超低频到底用几个周期的窗(固定秒数上限 vs 周期数上限 vs 理论值)
# ------------------------------------------------------------
def period_demo(dur: float = 2.0, pad_s: float = 2.3,
                f_start: float = 12.5, octaves: float = 10.0) -> None:
    """把测试 chirp 下探到 12.5 Hz, 让 20-40 Hz 那几个 ``N < 2 周期`` 的 bin 真正被激励。

    比较: C++ 的固定秒数上限(0.075 s) / 周期数上限 k=2,4,8 / 理论值(Q=44),
    全部走 NC 时频重分配(tf)。
    """
    import verify_log_chirp as V
    cfg_probe = nr.Config(max_window_s=1e3)          # 逐 bin 探针用的帧栅格
    x0, _, f_true0 = V.log_chirp(duration=dur, f0=f_start, octaves=octaves)
    pad = int(pad_s * FS)
    x = np.concatenate([np.zeros(pad), x0, np.zeros(int(0.2 * FS))])
    f_true = lambda tt: f_true0(np.asarray(tt) - pad_s)          # noqa: E731
    t_span = (pad_s + 0.15, pad_s + dur - 0.25)
    cfgs = [
        ("0.075s cap", nr.Config()),
        ("2 periods", nr.Config(max_window_s=1e3, max_periods=2)),
        ("4 periods", nr.Config(max_window_s=1e3, max_periods=4)),
        ("8 periods", nr.Config(max_window_s=1e3, max_periods=8)),
        ("theoretical", nr.Config(max_window_s=1e3)),
        ("floor 2 periods", nr.Config(min_periods=2)),
        ("floor 4 periods", nr.Config(min_periods=4)),
        ("floor 8 periods", nr.Config(min_periods=8)),
    ]
    bands = [("超低频 20-40Hz", (20.5, 40.0)), ("低频 20-120Hz", (20.5, 120.0)),
             ("全带 20Hz-20kHz", (20.5, 20000.0))]
    print(f"\n== 超低频窗长策略对比 (chirp {f_start:g}Hz -> {f_start * 2 ** octaves / 1000:g}kHz, "
          f"前置 {pad_s:.1f}s 静音, tf 重分配) ==")
    hdr = f"{'策略':<12} {'N@20Hz':>7} {'周期':>5} {'N@1kHz':>7} {'N_max':>7} {'延迟ms':>7} {'标定':>7} |"
    for bn, _ in bands:
        hdr += f" {bn:>16}"
    print(hdr)
    imgs = {}
    for tag, cfg in cfgs:
        bins = nr.build_nc_bins(cfg)
        n_frames = cfg.n_frames(x.size)
        n_max = max(b.N for b in bins)
        pad_cols = cfg.n_cols_sub + 8 + int(np.ceil(n_max / (2 * cfg.hop))) + 8
        pad_right = int(np.ceil(n_max / (2 * cfg.hop))) + 8
        cal = nr.calibrate(cfg, "nc", "tf", duration=0.5, bins=bins)
        series = nr.build_series(x, cfg, "nc", bins)
        table = nr.table_from_series(series, cfg, "tf")
        lin = nr.render(table, cfg, n_frames + pad_right, pad_cols)
        db = nr.to_db(lin, cal, cfg.db_floor)
        row_logf = np.log10(nr.log_row_centers(cfg))
        i1k = int(np.argmin(np.abs(np.array([b.f_center for b in bins]) - 1000.0)))
        line = (f"{tag:<12} {bins[0].N:7d} {bins[0].N * bins[0].f_center / FS:5.2f} "
                f"{bins[i1k].N:7d} {n_max:7d} {n_max / cfg.fs * 1000:7.0f} {cal:7.4f} |")
        for bn, bd in bands:
            m = im.describe(lin, db, row_logf, cfg, n_frames + pad_right, pad_cols,
                            f_true, t_span, f_band=bd)
            line += f" {m['on_line']:7.3f}/{m['err_cents_med']:5.1f}c"
            imgs.setdefault(tag, {})[bn] = (db, cfg, pad_cols, n_frames + pad_right)
        print(line)
    print("(每格 = 线上能量 / 脊线偏差中位 cent)")

    # ---- 逐 bin 数值探针: 低频 bin 在不同周期数下的估计器误差(与图像无关) ----
    print("\n== 逐 bin 数值探针(同一 chirp): N = k 个周期时, 低频 bin 的估计器质量 ==")
    print(f"{'f_c':>6} {'k周期':>6} {'N':>7} {'峰值增益':>9} {'IF误差med':>10} {'IF误差p90':>10} "
          f"{'落点误差med':>11} {'落点误差p90':>11}  (cent / ms; 参考点=窗中心)")
    for fc in (20.0, 30.0, 50.0, 100.0):
        for k in (1.5, 2.0, 4.0, 8.0, 44.0):
            N = int(max(8, round(k * FS / fc)))
            XL = nr.sliding_dft(x, fc - FS / (2 * N), N, FS)
            XR = nr.sliding_dft(x, fc + FS / (2 * N), N, FS)
            idx = np.arange(cfg_probe.n_frames(x.size)) * cfg_probe.hop + cfg_probe.fft_size - 1
            g = np.sqrt(np.maximum(-(XL.real * XR.real + XL.imag * XR.imag), 0.0))[idx] / N
            g[idx < N - 1] = 0.0
            sel = g > 0.5 * g.max()
            if not np.any(sel):
                continue
            i = idx[sel]
            IF = -FS / (2 * np.pi) * np.angle(XL[i] * np.conj(XL[i - 1]))
            t_cen = (i - (N - 1) / 2.0) / FS
            e_if = 1200 * np.log2(np.maximum(IF, 1e-9) / f_true(t_cen))
            tau = N * nr.wrap_pi(np.angle(XR[i] * np.conj(XL[i])) + np.pi) / (2 * np.pi)
            t_dep = (i - (N - 1) / 2.0 + tau) / FS
            t_cross = pad_s + np.log2(fc / f_start) / (octaves / dur)
            e_t = (t_dep - t_cross) * 1000.0
            print(f"{fc:6.0f} {k:6.1f} {N:7d} {g.max():9.4f} "
                  f"{np.median(np.abs(e_if)):10.2f} {np.percentile(np.abs(e_if), 90):10.2f} "
                  f"{np.median(np.abs(e_t)):11.2f} {np.percentile(np.abs(e_t), 90):11.2f}")

    fig, axes = plt.subplots(1, 4, figsize=(19.0, 4.6))
    ax = axes[0]
    for tag, cfg in cfgs:
        bins = nr.build_nc_bins(cfg)
        ax.semilogy([b.f_center for b in bins], [b.N for b in bins], lw=1.3, label=tag)
    ax.set_xscale("log"); ax.set_xlabel("bin center (Hz)"); ax.set_ylabel("N (samples)")
    ax.set_title("window length policies", fontsize=10)
    ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=7)
    for ax, tag in zip(axes[1:], ("0.075s cap", "floor 4 periods", "theoretical")):
        db, cfg, pad_cols, n_buf = imgs[tag]["超低频 20-40Hz"]
        times = cfg.col_times(db.shape[1], pad_cols)
        vmin = max(cfg.db_floor, float(db.max()) - 55.0)
        m = ax.pcolormesh(times, nr.log_row_centers(cfg), db, shading="auto", cmap=CMAP,
                          vmin=vmin, vmax=db.max(), rasterized=True)
        ax.set_yscale("log"); ax.set_ylim(20.0, 300.0)
        ax.set_xlim(t_span[0] - 0.05, t_span[0] + 0.85)
        ax.set_xlabel("Time (s)")
        ax.set_title(f"NC-tf, 20-300 Hz — {tag}", fontsize=10)
        tt = np.linspace(t_span[0] - 0.05, t_span[0] + 0.85, 200)
        ax.plot(tt, f_true(tt), color="cyan", ls="--", lw=0.9)
    axes[1].set_ylabel("Frequency (Hz)")
    fig.suptitle("ultra-low NC window: fixed-time cap vs period cap vs theoretical N", fontsize=12)
    out = os.path.join(OUT_DIR, "period_cap_compare.png")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("saved:", out)


def main() -> None:
    ap = argparse.ArgumentParser(description="NC bin 窗长 N 的下限复现实验")
    ap.add_argument("--quick", action="store_true", help="只打印表, 不出图")
    ap.add_argument("--only", choices=("tables", "boundary", "display", "period"), default=None,
                    help="只跑某一项(便于单独复跑)")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    if args.only in (None, "tables"):
        table_sweep_n(200.0)
        table_sweep_fc(256)
        hop_demo()
    if args.quick:
        return
    if args.only in (None, "boundary"):
        f1 = boundary_map()
        f1.savefig(os.path.join(OUT_DIR, "n_limit_map.png"), dpi=130, bbox_inches="tight")
        plt.close(f1)
        print("saved:", os.path.join(OUT_DIR, "n_limit_map.png"))
    if args.only in (None, "display"):
        f2 = display_demo()
        f2.savefig(os.path.join(OUT_DIR, "n_limit_display.png"), dpi=130, bbox_inches="tight")
        plt.close(f2)
        print("saved:", os.path.join(OUT_DIR, "n_limit_display.png"))
    if args.only in (None, "period"):
        period_demo()


if __name__ == "__main__":
    main()
