# -*- coding: utf-8 -*-
"""
image_metrics.py
================

把重分配谱图当**数字图像**做量化分析, 用于判定「对数 chirp 是否被重分配成
一条锐利的直线」, 以及比较无窗 NC 与标准 STFT 重分配的优劣。

全部指标只依赖 numpy(不引入 scikit-image):

1. ``ridge_stats``   逐列能量重心(脊线)与真值的偏差(cent) + 列内能量展宽
2. ``line_energy``   线上能量占比(±容差内 / 全部) 与虚假分量能量占比
3. ``hough_line``    霍夫变换找最强直线(斜率/支撑度/垂直 RMS), 检验"直线性"
4. ``peak_count``    每列超过阈值的局部极大值个数(虚假谱线的粗略计数)

坐标约定: 图像 ``img[y, c]``, y = log 频率行(0 = 最低频), c = 时间列。
"垂直距离"指垂直于 chirp 轨迹的距离: 轨迹在 (列, 行) 平面的斜率已知时,
行方向的展宽 w_rows 对应垂直距离 w_rows·cosθ。
"""
from __future__ import annotations

import numpy as np

LOG2_10 = np.log2(10.0)


# ------------------------------------------------------------
# 逐列脊线
# ------------------------------------------------------------
def ridge_stats(img_lin: np.ndarray, row_logf: np.ndarray, gt_logf: np.ndarray,
                cols: np.ndarray) -> dict:
    """逐列能量重心 vs 真值。

    参数
    ----
    img_lin  : (n_rows, n_cols) 线性幅度图
    row_logf : (n_rows,) 每行的 log10(中心频率)
    gt_logf  : (n_cols,) 每列的真值 log10(瞬时频率)
    cols     : 参与统计的列索引(真值落在显示范围内、且不是被垫的空列)

    返回 ``dict``:
      err_cents_med / err_cents_p90 : 脊线偏差(cent)的中位/90 分位绝对值
      width_oct_med / width_oct_p90 : 列内能量 10-90 展宽(octave, 未除法向)
      cols_used
    """
    p = np.square(img_lin[:, cols].astype(np.float64))
    tot = p.sum(axis=0)
    ok = tot > 0
    if not np.any(ok):
        return {"err_cents_med": np.nan, "err_cents_p90": np.nan,
                "width_oct_med": np.nan, "width_oct_p90": np.nan, "cols_used": 0}
    q = p[:, ok] / tot[ok]
    ridge = q.T @ row_logf
    err = (ridge - gt_logf[cols][ok]) * LOG2_10 * 1200.0            # cent
    # 10-90 分位展宽(对 log10 频率做线性插值, 避免行量化)
    cdf = np.cumsum(q, axis=0)
    width = np.empty(q.shape[1])
    for i in range(q.shape[1]):
        c = cdf[:, i]
        width[i] = (np.interp(0.90, c, row_logf) - np.interp(0.10, c, row_logf)) / np.log10(2.0)
    return {
        "err_cents_med": float(np.median(np.abs(err))),
        "err_cents_p90": float(np.percentile(np.abs(err), 90)),
        "width_oct_med": float(np.median(width)),
        "width_oct_p90": float(np.percentile(width, 90)),
        "cols_used": int(ok.sum()),
    }


# ------------------------------------------------------------
# 线上能量占比
# ------------------------------------------------------------
def line_energy(img_lin: np.ndarray, row_logf: np.ndarray, gt_logf: np.ndarray,
                cols: np.ndarray, tol_cents: float = 50.0,
                spur_cents: float = 200.0) -> dict:
    """能量(功率)落在真值线附近的占比, 以及远离线的虚假能量占比。"""
    p = np.square(img_lin[:, cols].astype(np.float64))
    tot = p.sum()
    if tot <= 0:
        return {"on_line": np.nan, "spurious": np.nan}
    dist = np.abs(row_logf[:, None] - gt_logf[cols][None, :]) * LOG2_10 * 1200.0
    near = dist <= tol_cents
    far = dist >= spur_cents
    return {
        "on_line": float(p[near].sum() / tot),
        "spurious": float(p[far].sum() / tot),
    }


# ------------------------------------------------------------
# 霍夫变换找直线
# ------------------------------------------------------------
def hough_line(img_db: np.ndarray, thr_db: float = -30.0, slope_range: tuple[float, float] = (0.0, 2.0),
               n_slopes: int = 121, col_shift: int = 0,
               cols_mask: np.ndarray | None = None) -> dict:
    """对二值化的谱图做霍夫变换, 找最强直线。

    参数
    ----
    img_db    : (n_rows, n_cols) dB 图(0 dB = 标定后的单位纯音, 或相对自身峰值)
    thr_db    : 二值化阈值(dB)
    slope_range : 候选斜率范围(行/列)
    col_shift : 图的左侧垫列数(用于把列还原成绝对列)
    cols_mask : 只在 True 的列上统计(限定在真值线所在的区间, 避免带外伪影稀释指标)

    返回 ``dict``: ``slope``(行/列)、``support``(落在最强直线 ±1 行的像素占比)、
    ``perp_rms_rows``(这些像素到直线的垂直 RMS, 行)、``theta_deg``(相对列轴的夹角)。
    """
    idx_all = np.arange(img_db.shape[1]) if cols_mask is None else np.nonzero(cols_mask)[0]
    if cols_mask is not None:
        img_db = img_db[:, cols_mask]
    ys, xs = np.nonzero(img_db > thr_db)
    if ys.size < 10:
        return {"slope": np.nan, "support": np.nan, "perp_rms_rows": np.nan,
                "theta_deg": np.nan, "pixels": int(ys.size)}
    x = idx_all[xs] - col_shift
    slopes = np.linspace(*slope_range, n_slopes)
    best = (None, -1.0, np.nan)
    for m in slopes:
        b = ys - m * x
        hist, edges = np.histogram(b, bins=np.arange(b.min() - 1.5, b.max() + 2.5, 1.0))
        i = int(np.argmax(hist))
        if hist[i] > best[1]:
            best = (m, float(hist[i]), 0.5 * (edges[i] + edges[i + 1]))
    m, cnt, b0 = best
    res = ys - (m * x + b0)
    sel = np.abs(res) <= 1.0
    perp = res[sel] * np.cos(np.arctan(m))
    return {
        "slope": float(m),
        "support": float(cnt / ys.size),
        "perp_rms_rows": float(np.sqrt(np.mean(perp ** 2))) if perp.size else np.nan,
        "theta_deg": float(np.degrees(np.arctan(m))),
        "pixels": int(ys.size),
    }


# ------------------------------------------------------------
# 每列局部极大值个数(虚假谱线计数)
# ------------------------------------------------------------
def peak_count(img_db: np.ndarray, cols: np.ndarray, thr_db: float = -30.0) -> dict:
    """统计每列超过 ``thr_db`` 的局部极大值个数(相邻行取 1)。"""
    sub = img_db[:, cols]
    if sub.size == 0:
        return {"peaks_med": np.nan, "peaks_mean": np.nan}
    mid = sub[1:-1]
    is_peak = (mid > sub[:-2]) & (mid >= sub[2:]) & (sub[1:-1] > thr_db)
    n = is_peak.sum(axis=0)
    return {"peaks_med": float(np.median(n)), "peaks_mean": float(np.mean(n))}


# ------------------------------------------------------------
# 汇总
# ------------------------------------------------------------
def describe(img_lin: np.ndarray, img_db: np.ndarray, row_logf: np.ndarray,
             cfg, n_cols: int, col_shift: int, f_true, t_span: tuple[float, float],
             tol_cents: float = 50.0, peak_rel_db: float = 25.0,
             f_band: tuple[float, float] = (20.0, 20000.0)) -> dict:
    """对一张图算一整套指标(自动挑出真值落在 ``t_span`` 与 ``f_band`` 内的列)。

    直线性/峰数这类阈值按**相对自身峰值**给(``peak_rel_db`` 低于峰值多少 dB),
    避免不同方法/变体标定电平不同时阈值失效。``f_band`` 用于排除低频段:
    NC 的窗长被 ``max_window_s`` 夹住后带宽远宽于低频行, 该段属于已知的
    分辨率退化区, 单独统计更公平。
    """
    times = cfg.col_times(n_cols + col_shift, col_shift)
    rows = np.arange(n_cols + col_shift)
    valid = (times >= t_span[0]) & (times <= t_span[1])
    gt_f = f_true(times)
    valid &= (gt_f >= max(cfg.f_min * 1.02, f_band[0])) & (gt_f <= min(cfg.f_max * 0.98, f_band[1]))
    valid &= rows >= col_shift
    cols = rows[valid]
    gt_logf = np.log10(gt_f)
    img_rel = img_db - np.max(img_db)
    out = {}
    out.update(ridge_stats(img_lin, row_logf, gt_logf, cols))
    out.update(line_energy(img_lin, row_logf, gt_logf, cols, tol_cents=tol_cents))
    out.update(hough_line(img_rel, thr_db=-peak_rel_db, col_shift=col_shift,
                          cols_mask=valid))
    out.update(peak_count(img_rel, cols, thr_db=-peak_rel_db))
    # 理论斜率(行/列): d(log10 f)/dt · (hop/fs) / 行间距
    row_step = np.median(np.diff(row_logf))
    slope_theory = (np.log10(f_true(t_span[0] + 0.5)) - np.log10(f_true(t_span[0]))) / 0.5 \
        if t_span[1] > t_span[0] else np.nan
    out["slope_theory"] = float(slope_theory * (cfg.hop / cfg.fs) / row_step)
    # 垂直(垂直于 chirp 轨迹)展宽
    out["width_perp_oct_med"] = float(out["width_oct_med"] * np.cos(np.arctan(out["slope_theory"])))
    return out
