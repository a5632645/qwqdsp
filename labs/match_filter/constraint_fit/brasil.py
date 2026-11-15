"""BRASIL 区间等化：通用的节点位置搜索。

节点把 ``[0, 1]`` 切成 ``n + 1`` 段。``measure(nodes)`` 返回每段的
max|ΔdB|；节点不可行时返回 ``None``。每轮按相对偏差给出目标段长，回溯
步长，接受让区间极值比或最大误差严格下降的候选，并保留 max 最小的解。
"""
from __future__ import annotations

import numpy as np


def segment_maxima(xgrid, err, nodes):
    """按 ``nodes`` 切段，返回每段的最大 ``err``。

    ``xgrid`` 必须单调递增，``err`` 与它等长。
    """
    edges = np.concatenate(([0.0], np.sort(np.asarray(nodes, float)), [1.0]))
    idx = np.searchsorted(xgrid, edges)
    out = np.empty(edges.size - 1)
    for i in range(out.size):
        lo, hi = idx[i], min(idx[i + 1] + 1, xgrid.size)
        out[i] = err[lo:hi].max() if hi > lo else 0.0
    return out


def equalize(nodes0, measure, *, max_iter=1000, step_max=0.15, tol=1e-6):
    """调整内部节点，使各段最大误差等幅。

    返回 ``(nodes, maxima, best_max)``；无可行初值返回
    ``(None, None, inf)``。
    """
    nodes = np.sort(np.asarray(nodes0, float))
    if nodes.size == 0:
        return None, None, np.inf
    edges0 = np.concatenate(([0.0], nodes, [1.0]))
    if np.min(np.diff(edges0)) <= 1e-12:
        return None, None, np.inf
    errs = measure(nodes)
    if errs is None:
        return None, None, np.inf
    errs = np.asarray(errs, float)
    if errs.size != nodes.size + 1 or not np.all(np.isfinite(errs)):
        return None, None, np.inf
    best = (float(errs.max()), nodes.copy(), errs.copy())
    best_spread = best[0] / max(float(errs.min()), 1e-300)
    stall = 0

    for _ in range(max_iter):
        top = float(errs.max())
        spread = top / max(float(errs.min()), 1e-300)
        if spread - 1.0 < tol:
            break
        mean = float(errs.mean())
        dev = errs - mean
        dmax = float(np.abs(dev).max())
        if dmax <= 0.0:
            break
        base0 = min(step_max, dmax / mean) if mean > 0.0 else step_max
        edges = np.concatenate(([0.0], nodes, [1.0]))
        improved = False
        for bt in range(24):
            base = base0 * 0.5 ** bt
            lengths = (1.0 + base) ** (-dev / dmax) * np.diff(edges)
            lengths /= lengths.sum()
            cand = np.cumsum(lengths)[:-1]
            if np.min(np.diff(np.concatenate(([0.0], cand, [1.0])))) <= 1e-12:
                continue
            cand_errs = measure(cand)
            if cand_errs is None:
                continue
            cand_errs = np.asarray(cand_errs, float)
            if (cand_errs.size != nodes.size + 1
                    or not np.all(np.isfinite(cand_errs))):
                continue
            cand_spread = cand_errs.max() / max(float(cand_errs.min()), 1e-300)
            if (cand_spread < spread - 1e-12
                    or cand_errs.max() < top - 1e-12):
                nodes, errs = cand, cand_errs
                improved = True
                break
        if not improved:
            break
        top_new = float(errs.max())
        spread_new = top_new / max(float(errs.min()), 1e-300)
        if (top_new < best[0] - 1e-12
                or (abs(top_new - best[0]) <= 1e-12
                    and spread_new < best_spread)):
            best = (top_new, nodes.copy(), errs.copy())
            best_spread = spread_new
        if top_new < top - 1e-12:
            stall = 0
        else:
            stall += 1
            if stall > 500:
                break
    return best[1], best[2], best[0]
