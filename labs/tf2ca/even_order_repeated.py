# -*- coding: utf-8 -*-
"""偶数阶、分母含重极点的设计（典型 Linkwitz–Riley）的双链判据与三链分解。

    LR_{2n} = (BW_n)^2

穷举 ``A*A^#`` 的全部谱因子 ``U``，用最小二乘检验 ``B`` 是否落在
``span{U + U^#, j*(U - U^#)}`` 内；LR 不满足该条件，正确结构是三条链
（权重 ``1/4, 1/2, 1/4``）。重极点下逐极点 ``alpha`` 判据会假阳性，脚本改用
系数空间残差，并同时展示该失效。

用法::

    python even_order_repeated.py                                  # 默认 LR4（半阶 2）, cutoff=0.25
    python even_order_repeated.py --half-order 4 --cutoff 0.35     # LR8
    python even_order_repeated.py --delays                         # 额外尝试 z^-k 延迟偏移
    python even_order_repeated.py --check                          # LR2..LR12 × 多截止 + 对照
"""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 无界面后端
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.signal import lfilter

from even_order import real_two_path_min_residual
from tf2ca import complex_allpass_response, conj_rev, design, split_by_radius, zfreq

SCRIPT_DIR = Path(__file__).resolve().parent


# ------------------------------------------------------------
# 复全通、完全平方与重复极点结构
# ------------------------------------------------------------


def design_lr(half_order: int, cutoff: float) -> tuple[np.ndarray, np.ndarray]:
    """设计 Linkwitz–Riley 低通 ``LR_{2n} = (n 阶 Butterworth)^2``，返回 (b, a)。"""
    b, a = signal.butter(half_order, cutoff, "low")
    return np.convolve(b, b), np.convolve(a, a)


def hp_partner(half_order: int, cutoff: float) -> np.ndarray:
    """LR 的高通分子（与低通同分母）：``LR 高通 = (n 阶 Butterworth 高通)^2``。"""
    b, _ = signal.butter(half_order, cutoff, "high")
    return np.convolve(b, b)


def lr_half_denominator(half_order: int, cutoff: float) -> np.ndarray:
    """``LR_{2n}`` 分母的平方根方向：``a = (半阶 BW 分母)^2`` 的那个半阶分母。

    重极点设计必须把这条**解析已知**的结构告诉判据；想靠数值手段判断「A 是不是完全平方」
    是不可靠的：butter 8 阶 @0.05 这种窄带高阶设计，A 根本不是平方，但系数级平方根递推的
    相对残差只有 2e-7（真平方是 1e-13 以下），启发式阈值会误判。
    """
    _, a_half = signal.butter(half_order, cutoff, "low")
    return a_half


def cluster_roots(poly: np.ndarray, tol: float = 1e-7):
    """把 ``poly`` 的根按共轭类归并，返回 ``[(rep, nu), ...]``（按 |p| 升序）。

    ``rep`` 取上半平面代表（实根取自身），``nu`` 为该共轭类在 ``poly`` 中的重数
    （复共轭对 ``(p, p*)`` 各出现 nu 次）。
    """
    poles = np.roots(np.asarray(poly, dtype=float))
    used = np.zeros(poles.size, bool)
    out: list[tuple[complex, int]] = []
    for i, p in enumerate(poles):
        if used[i]:
            continue
        rep = p if p.imag >= 0 else np.conjugate(p)
        idx = [j for j in range(poles.size)
               if not used[j] and (abs(poles[j] - rep) < tol or abs(poles[j] - rep.conjugate()) < tol)]
        used[idx] = True
        out.append((complex(rep), len(idx) if abs(rep.imag) < tol else len(idx) // 2))
    out.sort(key=lambda t: (abs(t[0]), -t[0].real))
    return out


def pole_classes(a: np.ndarray, tol: float = 1e-7, sqrt_hint: np.ndarray | None = None):
    """A 的极点结构 ``[(rep, nu), ...]``（``nu`` = 该代表根在 A 中的重数）。

    ``sqrt_hint`` 给出 ``A = (sqrt_hint)^2`` 的那个平方根多项式（LR 的情形，解析已知）：
    此时每个极点的重数是半阶设计的 2 倍。不给 hint 时直接对 A 聚类——这对「极点互异」
    的设计（四族）正确，对重极点设计会漏掉重数信息（脚本内部一律传 hint）。
    """
    if sqrt_hint is not None:
        return [(rep, 2 * nu) for rep, nu in cluster_roots(sqrt_hint, tol)]
    return cluster_roots(a, tol)


def describe_poles(a: np.ndarray, sqrt_hint: np.ndarray | None = None) -> str:
    """把 A 的极点结构压成一行（``|p|=0.5774(±0.3333j)×2``）。"""
    return "  ".join(
        f"|p|={abs(rep):.4f}({'实' if abs(rep.imag) < 1e-7 else f'±{abs(rep.imag):.4f}j'})×{nu}"
        for rep, nu in pole_classes(a, sqrt_hint=sqrt_hint)
    )


# ------------------------------------------------------------
# 判据：穷举全部谱因子 U
# ------------------------------------------------------------


def spectral_factors(a: np.ndarray, tol: float = 1e-7, sqrt_hint: np.ndarray | None = None):
    """穷举 A*A^# 的**全部**谱因子 U（``deg U = n`` 且 ``U * U^# = mu * A * A^#``）。

    谱因子的定义是「每个互反共轭根对取一个根」，按共轭类：

    - 复共轭类（代表 p、在 A 中的重数 nu）：根对有两族——``{p, 1/p*}`` 与 ``{p*, 1/p}``
      各 nu 份，每份二选一，于是 ``(nu+1)^2`` 种取法 ``(k, j)``（k 取 ``1/p*`` 的份数、
      j 取 ``1/p`` 的份数），U 的根为
      ``p^(nu-k) * (1/p*)^k * p*^(nu-j) * (1/p)^j``；
      对应 ``D1 = p^k * p*^j``、``D2 = p^(nu-k) * p*^(nu-j)``（满足 ``U = N1*D2``）。
    - 实根（重数 nu）：只有一族根对 ``{r, 1/r}`` 的 nu 份，``nu+1`` 种取法 ``m``（取 ``1/r``
      的份数），U 的根为 ``r^(nu-m) * (1/r)^m``；对应 ``D1 = r^m``、``D2 = r^(nu-m)``。

    极点互异情形的枚举只覆盖复类的对角线 ``j = nu - k``（等价于 ``D2 = conj(D1)``）；
    重极点下必须放开到全部 ``(nu+1)^2``（以及实根的 nu+1），否则会漏掉合法的谱因子。

    :return: ``[(U, D1, D2, 选择元组), ...]``；``U`` 与 ``U^#`` 给出同一个 2 维子空间
        （残差相同），故只枚举互补配对的一半。
    """
    groups = pole_classes(a, tol, sqrt_hint)

    def poly_of(roots):
        return np.poly(roots) if len(roots) else np.array([1.0])

    per_class = []
    for rep, nu in groups:
        entries = []
        if abs(rep.imag) < tol:  # 实根：根对 {r, 1/r}
            for m in range(nu + 1):
                entries.append(((m,), (nu - m,),
                                poly_of([rep] * (nu - m) + [1.0 / rep] * m),
                                poly_of([rep] * m), poly_of([rep] * (nu - m))))
        else:  # 复共轭类：根对 {p, 1/p*} 与 {p*, 1/p}
            conj, inv_rep, inv_conj = np.conjugate(rep), 1.0 / rep, 1.0 / np.conjugate(rep)
            for k in range(nu + 1):
                for j in range(nu + 1):
                    entries.append(((k, j), (nu - k, nu - j),
                                    poly_of([rep] * (nu - k) + [inv_conj] * k
                                            + [conj] * (nu - j) + [inv_rep] * j),
                                    poly_of([rep] * k + [conj] * j),
                                    poly_of([rep] * (nu - k) + [conj] * (nu - j))))
        per_class.append(entries)

    out = []
    aa = np.convolve(np.asarray(a, dtype=float), conj_rev(a))
    for combo in itertools.product(*per_class):
        sel = tuple(e[0] for e in combo)
        mate = tuple(e[1] for e in combo)
        if sel > mate:  # U 与 U^# 等价，只留字典序较小的一半
            continue
        u = np.array([1.0], complex)
        d1 = np.array([1.0], complex)
        d2 = np.array([1.0], complex)
        for entry in combo:
            u = np.convolve(u, entry[2])
            d1 = np.convolve(d1, entry[3])
            d2 = np.convolve(d2, entry[4])
        # 归一化：U*U^# 本来就是 A*A^# 的正实数倍，缩放到严格相等。
        # 只有这个归一化下「|alpha| = 0.5」才等价于「导出的互补支与低通功率互补」。
        gamma = np.sqrt(abs(aa[0]) / abs(np.convolve(u, conj_rev(u))[0]))
        u = u * gamma
        d1 = np.real_if_close(d1 * gamma)  # 保持 U = N1*D2（D1 同比缩放）
        out.append((u, d1, np.real_if_close(d2), sel))
    return out


def two_chain_fit(b: np.ndarray, a: np.ndarray, shifts=(0,), tol: float = 1e-7, keep: int = 8,
                  sqrt_hint: np.ndarray | None = None):
    """对每个谱因子（× 延迟偏移）做最小二乘，返回全部残差与最优解。

    无延迟（``j1 = j2 = 0``）时，两条链互为系数共轭 ⇒ 权重必为共轭对，
    于是待解的是 2 个实未知量 ``x, y``：``B ≈ x*(U + U^#) + y*j*(U - U^#)``
    = ``alpha*U + conj(alpha)*U^#``（``alpha = x + j*y``）。

    两处都用 SVD 最小二乘而不是手写正规方程：``U`` 与 ``U^#`` 近乎共线时（例如
    ``D1 = D2`` 的退化候选）手写解会给出巨大的 ``x, y``，残差在浮点下被抵消成「假精确」
    （实测 LR4@0.02 会出现 3e-16 的假阳性，而 60 位精度下真实残差是 0.93）。
    SVD 路径会自动截断退化方向，残差如实反映可实现的拟合精度。

    带延迟时链条形如 ``alpha*z^-j1*N1/D1 + beta*z^-j2*N2/D2``，两项不再互为共轭，
    改用 2 个**复**未知量（= 4 个实未知量）。

    :return: ``dict``，含
        ``best``        全局最优行 ``(D1, D2, counts, residual, alpha, (j1, j2))``；
        ``residuals``   无延迟全部候选的残差数组（与 :func:`spectral_factors` 顺序一致）；
        ``sample``      无延迟残差最小的 ``keep`` 行；
        ``delay_best``  带延迟时的最优行（``shifts`` 只有 (0,0) 时为 ``None``）；
        ``n_candidates``无延迟候选数。
    """
    b = np.asarray(b, dtype=float)
    scale = float(np.max(np.abs(b)))
    cands = spectral_factors(a, tol, sqrt_hint)
    rows: list[tuple] = []
    residuals = np.empty(len(cands))
    for idx, (u, d1, d2, counts) in enumerate(cands):
        # 无延迟：权重为共轭对，2 个实未知量（用 SVD 最小二乘，不用手写正规方程：
        # U 与 U^# 近乎共线时手写解会给出「假精确」的残差，见 docstring）
        n = max(u.size, b.size)
        cu = np.zeros(n, complex)
        cu[:u.size] = u
        m_ = conj_rev(u)
        cm = np.zeros(n, complex)
        cm[:m_.size] = m_
        col_u = (cu + cm).real                # U + U^#
        col_v = (1j * (cu - cm)).real         # j*(U - U^#)
        m = np.stack([col_u[:b.size], col_v[:b.size]], axis=1)
        sol, *_ = np.linalg.lstsq(m, b, rcond=None)
        delta = b - m @ sol
        residuals[idx] = float(np.max(np.abs(delta))) / scale
        rows.append((d1, d2, counts, float(residuals[idx]), sol[0] + 1j * sol[1], (0, 0)))
    order = np.argsort(residuals)[:keep]
    best = dict(zip(("d1", "d2", "counts", "residual", "alpha", "shift"),
                    rows[int(np.argmin(residuals))]))

    delay_best = None
    for j1 in shifts:
        for j2 in shifts:
            if (j1, j2) == (0, 0):
                continue
            for (u, d1, d2, counts) in cands:
                basis = (np.r_[np.zeros(j1, complex), u], np.r_[np.zeros(j2, complex), conj_rev(u)])
                n = max(basis[0].size, basis[1].size, b.size)
                m = np.zeros((n, 2), complex)
                for k, vec in enumerate(basis):
                    m[:vec.size, k] = vec
                rhs = np.r_[b, np.zeros(n - b.size)]
                sol, *_ = np.linalg.lstsq(m, rhs, rcond=None)
                delta = m @ sol - rhs
                res = float(np.max(np.abs(delta)) / scale)
                if delay_best is None or res < delay_best[3]:
                    delay_best = (d1, d2, counts, res, sol[0], (j1, j2))

    return {"best": best, "residuals": residuals,
            "sample": [rows[int(i)] for i in order],
            "delay_best": delay_best, "n_candidates": len(cands),
            "n_shifts": len(set(shifts))}


# ------------------------------------------------------------
# 重极点下逐极点判据的失效
# ------------------------------------------------------------


def naive_pole_spread(b: np.ndarray, d1: np.ndarray, d2: np.ndarray):
    """极点互异情形使用的逐极点判据：在 A 的每个极点独立解 alpha，返回相对离散度。

    判据只检查**值**条件，且默认「每个极点只属于一条链」。于是它对重极点有两重失效：
    ① D1 与 D2 共享极点时（``D1 = D2 = (z-p)``）分母为 0，判据直接无定义（返回 nan）；
    ② 即使它能算出值，重极点下值条件不足以定住恒等式，见 :func:`identity_multiple`。
    """
    n1, n2 = conj_rev(d1), conj_rev(d2)
    with np.errstate(divide="ignore", invalid="ignore"):
        vals = [np.polyval(b, p) / (np.polyval(n1, p) * np.polyval(d2, p)) for p in np.roots(d1)]
        vals += [np.conj(np.polyval(b, p) / (np.polyval(n2, p) * np.polyval(d1, p))) for p in np.roots(d2)]
    vals = np.asarray([v for v in vals if np.isfinite(v)])
    if vals.size == 0:
        return float("nan"), vals
    mean = vals.mean()
    return float(np.max(np.abs(vals - mean)) / abs(mean)), vals


def identity_multiple(b: np.ndarray, d1: np.ndarray, d2: np.ndarray):
    """重极点下逐极点判据的失效机制：``B - (a*N1*D2 + conj(a)*N2*D1) = c * A``。

    用逐极点判据解出的 ``alpha`` 代入，残差多项式在 A 的每个极点都取 0，但它不是零多项式
    ——它等于 ``c * A``（A 只有 n/2 个相异根，值条件定不住次数 n 的恒等式），
    于是残差量级是 ``|c| * max|A| / max|B|``，而逐极点判据完全看不见。

    :return: ``(c, rel_residual, spread)``。
    """
    spread, vals = naive_pole_spread(b, d1, d2)
    if vals.size == 0:
        return complex(np.nan), float("nan"), float("nan")
    alpha = np.mean(vals)
    resid = np.polyadd(
        np.asarray(b, dtype=float),
        -(alpha * np.convolve(conj_rev(d1), d2) + np.conjugate(alpha) * np.convolve(conj_rev(d2), d1)),
    )
    c, _ = np.polydiv(resid, np.convolve(d1, d2))
    return complex(np.mean(c)), float(np.max(np.abs(resid)) / np.max(np.abs(b))), spread


# ------------------------------------------------------------
# 三链结构：LR = (BW)^2 的精确分解
# ------------------------------------------------------------


def lr_three_chain(half_order: int, cutoff: float):
    """``LR_{2n}`` 的精确**三链**分解，返回 ``[(denom, 链内单位缩放, 权重), ...]``。

    半阶 n 为奇数：BW_n 是实系数双链 ``1/2 (A0 + A1)``（`tf2ca.split_by_radius`），
    平方得 ``1/4 A0^2 + 1/2 A0*A1 + 1/4 A1^2``（三条实链）。
    半阶 n 为偶数：BW_n 是复系数双链 ``alpha*A1 + conj(alpha)*A2``（`even_order`），
    实测 ``|alpha| = 0.5``，平方得 ``1/4 (2alpha)^2 A1^2 + 1/2 (A1*A2)
    + 1/4 (2conj(alpha))^2 A2^2``（外侧链复系数、互为系数共轭，中间链实系数）。
    """
    b_half, a_half = signal.butter(half_order, cutoff, "low")
    if half_order % 2 == 0:
        from even_order import decompose

        dc = decompose(b_half, a_half)
        d1, d2, alpha = dc["d1"], dc["d2"], dc["alpha"]
        triples = [(np.convolve(d1, d1), (2 * alpha) ** 2), (np.convolve(d1, d2), 1.0),
                   (np.convolve(d2, d2), (2 * np.conjugate(alpha)) ** 2)]
    else:
        d0, d1 = split_by_radius(a_half)
        triples = [(np.convolve(d0, d0), 1.0), (np.convolve(d0, d1), 1.0),
                   (np.convolve(d1, d1), 1.0)]
    weights = (0.25, 0.5, 0.25)
    return [(d, s, w) for (d, s), w in zip(triples, weights)]


def three_chain_metrics(b: np.ndarray, a: np.ndarray, chains, num: int = 4000, terms: int = 80):
    """校验三链分解：频响误差、冲激响应误差（含虚部泄漏）、链平坦度、``H(1)``。"""
    w = np.linspace(1e-6, np.pi - 1e-6, num)
    got = sum(wt * sc * complex_allpass_response(d, w) for d, sc, wt in chains)
    ref = zfreq(b, w) / zfreq(a, w)
    x = np.zeros(terms)
    x[0] = 1.0
    y_struct = sum(wt * sc * lfilter(conj_rev(d), d, x) for d, sc, wt in chains)
    y_direct = lfilter(b, a, x)
    return {
        "freq": float(np.max(np.abs(got - ref))),
        "imp": float(np.max(np.abs(y_struct - y_direct))),
        "imag": float(np.max(np.abs(y_struct.imag))),
        "flat": float(max(np.max(np.abs(np.abs(complex_allpass_response(d, w)) - 1)) for d, _, _ in chains)),
        "h1": float(abs(sum(wt * sc * complex_allpass_response(d, np.array([1e-9])) for d, sc, wt in chains)[0])),
        "orders": [len(d) - 1 for d, _, _ in chains],
    }


def allpass_defect(b_sum: np.ndarray, a: np.ndarray, num: int = 4000) -> float:
    """``|(b_sum/A)(e^jw)|`` 偏离 1 的最大量（0 表示该组合是全通）。"""
    w = np.linspace(1e-6, np.pi - 1e-6, num)
    return float(np.max(np.abs(np.abs(zfreq(b_sum, w) / zfreq(a, w)) - 1.0)))


# ------------------------------------------------------------
# 报告
# ------------------------------------------------------------


def report(b: np.ndarray, a: np.ndarray, fit: dict, chains: list, cm: dict,
           half_order: int, cutoff: float, sqrt_hint: np.ndarray) -> None:
    """打印结构、双链穷举结果、重极点假阳性、三链校验与同阶对照。"""
    best = fit["best"]
    d1, d2, counts, res, alpha, sh = (best["d1"], best["d2"], best["counts"],
                                      best["residual"], best["alpha"], best["shift"])
    b_half, a_half = signal.butter(half_order, cutoff, "low")
    print(f"  极点结构（{len(pole_classes(a, sqrt_hint=sqrt_hint))} 个共轭类，A 阶数 {a.size - 1}）：{describe_poles(a, sqrt_hint)}")
    print(f"  B/BW 镜像对称 {np.allclose(b, b[::-1], atol=1e-12)}"
          f"，B = (BW{half_order} 分子)^2 {np.allclose(b, np.convolve(b_half, b_half))}"
          f"，A = (BW{half_order} 分母)^2 {np.allclose(a, np.convolve(a_half, a_half))}")
    print()
    print(f"== 双链穷举（{fit['n_candidates']} 个谱因子，已按 U <-> U^# 去重"
          f"{'；每种再试 ' + str(fit['n_shifts']) + ' 种延迟偏移' if fit['delay_best'] is not None else ''}）==")
    print(f"  最优分法：residual = {res:.3e}   alpha = {alpha:.6g}  |alpha| = {abs(alpha):.6f}"
          f"   deg D1 = {d1.size - 1}, deg D2 = {d2.size - 1}   延迟偏移 = {sh}")
    print("  判据：B = (U + U^#)/2 且 U*U^# = A*A^#（共轭基下后者正是 |alpha| = 0.5）；"
          "四族偶数阶两者都满足，LR 连残差都下不去")
    print(f"  候选残差分布：min = {np.min(fit['residuals']):.3e}  "
          f"median = {np.median(fit['residuals']):.3e}  max = {np.max(fit['residuals']):.3e}")
    if fit["delay_best"] is not None:
        db = fit["delay_best"]
        print(f"  带延迟最优：residual = {db[3]:.3e}  延迟偏移 = {db[5]}"
              f"（延迟只改分母/分子对齐，救不回来）")

    c, rel, spread = identity_multiple(b, d1, d2)
    print()
    print("== 重极点下「逐极点 alpha 判据」为什么不可信 ==")
    print(f"  逐极点 alpha 离散度（even_order.py 的主判据） = {spread:.3e}   <- 看起来「一致」")
    print(f"  但同一组 alpha 代入后 B - (a*N1*D2 + conj(a)*N2*D1) = ({c:.6g}) * A"
          f"（相对残差 {rel:.3e}）")
    print(f"  A 有 {a.size - 1} 阶却只有 {len(pole_classes(a, sqrt_hint=sqrt_hint))} 个相异极点：值条件定不住恒等式，"
          f"判分法必须看系数空间残差")

    nyq = abs(zfreq(b, np.array([np.pi])) / zfreq(a, np.array([np.pi])))[0]
    print()
    print("== 实系数两路全通（1/2 (A0 + A1)）对照 ==")
    print(f"  最佳 B 残差 = {real_two_path_min_residual(b, a):.3e}    |H(-1)| = {nyq:.3e}")
    print(f"  证明：n 偶数 ⇒ 两链阶数同奇偶 ⇒ |H(-1)| = 1；LR 低通 |H(-1)| = {nyq:.1e} 故不可能"
          f"（加 z^-k 也不救：延迟只翻转符号）")

    print()
    print("== LR 的正确结构：三条链 ==")
    for k, (d, sc, wt) in enumerate(chains):
        kind = "复" if abs(complex(sc).imag) > 1e-12 or np.iscomplexobj(d) else "实"
        print(f"  链{k}: 阶 {len(d) - 1:2d}  权重 {wt:.2f}  链内单位缩放 {sc:.6g}（{kind}）")
    print(f"  频响最大偏差 = {cm['freq']:.3e}   冲激响应最大偏差 = {cm['imp']:.3e}"
          f"（虚部泄漏 {cm['imag']:.1e}）")
    print(f"  三链全通平坦度 = {cm['flat']:.3e}   H(1) = {cm['h1']:.9f}   链阶数 = {cm['orders']}")

    hbp = hp_partner(half_order, cutoff)
    print()
    print("== LR 的 LP / HP：和是全通 ≠ 能拆成耦合全通 ==")
    print(f"  |LP + HP| - 1 最大偏差 = {allpass_defect(b + hbp, a):.3e}")
    print(f"  |LP - HP| - 1 最大偏差 = {allpass_defect(b - hbp, a):.3e}")
    print("  （偶数半阶 LP+HP 全通、奇数半阶 LP-HP 全通；另一支不是全通，双链判据仍不满足）")

    print()
    print("== 同阶四族对照（同一套穷举代码）==")
    for kind in ("butter", "cheby1", "cheby2", "ellip"):
        bb, aa = design(kind, 2 * half_order, cutoff)
        cfit = two_chain_fit(bb, aa, keep=1)
        print(f"  {kind:8s} 阶 {2 * half_order:2d}：双链最优残差 = {cfit['best']['residual']:.3e}"
              f"  |alpha| = {abs(cfit['best']['alpha']):.6f}")


def run_check() -> None:
    """批量表：LR2..LR12 × 多截止频率 + 四族对照 + 四族「平方」对照 + 重极点构造。"""
    print("== LR（偶数阶 + 重极点）双链穷举（--check）==")
    print(f"{'LR':>5s}{'阶':>4s}{'截止':>7s}{'双链最优残差':>14s}{'朴素离散度':>12s}"
          f"{'实系数双链':>12s}{'|H(-1)|':>10s}{'三链误差':>11s}{'谱因子数':>10s}")
    for half_order in (1, 2, 3, 4, 5, 6):
        for cutoff in (0.1, 0.25, 0.45):
            b, a = design_lr(half_order, cutoff)
            fit = two_chain_fit(b, a, keep=1, sqrt_hint=lr_half_denominator(half_order, cutoff))
            spread, _ = naive_pole_spread(b, fit["best"]["d1"], fit["best"]["d2"])
            cm = three_chain_metrics(b, a, lr_three_chain(half_order, cutoff), num=2000, terms=60)
            nyq = abs(zfreq(b, np.array([np.pi])) / zfreq(a, np.array([np.pi])))[0]
            print(f"{f'LR{2 * half_order}':>5s}{2 * half_order:>4d}{cutoff:>7.2f}"
                  f"{fit['best']['residual']:>14.3e}{spread:>12.3e}"
                  f"{real_two_path_min_residual(b, a):>12.3e}{nyq:>10.1e}"
                  f"{cm['freq']:>11.3e}{fit['n_candidates']:>10d}")
    print()
    print("== 对照 A：四族偶数阶（无重极点）——同一套穷举代码应达机器精度 ==")
    print(f"{'族':10s}{'阶':>4s}{'截止':>7s}{'双链最优残差':>14s}{'朴素离散度':>12s}"
          f"{'|alpha|':>10s}{'谱因子数':>10s}{'实系数双链':>12s}")
    for kind in ("butter", "cheby1", "cheby2", "ellip"):
        for order in (2, 4, 6):
            for cutoff in (0.1, 0.3):
                b, a = design(kind, order, cutoff)
                fit = two_chain_fit(b, a, keep=1)
                spread, _ = naive_pole_spread(b, fit["best"]["d1"], fit["best"]["d2"])
                print(f"{kind:10s}{order:>4d}{cutoff:>7.2f}{fit['best']['residual']:>14.3e}"
                      f"{spread:>12.3e}{abs(fit['best']['alpha']):>10.6f}"
                      f"{fit['n_candidates']:>10d}{real_two_path_min_residual(b, a):>12.3e}")
    print()
    print("== 对照 B：四族的「平方」（同样是重极点设计）也做不到双链 ==")
    print(f"{'设计':14s}{'阶':>5s}{'双链最优残差':>14s}{'实系数双链':>12s}{'|H(-1)|':>10s}")
    for kind, order in (("butter", 2), ("cheby1", 2), ("cheby2", 2), ("ellip", 2),
                        ("butter", 4), ("ellip", 4)):
        b0, a0 = design(kind, order, 0.3)
        b, a = np.convolve(b0, b0), np.convolve(a0, a0)
        fit = two_chain_fit(b, a, keep=1, sqrt_hint=a0)
        nyq = abs(zfreq(b, np.array([np.pi])) / zfreq(a, np.array([np.pi])))[0]
        print(f"{f'({kind}{order})^2':14s}{2 * order:>5d}{fit['best']['residual']:>14.3e}"
              f"{real_two_path_min_residual(b, a):>12.3e}{nyq:>10.1e}")
    print()
    print("== 对照 C：重极点本身不是障碍（同一个 A_LR4，存在精确可双链的 B）==")
    _, a_lr = design_lr(2, 0.25)
    hint_lr = lr_half_denominator(2, 0.25)
    seen = set()
    for u, _, _, _ in spectral_factors(a_lr, sqrt_hint=hint_lr):
        b_alt = np.real(np.real_if_close((u + conj_rev(u)) / 2.0))
        key = tuple(np.round(b_alt, 9))
        if key in seen:
            continue
        seen.add(key)
        fit = two_chain_fit(b_alt, a_lr, keep=1, sqrt_hint=hint_lr)
        print(f"  B = (U + U^#)/2 = {np.array2string(b_alt, precision=6)}  残差 = "
              f"{fit['best']['residual']:.3e}  |H(1)| = "
              f"{abs(zfreq(b_alt, np.array([1e-9])) / zfreq(a_lr, np.array([1e-9])))[0]:.4f}")
    print("  -> 重极点并不让可双链的 B 变空，LR 的 B 恰好落在这些子空间之外"
          "（A、B 同时是 BW 的平方）")


# ------------------------------------------------------------
# 绘图
# ------------------------------------------------------------


def plot(b: np.ndarray, a: np.ndarray, fit: dict, chains, cm: dict,
         half_order: int, cutoff: float, out: Path) -> None:
    """六张图：双链对不上/三链对得上、重极点与谱因子、判据假阳性、冲激、阶数扫描、LP±HP。"""
    best = fit["best"]
    d1, d2, alpha, res = best["d1"], best["d2"], best["alpha"], best["residual"]
    w = np.linspace(1e-6, np.pi - 1e-6, 4000)
    ref = zfreq(b, w) / zfreq(a, w)
    got2 = alpha * complex_allpass_response(d1, w) + np.conjugate(alpha) * complex_allpass_response(d2, w)
    got3 = sum(wt * sc * complex_allpass_response(d, w) for d, sc, wt in chains)

    fig, ax = plt.subplots(2, 3, figsize=(16, 8.5))

    # (0,0) 幅频：直接型 vs 三链 vs 最优双链
    ax[0, 0].plot(w / np.pi, 20 * np.log10(np.abs(ref) + 1e-18), "k", lw=2, label="direct B/A")
    ax[0, 0].plot(w / np.pi, 20 * np.log10(np.abs(got3) + 1e-18), "r--", lw=1.2,
                  label=f"three chains (err {cm['freq']:.1e})")
    ax[0, 0].plot(w / np.pi, 20 * np.log10(np.abs(got2) + 1e-18), "b:", lw=1.2,
                  label=f"best two-chain (res {res:.2f})")
    ax[0, 0].set(title=f"LR{2 * half_order}: two chains fail, three match",
                 xlabel="$\\omega/\\pi$", ylabel="dB", ylim=(-100, 3))
    ax[0, 0].grid(alpha=0.3)
    ax[0, 0].legend(fontsize=8)

    # (0,1) z 平面：A 的重极点、被选中的谱因子 U
    theta = np.linspace(0, 2 * np.pi, 400)
    ax[0, 1].plot(np.cos(theta), np.sin(theta), "k:", lw=0.8)
    pa = np.roots(a)
    ax[0, 1].plot(pa.real, pa.imag, "x", color="#1f77b4", ms=9,
                  label="A poles (each with multiplicity 2)")
    pu = np.roots(conj_rev(np.convolve(conj_rev(d1), d2)))
    ax[0, 1].plot(pu.real, pu.imag, "o", mfc="none", mec="#d62728", ms=10,
                  label="roots of chosen U (one per pair)")
    mirrors = 1.0 / np.conjugate(pu)
    ax[0, 1].plot(mirrors.real, mirrors.imag, "+", color="gray", ms=8,
                  label="mirror roots (the other end of each pair)")
    for r, mr in zip(pu, mirrors):
        ax[0, 1].plot([r.real, mr.real], [r.imag, mr.imag], color="gray", lw=0.6, alpha=0.7)
    for p in pa:
        if p.imag > 1e-7:
            ax[0, 1].plot([p.real, p.real], [p.imag, -p.imag], color="gray", lw=0.6, alpha=0.6)
    ax[0, 1].axhline(0, color="gray", lw=0.5)
    ax[0, 1].axvline(0, color="gray", lw=0.5)
    ax[0, 1].set(title="repeated poles + spectral factor U: need B = (U+U$^\\#$)/2",
                 xlabel="Re", ylabel="Im", aspect="equal")
    ax[0, 1].grid(alpha=0.3)
    ax[0, 1].legend(loc="lower left", fontsize=8)

    # (0,2) 逐极点判据 vs 系数空间残差
    rs = np.sort(fit["residuals"])
    idx = np.arange(rs.size)
    ax[0, 2].semilogy(idx, np.maximum(rs, 1e-18), "b-", lw=1,
                      label="coeff-space residual (truth)")
    spread = [naive_pole_spread(b, s[0], s[1])[0] for s in fit["sample"]]
    spread = [s for s in spread if np.isfinite(s)]
    ax[0, 2].semilogy(np.arange(len(spread)), [max(s, 1e-18) for s in spread], "rx",
                      label="per-pole spread (naive, best few)")
    ax[0, 2].set(title=f"all {fit['n_candidates']} spectral splits, residual sorted",
                 xlabel="split (sorted by residual)", ylabel="relative")
    ax[0, 2].grid(alpha=0.3, which="both")
    ax[0, 2].legend(fontsize=8)

    # (1,0) 冲激响应
    x = np.zeros(80)
    x[0] = 1.0
    ax[1, 0].plot(lfilter(b, a, x), "k", lw=2, label="direct")
    ax[1, 0].plot(np.real(sum(wt * sc * lfilter(conj_rev(d), d, x) for d, sc, wt in chains)),
                  "r--", label="three chains")
    ax[1, 0].plot(np.real(alpha * lfilter(conj_rev(d1), d1, x)
                          + np.conjugate(alpha) * lfilter(conj_rev(d2), d2, x)), "b:",
                  label="best two-chain")
    ax[1, 0].set(title=f"impulse response (three-chain max|$\\Delta$|={cm['imp']:.1e})", xlabel="n")
    ax[1, 0].grid(alpha=0.3)
    ax[1, 0].legend(fontsize=8)

    # (1,1) 阶数扫描：LR vs 四族
    orders = [2, 4, 6, 8]
    ax[1, 1].semilogy(orders, [two_chain_fit(*design_lr(o // 2, cutoff), keep=1,
                                             sqrt_hint=lr_half_denominator(o // 2, cutoff))["best"]["residual"]
                               for o in orders], "ro-", label="LR (repeated poles)")
    for kind, style in (("butter", "s--"), ("cheby1", "^--"), ("cheby2", "v--"), ("ellip", "d--")):
        ax[1, 1].semilogy(orders, [two_chain_fit(*design(kind, o, cutoff), keep=1)["best"]["residual"]
                                   for o in orders], style, label=kind)
    ax[1, 1].set(title=f"two-chain residual vs order (cutoff={cutoff:g})", xlabel="order (even)",
                 ylabel="rel. residual")
    ax[1, 1].grid(alpha=0.3, which="both")
    ax[1, 1].legend(fontsize=8)

    # (1,2) LP±HP 的全通性
    hbp = hp_partner(half_order, cutoff)
    ax[1, 2].plot(w / np.pi, np.abs(zfreq(b + hbp, w) / zfreq(a, w)), "b", label="|LP + HP|")
    ax[1, 2].plot(w / np.pi, np.abs(zfreq(b - hbp, w) / zfreq(a, w)), "r--", label="|LP - HP|")
    ax[1, 2].axhline(1.0, color="k", lw=0.8, ls=":")
    ax[1, 2].set(title="LR crossover: one of LP$\\pm$HP is allpass, the other is not",
                 xlabel="$\\omega/\\pi$", ylabel="abs", ylim=(0, 2.2))
    ax[1, 2].grid(alpha=0.3)
    ax[1, 2].legend(fontsize=8)

    fig.suptitle(f"tf2ca even order with repeated poles: LR{2 * half_order} = (BW{half_order})^2, "
                 f"cutoff={cutoff:g}  ->  no two-chain form; needs three chains", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"图已保存：{out}")


# ------------------------------------------------------------
# main
# ------------------------------------------------------------


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="even-order Linkwitz-Riley (repeated poles) vs tf2ca")
    p.add_argument("--half-order", type=int, default=2,
                   help="BW 的半阶 n，滤波器为 LR_2n（默认 2 -> LR4）")
    p.add_argument("--cutoff", type=float, default=0.25, help="Nyquist 归一截止频率 (0..1)")
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--delays", action="store_true",
                   help="双链穷举时额外尝试链上的 z^-k 延迟偏移（j=0,1,2；明显更慢）")
    p.add_argument("--check", action="store_true", help="跑批量校验表")
    args = p.parse_args(argv)

    if args.half_order < 1:
        p.error("--half-order 必须 >= 1")

    if args.check:
        run_check()
        return

    shifts = (0, 1, 2) if args.delays else (0,)
    out = args.out or SCRIPT_DIR / "output" / f"lr_repeated_{2 * args.half_order}_{args.cutoff:g}.png"

    print(f"== Linkwitz–Riley LR{2 * args.half_order} = (BW{args.half_order} 阶)^2, "
          f"cutoff = {args.cutoff} ==")
    b, a = design_lr(args.half_order, args.cutoff)
    print(f"  B = {np.array2string(b, precision=6, suppress_small=True)}")
    print(f"  A = {np.array2string(a, precision=6, suppress_small=True)}")
    print()

    hint = lr_half_denominator(args.half_order, args.cutoff)
    fit = two_chain_fit(b, a, shifts, sqrt_hint=hint)
    chains = lr_three_chain(args.half_order, args.cutoff)
    cm = three_chain_metrics(b, a, chains)
    report(b, a, fit, chains, cm, args.half_order, args.cutoff, hint)
    if not args.delays:
        print("\n  提示：加 --delays 可额外尝试链上的 z^-k 延迟偏移（结论不变，但更慢）")
    plot(b, a, fit, chains, cm, args.half_order, args.cutoff, out)


if __name__ == "__main__":
    main()
