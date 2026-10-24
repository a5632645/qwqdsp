# -*- coding: utf-8 -*-
"""
even_order.py
=============

**偶数阶**实系数低通 ``H(z) = B(z)/A(z)`` 能不能写成「两条**复系数**全通链 + 复权重」之和？

    H(z) = alpha * A1(z) + beta * A2(z)          alpha, beta ∈ C
    A1   = conj(反序(D1)) / D1                   （复全通：|A1(e^jw)| ≡ 1）
    D2   = conj(D1)（系数逐个取共轭）

本脚本给出的结论是：

1. **可以**。四族（butter / cheby1 / cheby2 / ellip）偶数阶 2..8、多种截止频率，
   ``B`` 的还原残差、逐极点 ``alpha`` 离散度、频响误差通常都在 1e-12 以下；
   窄带高阶例外，原因见下（``--check`` 有完整表格）。
2. ``D2 = conj(D1)`` 时 ``A2`` 是 ``A1`` 的**系数共轭**（不是逐点共轭）：在单位圆上
   ``A2(e^jw) = conj(A1(e^-jw))``。于是 ``B = alpha*N1*D2 + beta*N2*D1`` 里两项互为系数共轭，
   实数 ``B`` 强制

       beta = conj(alpha)          =>   对实输入  y = 2 Re{ alpha * (A1 x) }

   （脚本用**不假设共轭**的 4 参数最小二乘独立求解 alpha、beta，再报告
   ``|beta - conj(alpha)|``；并用 ``2 Re{alpha*A1*x}`` 与直接型冲激响应对比，
   确认「只算一条复链」这个说法成立。）
3. 反过来说，**实系数**两路全通做不到偶数阶：``n = dA + dB`` 为偶数 ⇒ ``dA、dB`` 同奇偶，
   于是 ``H(-1) = ±(gA + gB)/2 = ±1``（由 ``H(1) = 1`` 得 ``gA + gB = 2``），
   而低通要求 ``|H(-1)| < 1``。放开复系数后复全通在 ``z = -1`` 处可以是模长 1 的任意复数，
   这个约束才消失。脚本对同一批设计同时跑「实系数最佳残差」作对照。

怎么判断「分对了」：``B(p) = alpha*N1(p)*D2(p)``（``p ∈ D1``）与
``B(p) = conj(alpha)*N2(p)*D1(p)``（``p ∈ D2``）逐极点代入即可独立解出 ``alpha``，
正确分组下这些 ``alpha`` 必须相等——这是 :func:`pole_alpha_spread`，只做点值代入、
不展开多项式，条件数远好于系数空间最小二乘，是本目录判定分组的主要依据。
窄带 / 高阶（如 ``cheby1`` 8 阶 ``cutoff=0.1``）时系数表示的消位会把残差抬到 1e-3 量级，
但**同样的**病态在奇数阶「已知精确解」上同样出现（``--check`` 第二张表），属输入数据
（双精度系数 -> 极点）的条件数限制，不是分组错了。

推导与数据见同目录 ``even_order.md``；奇数阶（实系数、``1/2 (A0 + A1)``）见 ``odd_order.py``。

用法
----
    python even_order.py                            # 默认 ellip 6 阶, cutoff=0.3
    python even_order.py --kind butter --order 4 --cutoff 0.25
    python even_order.py --exact                    # 额外用 60 位精度重算同一分组（看残差来自数据还是算术）
    python even_order.py --check                    # 四族 × 偶数阶 2/4/6/8 × 多截止频率
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 无界面后端
import matplotlib.pyplot as plt
import numpy as np

from tf2ca import design, split_by_radius, tf_response, zfreq

SCRIPT_DIR = Path(__file__).resolve().parent


# ------------------------------------------------------------
# 复全通与极点分组
# ------------------------------------------------------------


def conj_rev(d: np.ndarray) -> np.ndarray:
    """复全通的分子 ``D^# = conj(反序(D))``，使 ``D^#/D`` 在单位圆上模长恒为 1。"""
    return np.conj(np.asarray(d)[::-1])


def conjugate_splits(a: np.ndarray, tol: float = 1e-7):
    """枚举「每个共轭对取一个根给 D1、另一半由 ``D2 = conj(D1)`` 承担」的全部分法。

    共轭对按 ``|p|`` 升序编号，第 ``k`` 位为 1 表示取该对的**上半平面根**。
    实极点必须成对（重数为偶、数值相等）各占一条链；否则 ``D2 = conj(D1)`` 不可能成立，
    直接抛 ``ValueError``（四族偶数阶设计无实极点，不受此限）。

    :return: ``(candidates, n_pairs)``，``candidates`` 为 ``[(D1, bits), ...]``，共 ``2^n_pairs`` 项。
    """
    poles = np.roots(np.asarray(a, dtype=float))
    real = sorted(float(p.real) for p in poles if abs(p.imag) < tol)
    if len(real) % 2:
        raise ValueError(f"实极点共 {len(real)} 个，无法配成共轭共享对：D2 = conj(D1) 不成立")
    fixed: list[float] = []
    for k in range(0, len(real), 2):
        if abs(real[k] - real[k + 1]) > tol:
            raise ValueError(f"实极点 {real[k]:.6g} 与 {real[k + 1]:.6g} 不相等，无法配对")
        fixed.append(real[k])  # 一实根给 D1，另一份由 D2 = conj(D1) 承担

    upper = sorted((p for p in poles if p.imag >= tol), key=lambda p: abs(p))
    out: list[tuple[np.ndarray, int]] = []
    for bits in range(1 << len(upper)):
        chosen = [p if (bits >> k) & 1 else p.conjugate() for k, p in enumerate(upper)]
        out.append((np.poly(chosen + fixed), bits))
    return out, len(upper)


def complex_allpass_response(d: np.ndarray, w: np.ndarray) -> np.ndarray:
    """复全通 ``D^#/D`` 在 ``z = exp(jw)`` 上的复响应。"""
    return zfreq(conj_rev(d), w) / zfreq(d, w)


# ------------------------------------------------------------
# 复权重：最小二乘
# ------------------------------------------------------------


def weights_conjugate(b: np.ndarray, d1: np.ndarray):
    """在 ``D2 = conj(D1)``、``beta = conj(alpha)`` 下解 ``B = alpha*P + conj(alpha)*conj(P)``。

    ``P = N1*D2``（``N1 = conj(反序(D1))``），即 ``alpha`` 是**第一条链** ``A1 = N1/D1`` 的权重
    （``alpha*P`` 正是 ``H`` 通分后的第一个分子项 ``alpha*N1*D2``）。``U = P + conj(P)``、
    ``V = j(P - conj(P))`` 都是实向量，解 ``[U V] @ [x, y]^T = B`` 得 ``alpha = x + j*y``。

    :return: ``(alpha, residual, imag_leak)``，残差为 ``max|B_还原 - B| / max|B|``，
        ``imag_leak`` 为还原分子虚部的相对泄漏（应为机器精度，说明 B 自动是实的）。
    """
    b = np.asarray(b, dtype=float)
    d2 = np.conj(d1)
    p = np.convolve(conj_rev(d1), d2)
    u, v = p + np.conj(p), 1j * (p - np.conj(p))
    scale = float(np.max(np.abs(b)))
    m = np.stack([u.real, v.real], axis=1)
    xy, *_ = np.linalg.lstsq(m, b, rcond=None)
    alpha = xy[0] + 1j * xy[1]
    recon = alpha * np.convolve(conj_rev(d1), d2) + np.conj(alpha) * np.convolve(conj_rev(d2), d1)
    return alpha, float(np.max(np.abs(recon - b)) / scale), float(np.max(np.abs(recon.imag)) / scale)


def weights_free(b: np.ndarray, d1: np.ndarray):
    """**不假设** ``beta = conj(alpha)``：4 个实未知量独立最小二乘解 ``alpha*P + beta*Q``。

    :return: ``(alpha, beta, residual)``，用来验证共轭权重是自动出现的而非强加。
    """
    b = np.asarray(b, dtype=float)
    d2 = np.conj(d1)
    p = np.convolve(conj_rev(d1), d2)
    q = np.convolve(conj_rev(d2), d1)
    scale = float(np.max(np.abs(b)))
    m = np.stack(
        [np.concatenate([c.real, c.imag]) for c in (p, 1j * p, q, 1j * q)], axis=1
    )
    rhs = np.concatenate([b, np.zeros_like(b)])
    sol, *_ = np.linalg.lstsq(m, rhs, rcond=None)
    alpha, beta = sol[0] + 1j * sol[1], sol[2] + 1j * sol[3]
    recon = alpha * p + beta * q
    return alpha, beta, float(np.max(np.abs(recon - b)) / scale)


# ------------------------------------------------------------
# 分解与校验
# ------------------------------------------------------------


def decompose(b: np.ndarray, a: np.ndarray) -> dict:
    """穷举 ``2^(n/2)`` 种共轭对分法，取 ``B`` 还原残差最小的复全通解。

    :return: 含 ``d1``、``d2``、``alpha``、``beta``、``residual``、``bits``、
        ``n_pairs``、``table``（每种分法的残差）的字典。
    """
    cands, n_pairs = conjugate_splits(a)
    table: list[tuple[int, float]] = []
    best: dict | None = None
    for d1, bits in cands:
        alpha, res, leak = weights_conjugate(b, d1)
        table.append((bits, res))
        if best is None or res < best["residual"]:
            best = {
                "d1": d1,
                "d2": np.conj(d1),
                "alpha": alpha,
                "beta": np.conj(alpha),
                "residual": res,
                "imag_leak": leak,
                "bits": bits,
            }
    assert best is not None
    best["n_pairs"] = n_pairs
    best["table"] = table
    return best


def metrics(b: np.ndarray, a: np.ndarray, decomp: dict, num: int = 4000) -> dict:
    """数值校验：频响误差、两链全通平坦度、共轭对称性、``A`` 还原误差、``|H(-1)|``。"""
    b = np.asarray(b, dtype=float)
    a = np.asarray(a, dtype=float)
    d1, d2 = decomp["d1"], decomp["d2"]
    alpha, beta = decomp["alpha"], decomp["beta"]
    w = np.linspace(1e-6, np.pi - 1e-6, num)

    a1 = complex_allpass_response(d1, w)
    a2 = complex_allpass_response(d2, w)
    # 「A2 = conj(A1)」是**系数**共轭：即 A2(e^jw) = conj(A1(e^-jw))，不是逐点 conj。
    a1_mirror = complex_allpass_response(d1, -w)
    got = alpha * a1 + beta * a2
    ref = tf_response(b, a, w)

    a_hat = np.convolve(d1, d2)  # 应与 A 相等（复共轭分母对相乘 -> 实）
    a_scale = float(np.max(np.abs(a)))
    return {
        "freq_err": float(np.max(np.abs(got - ref))),
        "flat": float(max(np.max(np.abs(np.abs(a1) - 1.0)), np.max(np.abs(np.abs(a2) - 1.0)))),
        "conj_err": float(np.max(np.abs(a2 - np.conj(a1_mirror)))),
        "a_err": float(np.max(np.abs(a_hat - a))) / a_scale,
        "h_nyquist": float(abs(tf_response(b, a, np.array([np.pi]))[0])),
        "freq": w,
        "ref": ref,
        "got": got,
        "a1": a1,
        "a2": a2,
        "a1_mirror": a1_mirror,
    }


def pole_alpha_spread(b: np.ndarray, d1: np.ndarray, d2: np.ndarray | None = None):
    """在 ``A`` 的每个极点处**独立**解出 ``alpha``，返回这些解的相对离散度。

    对 ``p ∈ D1``：``alpha = B(p) / (N1(p)*D2(p))``；对 ``p ∈ D2``：
    ``alpha = conj( B(p) / (N2(p)*D1(p)) )``。分组正确时这些值必须相等，所以离散度
    是比系数空间最小二乘**条件数好得多**的判据（只做点值代入，不展开多项式、不整体拟合）——
    本目录主要靠它判断「分对了没有」（错的分法会离散到 O(1)）。

    :return: ``(spread, values)``，``spread = max|a_i - mean| / |mean|``。
    """
    d2 = np.conj(d1) if d2 is None else d2
    n1, n2 = conj_rev(d1), conj_rev(d2)
    vals = [np.polyval(b, p) / (np.polyval(n1, p) * np.polyval(d2, p)) for p in np.roots(d1)]
    vals += [np.conj(np.polyval(b, p) / (np.polyval(n2, p) * np.polyval(d1, p)))
             for p in np.roots(d2)]
    vals = np.asarray(vals)
    mean = vals.mean()
    return float(np.max(np.abs(vals - mean)) / abs(mean)), vals


def single_chain_form(b: np.ndarray, a: np.ndarray, d1: np.ndarray, alpha: complex,
                      num: int = 200) -> float:
    """验证「一条复全通 + 复增益 + 取实部」：``y = 2 Re{ alpha * (A1 x) }``。

    系数共轭的 ``A2`` 对**实输入** ``x`` 等价于对滤波输出取共轭
    （``conj(A1*x) = A2*x``），所以 ``alpha*A1 + conj(alpha)*A2`` 可以只算一条复链再取实部。
    返回该结构与直接型冲激响应的最大偏差。
    """
    from scipy.signal import lfilter

    x = np.zeros(num)
    x[0] = 1.0
    return float(np.max(np.abs(lfilter(b, a, x)
                                     - 2.0 * np.real(alpha * lfilter(conj_rev(d1), d1, x)))))


def exact_residual(b: np.ndarray, d1: np.ndarray, digits: int = 60):
    """用 ``mpmath`` 在 ``digits`` 位精度重算**同一分组**的 ``B`` 还原残差（QR 最小二乘）。

    与双精度结果对比即可判断残差来自哪里：80 位下残差**不变**说明误差在**输入数据**里
    （被舍入到双精度的 ``A``、``B`` 系数），与浮点算术无关；也就是说这些设计的
    ``B = alpha*N1*D2 + beta*N2*D1`` 恒等式条件数太大。见 ``even_order.md``。

    :return: ``(residual, alpha, beta)``，``alpha/beta`` 为高精度解（``mpmath`` 复数）。
    """
    import mpmath as mp

    d1 = np.asarray(d1)
    d2 = np.conj(d1)

    def mp_poly(roots):
        poly = [mp.mpc(1)]
        for r in roots:
            new = [mp.mpc(0)] * (len(poly) + 1)
            for k, c in enumerate(poly):
                new[k] += c
                new[k + 1] += -c * r
            poly = new
        return poly

    def mp_conv(x, y):
        out = [mp.mpc(0)] * (len(x) + len(y) - 1)
        for i, xi in enumerate(x):
            for j, yj in enumerate(y):
                out[i + j] += xi * yj
        return out

    def to_mp(p):
        return [mp.conj(c) for c in p[::-1]]

    with mp.workdps(digits):
        mkr = lambda p: mp.mpc(float(p.real), float(p.imag))  # noqa: E731
        dm = [mp_poly([mkr(p) for p in np.roots(d)]) for d in (d1, d2)]
        P = mp_conv(to_mp(dm[0]), dm[1])
        Q = mp_conv(to_mp(dm[1]), dm[0])
        rhs = mp.matrix([mp.mpf(float(v)) for v in b])
        A = mp.matrix(len(P), 2)
        for k in range(len(P)):
            A[k, 0], A[k, 1] = P[k], Q[k]
        (sol, _) = mp.qr_solve(A, rhs)  # 复 2 未知量最小二乘；实嵌入矩阵在共轭分组下秩亏
        y = [sol[0] * P[k] + sol[1] * Q[k] for k in range(len(P))]
        scale = max(abs(v) for v in rhs)
        res = float(max(abs(y[k] - rhs[k]) for k in range(len(P))) / scale)
        return res, complex(sol[0]), complex(sol[1])


def real_two_path_min_residual(b: np.ndarray, a: np.ndarray) -> float:
    """实系数两路全通 ``1/2 [gA*AA + gB*AB]`` 的最佳 ``B`` 残差（穷举所有分法）。

    共轭对必须**整对**归同一条链（否则该链系数不为实），实极点各自可归任一条链。
    用来做对照：奇数阶应达机器精度，偶数阶做不到（见 ``even_order.md``）。
    """
    b = np.asarray(b, dtype=float)
    a = np.asarray(a, dtype=float)
    poles = np.roots(a)
    real = sorted(float(p.real) for p in poles if abs(p.imag) < 1e-7)
    upper = sorted((p for p in poles if p.imag >= 1e-7), key=lambda p: abs(p))

    scale = float(np.max(np.abs(b)))
    best = np.inf
    for bits in range(1 << len(upper)):
        for rbits in range(1 << len(real)):
            da, db = np.array([1.0]), np.array([1.0])
            for k, p in enumerate(upper):
                pair = np.poly([p, p.conjugate()]).real  # 共轭对整对归一条链 -> 系数为实
                if (bits >> k) & 1:
                    da = np.convolve(da, pair)
                else:
                    db = np.convolve(db, pair)
            for k, r in enumerate(real):
                if (rbits >> k) & 1:
                    da = np.convolve(da, [1.0, -r])
                else:
                    db = np.convolve(db, [1.0, -r])
            na, nb = da[::-1], db[::-1]
            m = np.stack([np.convolve(na, db) / 2.0, np.convolve(nb, da) / 2.0], axis=1)
            g, *_ = np.linalg.lstsq(m, b, rcond=None)
            res = float(np.max(np.abs(m @ g - b))) / scale
            best = min(best, res)
    return best


def describe_chain(name: str, d: np.ndarray) -> None:
    """打印一条复全通链的阶数与分母系数。"""
    print(f"  {name}(z): 阶数 = {d.size - 1}，复系数分母 = "
          f"{np.array2string(d, precision=6, suppress_small=True)}")


def report(b: np.ndarray, a: np.ndarray, decomp: dict, m: dict, exact: bool = False) -> None:
    """打印分组、权重与全部校验量。"""
    print(f"  共轭对 {decomp['n_pairs']} 个 -> 候选分法 {1 << decomp['n_pairs']} 种，"
          f"最优 bits = {decomp['bits']:0{decomp['n_pairs']}b}（按 |p| 升序，1 = 取上半平面根）")
    describe_chain("D1", decomp["d1"])
    describe_chain("D2", decomp["d2"])
    print(f"  alpha = {decomp['alpha']:.9g}   beta = {decomp['beta']:.9g}")
    print(f"  频响最大偏差 |H_struct - H_direct| = {m['freq_err']:.3e}")
    print(f"  B 还原残差（相对）              = {decomp['residual']:.3e}")
    print(f"  全通链平坦度 max||A_i| - 1|     = {m['flat']:.3e}")
    print(f"  共轭对称性 max|A2 - conj(A1)|   = {m['conj_err']:.3e}")
    print(f"  A 还原误差 max|D1*D2 - A|/|A|   = {m['a_err']:.3e}")

    alpha, beta, res_free = weights_free(b, decomp["d1"])
    scale = max(abs(alpha), abs(beta))
    print(f"  独立求解 alpha/beta（不设共轭）：残差 {res_free:.3e}，"
          f"|beta - conj(alpha)|/max|.| = {abs(beta - np.conj(alpha)) / scale:.3e}")
    spread, _ = pole_alpha_spread(b, decomp["d1"], decomp["d2"])
    print(f"  逐极点独立解 alpha 的离散度      = {spread:.3e}")
    print(f"  单链形式 max|2Re(alpha*A1*x) - y| = "
          f"{single_chain_form(b, a, decomp['d1'], decomp['alpha']):.3e}")
    print(f"  |H(-1)|：设计 = {m['h_nyquist']:.3e}，"
          f"实系数两路全通此时恒为 1（故偶数阶必须复系数）")

    if exact:
        res_mp, al_mp, be_mp = exact_residual(b, decomp["d1"])
        print(f"  60 位精度重算同一分组：B 还原残差 = {res_mp:.3e}"
              f"（双精度 {decomp['residual']:.3e}），"
              f"|beta - conj(alpha)| = {abs(be_mp - np.conj(al_mp)):.1e}")


def run_check() -> None:
    """批量校验：四族 × 偶数阶 × 多截止频率，并与实系数最佳残差对照。"""
    print("== 偶数阶复全通分解（--check）==")
    print(f"{'family':8s}{'n':>3s}{'cutoff':>8s}{'B残差':>11s}{'极点离散度':>12s}"
          f"{'频响误差':>11s}{'|A_i|-1':>10s}{'A还原':>10s}{'实系数最佳':>12s}{'分法':>7s}")
    worst = 0.0
    for kind in ("butter", "cheby1", "cheby2", "ellip"):
        for order in (2, 4, 6, 8):
            for cutoff in (0.1, 0.2, 0.35, 0.5):
                b, a = design(kind, order, cutoff)
                dc = decompose(b, a)
                m = metrics(b, a, dc, num=2000)
                spread, _ = pole_alpha_spread(b, dc["d1"], dc["d2"])
                real_res = real_two_path_min_residual(b, a)
                worst = max(worst, m["freq_err"])
                bits = f"{dc['bits']:0{dc['n_pairs']}b}"
                print(f"{kind:8s}{order:3d}{cutoff:8.2f}{dc['residual']:11.3e}{spread:12.3e}"
                      f"{m['freq_err']:11.3e}{m['flat']:10.3e}{m['a_err']:10.3e}"
                      f"{real_res:12.3e}{bits:>7s}")
    print(f"  最差频响偏差 = {worst:.3e}")
    print()
    print("== 对照：奇数阶（实系数两路全通即可，见 odd_order.py）==")
    print(f"{'family':8s}{'n':>3s}{'cutoff':>8s}{'实系数最佳':>13s}{'逐极点alpha离散度':>18s}{'复共轭式':>12s}")
    for kind in ("butter", "cheby1", "cheby2", "ellip"):
        for order in (3, 5, 7):
            for cutoff in (0.2, 0.35):
                b, a = design(kind, order, cutoff)
                d0, d1 = split_by_radius(a)  # 已知精确的实分组：alpha = beta = 1/2
                spread, _ = pole_alpha_spread(b, d0, d1)
                try:
                    dc = decompose(b, a)
                    cl = f"{dc['residual']:12.3e}"
                except ValueError:
                    cl = "不适用(实根奇重)"
                print(f"{kind:8s}{order:3d}{cutoff:8.2f}"
                      f"{real_two_path_min_residual(b, a):13.3e}{spread:18.3e}{cl:>12s}")


# ------------------------------------------------------------
# 绘图
# ------------------------------------------------------------


def plot(b: np.ndarray, a: np.ndarray, decomp: dict, m: dict, out: Path) -> None:
    """六张图：幅频对比、z 平面分组、全通平坦度、分法残差、冲激响应、频响误差。"""
    w = m["freq"]
    d1, d2 = decomp["d1"], decomp["d2"]
    fig, ax = plt.subplots(2, 3, figsize=(16, 8.5))

    # (0,0) 幅频：直接型 vs 复全通结构
    ax[0, 0].plot(w / np.pi, 20 * np.log10(np.abs(m["ref"]) + 1e-18), "k", lw=2, label="direct B/A")
    ax[0, 0].plot(w / np.pi, 20 * np.log10(np.abs(m["got"]) + 1e-18), "r--", lw=1,
                  label="alpha*A1 + beta*A2")
    ax[0, 0].set(title="mag: direct-form vs complex two-path", xlabel="$\\omega/\\pi$",
                 ylabel="dB", ylim=(-80, 3))
    ax[0, 0].grid(alpha=0.3)
    ax[0, 0].legend()

    # (0,1) z 平面：共轭对拆成两条复链
    theta = np.linspace(0, 2 * np.pi, 400)
    ax[0, 1].plot(np.cos(theta), np.sin(theta), "k:", lw=0.8)
    p1, p2 = np.roots(d1), np.roots(d2)
    ax[0, 1].plot(p1.real, p1.imag, "x", color="#1f77b4", ms=9, label="D1 poles")
    ax[0, 1].plot(p2.real, p2.imag, "x", color="#d62728", ms=9, label="D2 = conj(D1) poles")
    for p in p1:  # 每个极点与它的共轭（在另一条链）连线
        ax[0, 1].plot([p.real, p.real], [p.imag, -p.imag], color="gray", lw=0.6, alpha=0.7)
    z = np.roots(b)
    ax[0, 1].plot(z.real, z.imag, "o", mfc="none", mec="g", ms=7, label="B zeros")
    ax[0, 1].axhline(0, color="gray", lw=0.5)
    ax[0, 1].axvline(0, color="gray", lw=0.5)
    ax[0, 1].set(title="z-plane: conjugate pairs split across chains", xlabel="Re", ylabel="Im",
                 aspect="equal")
    ax[0, 1].grid(alpha=0.3)
    ax[0, 1].legend(loc="upper left", fontsize=8)

    # (0,2) 两条链的全通平坦度
    ax[0, 2].semilogy(w / np.pi, np.abs(np.abs(m["a1"]) - 1) + 1e-18, label="$|A_1|-1$")
    ax[0, 2].semilogy(w / np.pi, np.abs(np.abs(m["a2"]) - 1) + 1e-18, "--", label="$|A_2|-1$")
    ax[0, 2].set(title="allpass flatness (both chains)", xlabel="$\\omega/\\pi$", ylabel="|A|-1")
    ax[0, 2].grid(alpha=0.3, which="both")
    ax[0, 2].legend(fontsize=8)

    # (1,0) 各分法的残差：只有极少数分法成立
    bits = [t[0] for t in decomp["table"]]
    res = [t[1] for t in decomp["table"]]
    colors = ["#d62728" if b_ == decomp["bits"] else "#1f77b4" for b_ in bits]
    ax[1, 0].bar(range(len(bits)), [max(r, 1e-18) for r in res], color=colors)
    ax[1, 0].set_yscale("log")
    ax[1, 0].set(title="B residual of every $2^{n/2}$ conjugate split (red = best)",
                 xlabel="split index", ylabel="rel. residual")
    ax[1, 0].grid(alpha=0.3, axis="y")

    # (1,1) 冲激响应（结构输出应为实）
    from scipy.signal import lfilter

    n_imp = 200
    x = np.zeros(n_imp)
    x[0] = 1.0
    y_direct = lfilter(b, a, x)
    y_struct = decomp["alpha"] * lfilter(conj_rev(d1), d1, x) + decomp["beta"] * lfilter(conj_rev(d2), d2, x)
    y_single = 2.0 * np.real(decomp["alpha"] * lfilter(conj_rev(d1), d1, x))
    ax[1, 1].plot(y_direct, "k", lw=2, label="direct")
    ax[1, 1].plot(y_struct.real, "r--", label="alpha*A1 + beta*A2")
    ax[1, 1].plot(y_single, "b:", label="2 Re(alpha*A1)  [single chain]")
    ax[1, 1].set(title=f"impulse response (max|$\\Delta$|={np.max(np.abs(y_direct - y_struct)):.1e})",
                 xlabel="n")
    ax[1, 1].grid(alpha=0.3)
    ax[1, 1].legend(fontsize=8)

    # (1,2) 频响误差与系数共轭关系
    ax[1, 2].semilogy(w / np.pi, np.abs(m["got"] - m["ref"]) + 1e-18, label="$|\\Delta H|$")
    ax[1, 2].semilogy(w / np.pi, np.abs(m["a2"] - np.conj(m["a1_mirror"])) + 1e-18, "--",
                      label="$|A_2-\\overline{A_1}|$ (coeff conj)")
    ax[1, 2].set(title="reconstruction error", xlabel="$\\omega/\\pi$", ylabel="abs")
    ax[1, 2].grid(alpha=0.3, which="both")
    ax[1, 2].legend(fontsize=8)

    fig.suptitle("tf2ca even order: H = alpha*A1 + beta*A2  (complex allpass, beta = conj(alpha))",
                 fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"图已保存：{out}")


# ------------------------------------------------------------
# main
# ------------------------------------------------------------


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="even-order IIR lowpass -> complex two-path allpass")
    p.add_argument("--kind", default="ellip", choices=["butter", "cheby1", "cheby2", "ellip"])
    p.add_argument("--order", type=int, default=6, help="偶数阶")
    p.add_argument("--cutoff", type=float, default=0.3, help="Nyquist 归一截止频率 (0..1)")
    p.add_argument("--ripple", type=float, default=0.2, help="通带波纹 dB (cheby1/ellip)")
    p.add_argument("--stop", type=float, default=40.0, help="阻带衰减 dB (cheby2/ellip)")
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--check", action="store_true", help="跑批量校验表")
    p.add_argument("--exact", action="store_true",
                   help="用 mpmath 在 60 位精度重算同一分组的 B 还原残差（判断残差来自数据还是算术）")
    args = p.parse_args(argv)

    if args.order % 2 != 0:
        p.error("本脚本只处理偶数阶：--order 必须为偶数（奇数阶见 odd_order.py）")

    if args.check:
        run_check()
        return

    out = args.out or SCRIPT_DIR / "output" / f"even_order_{args.kind}_{args.order}_{args.cutoff:g}.png"

    print(f"== 设计 {args.kind} 低通（偶数阶）：order={args.order}, cutoff={args.cutoff} ==")
    b, a = design(args.kind, args.order, args.cutoff, args.ripple, args.stop)
    print(f"  B = {np.array2string(b, precision=6, suppress_small=True)}")
    print(f"  A = {np.array2string(a, precision=6, suppress_small=True)}")
    print()

    print("== 复全通分解：H = alpha*A1 + beta*A2，D2 = conj(D1) ==")
    dc = decompose(b, a)
    m = metrics(b, a, dc)
    report(b, a, dc, m, args.exact)
    print()

    print(f"  实系数对照：两路全通最佳 B 残差 = {real_two_path_min_residual(b, a):.3e}"
          f"（偶数阶做不到，见 even_order.md）")
    plot(b, a, dc, m, out)


if __name__ == "__main__":
    main()
