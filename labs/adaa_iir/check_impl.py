"""交叉验证 ``aaiir.py`` 的实现：三条互相独立的证据链。

1. **作者 MATLAB 参考代码逐行移植**（``AA_osc_cplx.m``）——验证公式转录与递推约定；
2. **卷积积分的高精度数值求积**（论文式 (6)：``y_n = ∫_0^n 2Re(B e^{βt}) f(x̃(n-t)) dt``）
   ——验证数学推导本身，不依赖任何索引技巧；
3. **方法间恒等式 / 连续卷积定义**：
   * ``DPW-2 == AA-FIR-1``、``DPW-3 == AA-FIR-2``（论文 §III-A）；
   * AA-FIR-1/2 与「矩形核 / 三角核连续卷积」逐点比对（验证核的宽度与归一化）；
   * 1 阶 Butterworth（实极点路径）与 AA 卷积求积比对（验证式 (20)–(25)）。

任一项超差即以非零退出码结束。

用法: python qwqdsp/labs/adaa_iir/check_impl.py
"""
from __future__ import annotations

import sys

import numpy as np
import scipy.integrate as si
import scipy.signal as ss

import aaiir as A

FS = 44100.0
TOL_REF = 1e-9      # 与参考代码逐点比（同一套公式，只差运算次序）
TOL_QUAD = 1e-8     # 与数值求积比
TOL_IDENT = 1e-12   # 恒等式（纯代数）
TOL_DPW3 = 1e-10    # DPW-3 与 AA-FIR-2 相差一个「线性于 n」项，抵消量级更大

failures = []


def report(name: str, err: float, tol: float) -> None:
    ok = err <= tol
    print(f"[{'PASS' if ok else 'FAIL'}] {name:<58} max|diff| = {err:.3e} (tol {tol:.0e})")
    if not ok:
        failures.append(name)


def report_pred(name: str, ok: bool, detail: str = "") -> None:
    print(f"[{'PASS' if ok else 'FAIL'}] {name:<58} {detail}")
    if not ok:
        failures.append(name)


# ------------------------------------------------------------
# 1. 作者参考代码 AA_osc_cplx.m 的逐行移植（1-based MATLAB 索引 → 0-based）
# ------------------------------------------------------------

def ref_aa_osc_cplx(x, B, beta, tab: A.PwlTable):
    """``shareable/AA_osc_cplx.m`` 的等价转录。

    两处 ``binary_search_*`` 用 ``np.searchsorted`` 等价替换（已按 x0 在左端 /
    右端 / 内部三种情形核对过），其余逐行照抄。
    """
    X, m, q = tab.X, tab.m, tab.q
    T, k = tab.T, tab.k
    m_diff = np.zeros(k)
    q_diff = np.zeros(k)
    m_diff[:k - 1] = m[1:] - m[:-1]
    q_diff[:k - 1] = q[1:] - q[:-1]
    m_diff[k - 1] = m[0] - m[k - 1]
    q_diff[k - 1] = q[0] - q[k - 1] - m[0] * T

    def mod_bar(v, kk):
        r = v % kk
        return r if r != 0 else kk

    def bsearch_down(x0):        # 最后一个 <= x0 的 1-based 下标，x0 < X(1) 时返回 0
        return int(np.searchsorted(X, x0, side="right"))

    def bsearch_up(x0):          # 第一个 >= x0 的 1-based 下标，x0 > X(end) 时返回 len+1
        return int(np.searchsorted(X, x0, side="left")) + 1

    expbeta = np.exp(beta)
    x_vz1, y_hat_vz1, x_diff_vz1 = 0.0, 0j, 0.0
    x_red = x_vz1 % T
    j_red = bsearch_down(x_red)
    j = k * int(np.floor(x_vz1 / T)) + j_red - 1
    y = np.zeros(len(x))
    for n in range(1, len(x)):
        x_diff = x[n] - x_vz1
        j_vz1 = j
        if (x_diff >= 0 and x_diff_vz1 >= 0) or (x_diff < 0 and x_diff_vz1 <= 0):
            j_vz1 = j + int(np.sign(x_red - X[j_red - 1]))
        x_red = x[n] % T
        if x_diff >= 0:
            j_red = bsearch_down(x_red)
            j = k * int(np.floor(x[n] / T)) + j_red - 1
            j_min, j_max = j_vz1, j
        else:
            j_red = bsearch_up(x_red)
            j = k * int(np.floor(x[n] / T)) + j_red - 1
            j_min, j_max = j, j_vz1
        j_min_bar = mod_bar(j_min, k)
        j_max_p_bar = mod_bar(j_max + 1, k)
        if x_diff >= 0:
            I = (expbeta * (m[j_min_bar - 1] * x_diff
                            + beta * (m[j_min_bar - 1] * (x_vz1 - T * np.floor((j_min - 1) / k))
                                      + q[j_min_bar - 1]))
                 - m[j_max_p_bar - 1] * x_diff
                 - beta * (m[j_max_p_bar - 1] * (x[n] - T * np.floor(j_max / k))
                           + q[j_max_p_bar - 1]))
        else:
            I = (expbeta * (m[j_max_p_bar - 1] * x_diff
                            + beta * (m[j_max_p_bar - 1] * (x_vz1 - T * np.floor(j_max / k))
                                      + q[j_max_p_bar - 1]))
                 - m[j_min_bar - 1] * x_diff
                 - beta * (m[j_min_bar - 1] * (x[n] - T * np.floor((j_min - 1) / k))
                           + q[j_min_bar - 1]))
        I_sum = 0j
        for l in range(j_min, j_max + 1):
            l_bar = mod_bar(l, k)
            I_sum += (np.exp(beta * (x[n] - X[l_bar] - T * np.floor((l - 1) / k)) / x_diff)
                      * (beta * q_diff[l_bar - 1]
                         + m_diff[l_bar - 1] * (x_diff + beta * X[l_bar])))
        I = (I + np.sign(x_diff) * I_sum) / beta ** 2
        y_hat = expbeta * y_hat_vz1 + 2.0 * B * I
        y[n] = y_hat.real
        x_vz1, y_hat_vz1, x_diff_vz1 = x[n], y_hat, x_diff
    return y


# ------------------------------------------------------------
# 2. 连续卷积（式 (6)）的高精度数值求积
# ------------------------------------------------------------

def ref_convolution(x, B, beta, tab: A.PwlTable, n_max: int):
    """``y_n = ∫_0^n 2Re(B·e^{βt})·f(x̃(n-t)) dt``，逐采样区间自适应求积。

    ``quad`` 不接受复被积函数，故实部/虚部分开积分后再合成。
    """
    idx = np.arange(len(x), dtype=float)

    def make(n):
        def integrand(t):
            xi = np.interp(n - t, idx, x)
            return np.exp(beta * t) * A.pwl_value(xi, tab)
        return integrand

    def cquad(f, a, b):
        return (si.quad(lambda t: f(t).real, a, b, epsabs=1e-13, epsrel=1e-13, limit=200)[0]
                + 1j * si.quad(lambda t: f(t).imag, a, b, epsabs=1e-13, epsrel=1e-13, limit=200)[0])

    y = np.zeros(n_max)
    for n in range(n_max):
        val = sum(cquad(make(n), s, s + 1) for s in range(n))
        y[n] = 2.0 * np.real(B * val)
    return y


def ref_kernel_conv(x, tab: A.PwlTable, n_max: int, kernel: str):
    """连续卷积：``y_n = ∫ h(t) f(x̃(n-t)) dt``。

    矩形核 ``h(t) = 1, t ∈ [0,1]``（论文式 (3)）；三角核 ``h(t) = 1-|t-1|, t ∈ [0,2]``
    （[Parker 2016] 的三角核，宽 2 采样）。起点 ``n < span`` 无定义（返回 0），
    比对时要跳过。
    """
    idx = np.arange(len(x), dtype=float)
    span = 1 if kernel == "rect" else 2
    y = np.zeros(n_max)
    for n in range(span, n_max):
        def integrand(t, n=n):
            xi = np.interp(n - t, idx, x)
            w = 1.0 if kernel == "rect" else 1.0 - abs(t - 1.0)
            return w * A.pwl_value(xi, tab)
        y[n] = sum(si.quad(integrand, s, s + 1, epsabs=1e-13, epsrel=1e-13, limit=200)[0]
                   for s in range(span))
    return y


# ------------------------------------------------------------
# 主流程
# ------------------------------------------------------------

def main() -> int:
    cases = [
        ("saw 1 kHz", A.saw_table(), 1000.0),
        ("saw 5 kHz", A.saw_table(), 5000.0),
        ("Escalation 1 kHz", A.escalation_ii_w3_table(), 1000.0),
        ("Escalation 7 kHz", A.escalation_ii_w3_table(), 7000.0),
    ]
    filters = {
        "AA-IIR-1": A.aa_iir_1(),
        "AA-IIR-2": A.aa_iir_2(),
        "butter-1": A.aa_filter_butter(1, A.AAIIR1_WN),
    }

    print("== 滤波器设计 ==")
    for name, f in filters.items():
        pairs, reals, direct = f
        direct_pole_pairs = len(pairs)
        print(f"  {name}: 共轭极点对 = {direct_pole_pairs}, 实极点 = {len(reals)}, "
              f"直接项 = {direct:.3e}")
    assert len(filters["AA-IIR-1"][0]) == 1
    assert len(filters["AA-IIR-2"][0]) == 5
    assert len(filters["butter-1"][1]) == 1
    # 2 阶 Butterworth 截止点 = 半功率点（-3.0103 dB）
    b, a = ss.zpk2tf(*ss.butter(2, A.AAIIR1_WN, btype="low", analog=True, output="zpk"))
    _, h = ss.freqs(b, a, np.array([A.AAIIR1_WN]))
    report("2 阶 Butterworth 截止点 = 半功率", abs(abs(h[0]) - 1 / np.sqrt(2.0)), 1e-12)
    # 10 阶 Chebyshev II：阻带内波纹峰值应恰为 -rs dB
    b2, a2 = ss.zpk2tf(*ss.cheby2(10, 60.0, A.AAIIR2_WN, btype="low", analog=True, output="zpk"))
    wstop = np.linspace(A.AAIIR2_WN, 200.0, 200001)
    _, h2 = ss.freqs(b2, a2, wstop)
    report("10 阶 Chebyshev II 阻带波纹峰值 = -60 dB",
           abs(20 * np.log10(np.max(np.abs(h2))) + 60.0), 1e-3)
    report("10 阶 Chebyshev II 直接项 = 阻带增益",
           abs(filters["AA-IIR-2"][2] - 10 ** (-60 / 20)), 1e-12)

    print("== 1) 向量化实现 vs 作者 MATLAB 参考代码 ==")
    for cname, tab, f0 in cases:
        x = A.phase_ramp(400, f0, FS)
        for fname in ("AA-IIR-1", "AA-IIR-2"):
            pairs, reals, _ = filters[fname]
            err = 0.0
            for B, beta in pairs:
                err = max(err, float(np.max(np.abs(
                    ref_aa_osc_cplx(x, B, beta, tab) - _single_pair(x, B, beta, tab)))))
            report(f"ref vs 本实现  {cname}, {fname}", err, TOL_REF)
        # 实极点路径（参考代码只处理共轭对，取 B/2 抵消它的 2B 因子）
        A_, alpha = filters["butter-1"][1][0]
        err = float(np.max(np.abs(ref_aa_osc_cplx(x, A_ / 2.0, alpha, tab)
                                  - A.aa_iir_osc(x, filters["butter-1"], tab))))
        report(f"ref(/2) vs 本实现  {cname}, butter-1（实极点）", err, TOL_REF)

    print("== 2) 卷积积分数值求积（式 (6)） ==")
    for cname, tab, f0 in cases:
        x = A.phase_ramp(24, f0, FS)
        for fname in ("AA-IIR-1", "butter-1"):
            pairs, reals, _ = filters[fname]
            if pairs:
                B, beta = pairs[0]
                y_ref = ref_convolution(x, B, beta, tab, len(x))
            else:
                A_, alpha = reals[0]
                y_ref = ref_convolution(x, A_ / 2.0, alpha + 0j, tab, len(x))
            y_got = A.aa_iir_osc(x, filters[fname], tab)
            report(f"quadrature vs 本实现  {cname}, {fname}",
                   float(np.max(np.abs(y_ref - y_got))), TOL_QUAD)

    print("== 3) 方法恒等式（论文 §III-A） ==")
    tab = A.saw_table()
    for f0 in (1000.0, 5000.0, 12000.0):
        n = 4096
        e = float(np.max(np.abs(A.dpw(tab, n, f0, FS, 2) - A.aa_fir(tab, n, f0, FS, 1))))
        report(f"DPW-2 == AA-FIR-1 @ {f0:.0f} Hz", e, TOL_IDENT)
        e = float(np.max(np.abs(A.dpw(tab, n, f0, FS, 3) - A.aa_fir(tab, n, f0, FS, 2))))
        report(f"DPW-3 == AA-FIR-2 @ {f0:.0f} Hz", e, TOL_DPW3)

    print("== 4) AA-FIR 核宽度/归一化（连续卷积定义） ==")
    tab = A.saw_table()
    x = A.phase_ramp(24, 3000.0, FS)
    for order, kernel in ((1, "rect"), (2, "tri")):
        y_ref = ref_kernel_conv(x, tab, len(x), kernel)
        y_got = A.aa_fir(tab, len(x), 3000.0, FS, order)
        report(f"连续 {kernel} 核卷积 vs AA-FIR-{order}",
               float(np.max(np.abs(y_ref[order:] - y_got[order:]))), 1e-7)

    print("== 5) AA-IIR 通带增益不变量（滤波器归一化 / 2B 因子） ==")
    # (a) 常数输入：稳态输出 = 输入 × (1 - 被丢弃的直接项)
    #     参考实现忽略了 residue 的直接项 A0（Cheby2 的 A0 = 1e-3 = 阻带增益），
    #     因此 AA-IIR-2 的直流增益是 1 - 1e-3 而非 1（见 README「与参考实现的差异」）。
    const = A.PwlTable(X=[0.0, 1.0], m=[0.0], q=[1.0])
    for fname in ("AA-IIR-1", "AA-IIR-2"):
        pairs, reals, direct = filters[fname]
        y = A.aa_iir_osc(A.phase_ramp(4096, 100.0, FS), filters[fname], const)
        report(f"DC 增益 = 1 - A0  {fname}（A0 = {direct:.1e}）",
               abs(float(y[1024:].mean()) - (1.0 - direct)), 1e-9)
        yi = A.aa_iir_osc(A.phase_ramp(4096, 100.0, FS), filters[fname], const,
                          include_direct=True)
        report(f"含直接项时 DC 增益 = 1  {fname}",
               abs(float(yi[1024:].mean()) - 1.0), 1e-9)
    # (b) 正弦表（256 段线性插值）@ 441 Hz：基波幅度不应被改动（同上，差 A0）
    k = 256
    Xs = np.arange(k + 1) / k
    ys = np.sin(2 * np.pi * Xs)
    m = (ys[1:] - ys[:-1]) * k
    q = ys[:-1] - m * Xs[:-1]
    sine = A.PwlTable(X=Xs, m=m, q=q, wt=np.sin(2 * np.pi * np.arange(2048) / 2048))
    ns = 8192
    f0s = 441.0
    w = ss.windows.blackmanharris(ns)
    ref = np.abs(np.fft.rfft(np.sin(2 * np.pi * f0s / FS * np.arange(ns)) * w)).max()
    for fname in ("AA-IIR-1", "AA-IIR-2"):
        direct = filters[fname][2]
        y = A.aa_iir_osc(A.phase_ramp(ns, f0s, FS), filters[fname], sine)
        amp = np.abs(np.fft.rfft(y * w)).max()
        report(f"441 Hz 基波增益 = 1 - A0  {fname}",
               abs(20 * np.log10(amp / ref) - 20 * np.log10(1.0 - direct)), 0.002)

    print("== 6) SNR 口径的解析锚点（trivial 锯齿波，理想折叠分量功率比） ==")
    tab = A.saw_table()
    for f0 in (440.0, 1000.0, 2000.0, 3000.0, 5000.0):
        y = A.trivial(tab, 2 ** 16, f0, FS)
        got = A.snr_db(y, FS, f0, guard_bins=4, warmup=4096)
        ref = A.analytic_saw_snr(f0, FS)
        report(f"trivial SNR 实测 vs 解析 @ {f0:.0f} Hz", abs(got - ref), 1.0)
    # 病态音（折叠分量落在谐波栅格上）必须被门限识别出来
    bad = A.trivial(A.saw_table(), 2 ** 16, 49.0, FS)
    got = A.snr_db(bad, FS, 49.0, guard_bins=4, warmup=4096)
    dev = got - A.analytic_saw_snr(49.0, FS)
    report_pred("病态音 49 Hz 被门限识别（偏差 > 2 dB）", dev > 2.0,
                f"偏差 {dev:+.1f} dB")

    print()
    if failures:
        print(f"FAILED: {len(failures)} 项 -> {failures}")
        return 1
    print("全部通过。")
    return 0


def _single_pair(x, B, beta, tab):
    """只取单个共轭极点对的 AA-IIR 递推（与 ``aa_iir_osc`` 内部同一路径）。"""
    return A.aa_iir_osc(x, ([(B, beta)], [], 0.0), tab)


if __name__ == "__main__":
    sys.exit(main())
