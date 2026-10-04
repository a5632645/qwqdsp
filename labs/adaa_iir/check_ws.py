"""交叉验证 `ws.py`（ADAA-IIR waveshaper）。

四条独立证据链 + 两条鲁棒性检查：

1. **作者参考代码 `aaiir_hard_complex.m` 逐行移植** —— 硬削波 + 复共轭极点对；
2. **连续卷积的数值求积**（式 (5)：``y_n = ∫_0^n h(t)·f(x̃(n-t))dt``），不依赖任何索引技巧；
3. **原函数 vs 解析积分**（`scipy.integrate.quad`）；
4. **附录 B 的线性恒等式**：f(x)=x 时 AA-IIR 必须**严格等于**一阶低通的冲激不变滤波器
   ``H(z) = -1/α · [ (e^α-α-1) + ((α-1)e^α+1)z^{-1} ] / (1 - e^α z^{-1})``（论文原文给出）；
5. 鲁棒性：Δ=0（输入停滞）、Δ 变号、Δ→0 的连续性；
6. 一致性：与已验证的振荡器实现 `aaiir.aa_iir_osc` 在同一输入下逐点一致。

用法: python qwqdsp/labs/adaa_iir/check_ws.py
"""
from __future__ import annotations

import sys

import numpy as np
import scipy.integrate as si
import scipy.signal as ss

import aaiir as A
import ws as W

FS = 44100.0
failures = []


def report(name: str, err: float, tol: float) -> None:
    ok = err <= tol
    print(f"[{'PASS' if ok else 'FAIL'}] {name:<56} max|diff| = {err:.3e} (tol {tol:.0e})")
    if not ok:
        failures.append(name)


def report_pred(name: str, ok: bool, detail: str = "") -> None:
    print(f"[{'PASS' if ok else 'FAIL'}] {name:<56} {detail}")
    if not ok:
        failures.append(name)


# ------------------------------------------------------------
# 1. 作者参考代码 aaiir_hard_complex.m 的逐行移植
# ------------------------------------------------------------

def _clip(x):
    return min(max(x, -1.0), 1.0)


def _int_ref(x_z1, x, s):
    """参考代码里的 `int(x_z1, x, s)`（四分五裂的分支照抄）。"""
    if x_z1 < x:
        if x <= -1:
            return (1 - np.exp(s)) / s
        if -1 < x <= 1 and x_z1 <= -1:
            return ((x - x_z1) * np.exp(s * (x + 1) / (x - x_z1))
                    - (s + 1) * x + x_z1 - s * np.exp(s)) / s ** 2
        if x > 1 and x_z1 <= -1:
            return ((x - x_z1) * (np.exp(s * (x + 1) / (x - x_z1))
                                  - np.exp(s * (x - 1) / (x - x_z1)))
                    - s * (np.exp(s) + 1)) / s ** 2
        if x_z1 > -1 and x <= 1:
            return (np.exp(s) * (x + (s - 1) * x_z1) - (s + 1) * x + x_z1) / s ** 2
        if -1 < x_z1 <= 1 and x > 1:
            return ((x_z1 - x) * np.exp(s * (x - 1) / (x - x_z1))
                    + np.exp(s) * (x + (s - 1) * x_z1) - s) / s ** 2
        return (np.exp(s) - 1) / s
    else:
        if x_z1 <= -1:
            return (1 - np.exp(s)) / s
        if -1 < x_z1 <= 1 and x <= -1:
            return ((x_z1 - x) * np.exp(s * (x + 1) / (x - x_z1))
                    + np.exp(s) * (x + (s - 1) * x_z1) + s) / s ** 2
        if x_z1 > 1 and x <= -1:
            return ((x_z1 - x) * (np.exp(s * (x + 1) / (x - x_z1))
                                  - np.exp(s * (x - 1) / (x - x_z1)))
                    + s * (np.exp(s) + 1)) / s ** 2
        if x > -1 and x_z1 <= 1:
            return (np.exp(s) * (x + (s - 1) * x_z1) - (s + 1) * x + x_z1) / s ** 2
        if -1 < x <= 1 and x_z1 > 1:
            return ((x - x_z1) * np.exp(s * (x - 1) / (x - x_z1))
                    - (s + 1) * x + x_z1 + s * np.exp(s)) / s ** 2
        return (np.exp(s) - 1) / s


def ref_hard_complex(x, B, beta, tol=1e-9):
    """`aaiir_hard_complex.m` 的等价转录（含 x_z1 初值 0 的约定）。"""
    E = np.exp(beta)
    F = (E - 1.0) / beta
    x_z1, Y_z1 = 0.0, 0j
    y = np.zeros(len(x))
    for n in range(len(x)):
        if abs(x[n] - x_z1) < tol:
            I = _clip(0.5 * (x[n] + x_z1)) * F
        else:
            I = _int_ref(x_z1, x[n], beta)
        Y = E * Y_z1 + 2.0 * B * I
        y[n] = Y.real
        x_z1, Y_z1 = x[n], Y
    return y


# ------------------------------------------------------------
# 2. 连续卷积数值求积（式 (5)），f 任意
# ------------------------------------------------------------

def ref_convolution(x, B, beta, f, n_max: int):
    """``y_n = ∫_0^n 2Re(B e^{βt}) f(x̃(n-t)) dt``（实/虚部各自自适应求积）。"""
    idx = np.arange(len(x), dtype=float)

    def make(n):
        def g(t):
            xi = np.interp(n - t, idx, x)
            return np.exp(beta * t) * W.trivial(xi, f)
        return g

    def cquad(g, a, b):
        return (si.quad(lambda t: g(t).real, a, b, epsabs=1e-13, epsrel=1e-13, limit=300)[0]
                + 1j * si.quad(lambda t: g(t).imag, a, b, epsabs=1e-13, epsrel=1e-13, limit=300)[0])

    y = np.zeros(n_max)
    for n in range(n_max):
        val = sum(cquad(make(n), s, s + 1) for s in range(n))
        y[n] = 2.0 * np.real(B * val)
    return y


# ------------------------------------------------------------

def main() -> int:
    clip = W.hard_clip()
    filt1 = W.aa_iir_1()
    filt2 = W.aa_iir_2()
    filt_r = A.aa_filter_butter(1, 2.0 * np.pi * 0.45)
    print(f"AA-IIR-1: {len(filt1[0])} 对极点, 实极点 {len(filt1[1])}; "
          f"AA-IIR-2: {len(filt2[0])} 对极点")

    # ---- 1. 参考实现逐行移植 ----
    print("== 1) 本实现 vs 作者 MATLAB 参考代码（硬削波） ==")
    rng = np.random.default_rng(0)
    n = 300
    sigs = {
        "10·sin(1 kHz)": 10.0 * np.sin(2 * np.pi * 1000.0 * np.arange(n) / FS),
        "白噪声 ×3": 3.0 * rng.standard_normal(n),
        "直流 0.3": np.full(n, 0.3),
        "缓扫（过零多次）": 8.0 * np.sin(2 * np.pi * np.cumsum(np.linspace(20, 4000, n)) / FS),
    }
    for name, x in sigs.items():
        xa = np.concatenate(([0.0], x))
        err = 0.0
        for B, beta in filt2[0]:
            ref = ref_hard_complex(x, B, beta)
            got = W.aa_iir_pwl(xa, ([(B, beta)], [], 0.0), clip)[1:]
            err = max(err, float(np.max(np.abs(ref - got))))
        report(f"ref vs 本实现  {name}", err, 1e-8)

    # ---- 2. 连续卷积数值求积 ----
    print("== 2) 连续卷积数值求积（式 (5)） ==")
    for name, x in (("正弦 ×10", sigs["10·sin(1 kHz)"][:24]), ("过零噪声", sigs["白噪声 ×3"][:24])):
        for fname, filt in (("AA-IIR-1", filt1), ("butter-1（实极点）", filt_r)):
            if filt[0]:
                B, beta = filt[0][0]
                y_ref = ref_convolution(x, B, beta, clip, len(x))
            else:
                A_, alpha = filt[1][0]
                y_ref = ref_convolution(x, A_ / 2.0, alpha + 0j, clip, len(x))
            y_got = W.aa_iir_pwl(x, filt, clip)
            report(f"quadrature vs 本实现  {name}, {fname}",
                   float(np.max(np.abs(y_ref - y_got))), 1e-7)

    # ---- 3. 原函数 ----
    print("== 3) 原函数 F1/F2 vs 解析积分 ==")
    for name, f in (("hard_clip", clip), ("wavefold", W.wavefold(0.7)),
                    ("algebraic_approx(16)", W.algebraic_approx(16))):
        for order in (1, 2):
            xs = np.linspace(f.X[0] - 2.0, f.X[-1] + 2.0, 37)
            got = f.antideriv(xs, order)
            ref = np.array([
                si.quad(lambda v: float(np.atleast_1d(f.antideriv(np.array([v]), order - 1))[0])
                        if order == 2 else float(np.atleast_1d(f.value(np.array([v])))[0]),
                        f.X[0], xx, epsabs=1e-12, epsrel=1e-12, limit=400,
                        points=list(f.X))[0]
                for xx in xs])
            err = float(np.max(np.abs((got - got[0]) - (ref - ref[0]))))
            report(f"antideriv{order}  {name}", err, 1e-8)

    # ---- 4. 附录 B 的线性恒等式 ----
    print("== 4) 线性情形 f(x)=x：AA-IIR ≡ 冲激不变一阶低通（论文附录 B） ==")
    A_, alpha = filt_r[1][0]
    bz = (-1.0 / alpha) * np.array([np.exp(alpha) - alpha - 1.0,
                                    (alpha - 1.0) * np.exp(alpha) + 1.0])
    az = np.array([1.0, -np.exp(alpha)])
    lin = W.Pwl([-1e6, 1e6], [1.0], [0.0])
    x = rng.standard_normal(2000)
    # 本实现的 y[0] 是递推的占位（第 0 个区间不存在）；参考 H(z) 的 y[0] 已用「x[-1]=0」。
    # 故前置一个 0 采样对齐（与 check 第 1 节同一约定）。
    y_lin = W.aa_iir_pwl(np.concatenate(([0.0], x)), filt_r, lin)
    y_ref = ss.lfilter(bz, az, x)
    report("AA-IIR(f=x)  vs  H(z)（附录 B）", float(np.max(np.abs(y_lin[1:] - y_ref))), 1e-9)

    # ---- 5. 鲁棒性 ----
    print("== 5) 鲁棒性：Δ=0 / Δ 变号 / Δ→0 连续 ==")
    x_stall = np.concatenate(([0.0], np.tile(np.linspace(0, 10, 50), 4)))
    y_stall = W.aa_iir_pwl(x_stall, filt2, clip)
    report_pred("停滞/过零输入无 NaN/Inf",
                bool(np.all(np.isfinite(y_stall))), f"max|y| = {np.abs(y_stall).max():.3f}")
    base = 5.0 * np.sin(2 * np.pi * 3000.0 * np.arange(400) / FS)
    y0 = W.aa_iir_pwl(base, filt1, clip)
    eps = 1e-9
    y1 = W.aa_iir_pwl(base + eps, filt1, clip)
    report("Δ→0 连续（扰动 ε=1e-9）", float(np.max(np.abs(y1 - y0))), 1e-8)

    # ---- 6. 与振荡器实现一致 ----
    print("== 6) 与已验证的振荡器实现一致（同一锯齿输入） ==")
    f0, m = 1000.0, 64
    xr = A.phase_ramp(m, f0, FS)
    saw = A.saw_table()
    # 把周期锯齿写成一个非周期 Pwl（覆盖 0..m*Δ 的每个周期）
    nper = int(np.ceil(xr[-1])) + 1
    X, mm, qq = [], [], []
    for c in range(nper):
        X.append(float(c))
        mm.append(2.0)
        qq.append(-1.0 - 2.0 * c)
    X.append(float(nper))
    pwl_saw = W.Pwl(X, mm, qq)
    y_osc = A.aa_iir_osc(xr, filt1, saw)
    y_ws = W.aa_iir_pwl(xr, filt1, pwl_saw)
    report("ws vs aaiir.aa_iir_osc（同一锯齿）", float(np.max(np.abs(y_osc - y_ws))), 1e-9)

    # ---- 7. AA-FIR-2 与三角核卷积 ----
    print("== 7) AA-FIR-2 vs 三角核连续卷积（等步长应严格一致） ==")
    xc = np.linspace(-2.0, 3.0, 40)
    y_fir = W.aa_fir2(xc, clip)
    K = W.mean_integral_pwl  # noqa: F841  (仅提示：下面用求积做独立对照)
    y_ref = _tri_kernel_conv(xc, clip, 40)
    report("等步长 Δ：AA-FIR-2 vs 三角核卷积", float(np.max(np.abs(y_fir[3:] - y_ref[3:]))), 1e-7)

    # ---- 8. 线性情形解析传递函数 vs 时域实现 ----
    print("== 8) 线性情形解析 H(e^{jω}) vs 时域实现（复状态 → 取实部） ==")
    lin = W.Pwl([-1e6, 1e6], [1.0], [0.0])
    n_t, skip = 20000, 6000
    for fname in ("AA-IIR-1", "AA-IIR-2"):
        filt = W.aa_iir_1() if fname == "AA-IIR-1" else W.aa_iir_2()
        err = 0.0
        for f0 in (100.0, 1000.0, 5000.0, 10000.0):
            w0 = 2.0 * np.pi * f0 / FS
            tt = np.arange(n_t)
            yc = W.aa_iir_pwl(np.cos(w0 * tt), filt, lin)[skip:]
            ys = W.aa_iir_pwl(np.sin(w0 * tt), filt, lin)[skip:]
            ref = np.exp(1j * w0 * tt[skip:])
            got = np.mean((yc + 1j * ys) * np.conj(ref))
            err = max(err, float(abs(got - W.linear_transfer(filt, np.array([w0]))[0])))
        report(f"解析 H vs 时域  {fname}", err, 1e-6)

    print()
    if failures:
        print(f"FAILED: {len(failures)} 项 -> {failures}")
        return 1
    print("全部通过。")
    return 0


def _tri_kernel_conv(x, f, n_max):
    """三角核（宽 2 采样）对线性插值信号的连续卷积（等步长输入的独立对照）。"""
    idx = np.arange(len(x), dtype=float)
    y = np.zeros(n_max)
    for n in range(2, n_max):
        def g(t, n=n):
            xi = np.interp(n - t, idx, x)
            return (1.0 - abs(t - 1.0)) * W.trivial(xi, f)
        y[n] = sum(si.quad(g, s, s + 1, epsabs=1e-13, epsrel=1e-13, limit=200)[0] for s in (0, 1))
    return y


if __name__ == "__main__":
    sys.exit(main())
