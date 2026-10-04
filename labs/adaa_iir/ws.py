"""ADAA-IIR 用于静态非线性（waveshaper）——探索实现。

对应论文
--------
P. P. La Pastina, S. D'Angelo, L. Gabrielli, "Arbitrary-Order IIR Antiderivative
Antialiasing", Proc. 24th Int. Conf. Digital Audio Effects (DAFx20in21), 2021, pp. 9-16.
官方参考实现: https://dangelo.audio/dafx21-aa-iir （硬削波 + 复共轭极点对）

与振荡器版（`aaiir.py`，TASLP 那篇）的本质区别
---------------------------------------------
振荡器里非线性函数的输入是「时间/相位」，每采样前进量 Δ 恒定且单调；waveshaper 的输入是
音频信号本身，Δ = x_{n+1} - x_n **逐样本变化、可正可负、可为零**。因此：

* 积分核依赖 Δ（``γ = β/Δ``），**无法**像 FIR 版（ADAA-LUT，DAFx25）那样预计算一张与 Δ
  无关的原函数表 —— 这是把 ADAA-IIR 用到 waveshaper 上最核心的工程约束；
* 本实现改用**区间参数化** u ∈ [0,1]（ξ = x_n + uΔ、权重 e^{βw}，w = 1-u）：

      I_n = ∫_0^1 f(x_n + uΔ) · e^{β(1-u)} du

  Δ 只通过 ``α = -m·Δ`` 进入斜率项，指数部分的自变量始终有界（``βw``, w ∈ [0,1]），
  于是 Δ→0 时数值不炸（Δ=0 单独取极限），且与作者参考代码逐点一致（见 `check_ws.py`）。

两条计算路线
------------
* **Pwl（分段线性，闭式）**：段内 f 是直线，积分有解析式

      ∫_0^h (γ₀ + αz)e^{βz}dz = γ₀(e^{βh}-1)/β + α((βh-1)e^{βh}+1)/β²

  代价 ∝ 每个采样区间跨过的段数。硬削波/波折叠只有 2 个折点；任意波表（LUT）则可能几十上百，
  这正是 [19] §III-F 说的 ``l = i_max - i_min`` 代价，在 waveshaper 上比振荡器严重得多。
* **任意可调用 f + Gauss-Legendre 求积**：区间内被积函数光滑，固定点数即可高精度，
  代价与段数无关（适合 tanh 这类无闭式原函数的 f）。

用法见 `exp_ws_*.py` / `check_ws.py`。
"""

from __future__ import annotations

import numpy as np
import scipy.signal as ss

import aaiir as A


# ============================================================
# 分段线性非线性函数（定义域外：钳位 clamp 或线性外推）
# ============================================================

class Pwl:
    """分段线性 f，定义域 [X[0], X[-1]]；段 j 覆盖 [X[j], X[j+1]]，f = m[j]x + q[j]。"""

    def __init__(self, X, m, q, clamp: bool = True, left=None, right=None):
        """
        ``left`` / ``right``：``(m, q)`` 元组，显式指定定义域外的延长线
        （默认：``clamp=True`` 用水平线，否则用端段直线）。波折叠这类函数需要它
        （域外的斜率是折回段，与域内端段不同）。
        """
        X = np.asarray(X, dtype=float)
        m = np.asarray(m, dtype=float)
        q = np.asarray(q, dtype=float)
        assert len(X) == len(m) + 1 and np.all(np.diff(X) > 0), "断点需严格递增"
        self.X, self.m, self.q = X, m, q
        self.clamp = clamp
        if left is not None:
            self.ml, self.ql = float(left[0]), float(left[1])
        elif clamp:
            self.ml, self.ql = 0.0, float(m[0] * X[0] + q[0])
        else:
            self.ml, self.ql = float(m[0]), float(q[0])
        if right is not None:
            self.mr, self.qr = float(right[0]), float(right[1])
        elif clamp:
            self.mr, self.qr = 0.0, float(m[-1] * X[-1] + q[-1])
        else:
            self.mr, self.qr = float(m[-1]), float(q[-1])
        self.xe = np.concatenate(([-np.inf], X, [np.inf]))
        self.me = np.concatenate(([self.ml], m, [self.mr]))
        self.qe = np.concatenate(([self.ql], q, [self.qr]))
        self._c1, self._c2 = self._cum()

    @property
    def n_seg(self) -> int:
        return len(self.me)

    def seg_index(self, x):
        """x 所在（扩展）段的索引，供积分路径切分使用。"""
        return np.clip(np.searchsorted(self.xe, x, side="right") - 1, 0, self.n_seg - 1)

    def value(self, x):
        x = np.asarray(x, dtype=float)
        i = self.seg_index(x)
        return self.me[i] * x + self.qe[i]

    def _cum(self):
        """域内断点处的 P1 = ∫_{X[0]}^x f、P2 = ∫_{X[0]}^x P1。"""
        X, m, q = self.X, self.m, self.q
        c1 = np.zeros(len(X))
        c2 = np.zeros(len(X))
        for j in range(len(m)):
            a, b = X[j], X[j + 1]
            c1[j + 1] = c1[j] + m[j] * (b * b - a * a) / 2.0 + q[j] * (b - a)
            c2[j + 1] = (c2[j] + c1[j] * (b - a)
                         + m[j] * ((b ** 3 - a ** 3) / 3.0 - a * a * (b - a)) / 2.0
                         + q[j] * ((b * b - a * a) / 2.0 - a * (b - a)))
        return c1, c2

    def antideriv(self, x, order: int = 1):
        """``order`` 阶连续原函数 F（归一化 F(0) = 0）；段内解析积分，域外按延长线解析。

        ⚠ ``_eval_raw`` 给出的是以 ``X[0]`` 为参考的原函数 ``P₁/P₂``；平移到 F(0)=0 时，
        ``F₂`` 必须同时补掉一次项：``F₂(x) = P₂(x) - P₂(0) - P₁(0)·x``
        （否则 ``F₂' ≠ F₁``，AA-FIR-2 会错）。
        """
        if not hasattr(self, "_f0"):
            self._f0 = (float(self._eval_raw(np.array([0.0]), 1)[0]),
                        float(self._eval_raw(np.array([0.0]), 2)[0]))
        x = np.asarray(x, dtype=float)
        if order == 1:
            return self._eval_raw(x, 1) - self._f0[0]
        if order == 2:
            return self._eval_raw(x, 2) - self._f0[1] - self._f0[0] * x
        raise ValueError("只实现 1/2 阶原函数")

    def _eval_raw(self, x, order):
        """未归一化的原函数（P(X[0]) = 0）。"""
        x = np.asarray(x, dtype=float)
        X, m, q, c1, c2 = self.X, self.m, self.q, self._c1, self._c2
        K = len(m)
        j = np.clip(np.searchsorted(X, x, side="right") - 1, 0, K - 1)
        a = X[j]

        def d1(a_, mm, qq):
            return mm * (x * x - a_ * a_) / 2.0 + qq * (x - a_)

        def d2(a_, mm, qq):
            return (mm * ((x ** 3 - a_ ** 3) / 3.0 - a_ * a_ * (x - a_)) / 2.0
                    + qq * ((x * x - a_ * a_) / 2.0 - a_ * (x - a_)))

        inside = (x >= X[0]) & (x <= X[-1])
        v1 = np.where(inside, c1[j] + d1(a, m[j], q[j]), 0.0)
        v2 = np.where(inside, c2[j] + c1[j] * (x - a) + d2(a, m[j], q[j]), 0.0)

        # 右端延长：从 X[-1] 起用右延长线
        xr = X[-1]
        e1r = self.mr * (x * x - xr * xr) / 2.0 + self.qr * (x - xr)
        e2r = (self.mr * ((x ** 3 - xr ** 3) / 3.0 - xr * xr * (x - xr)) / 2.0
               + self.qr * ((x * x - xr * xr) / 2.0 - xr * (x - xr)))
        right = x > X[-1]
        v1 = np.where(right, c1[-1] + e1r, v1)
        v2 = np.where(right, c2[-1] + c1[-1] * (x - xr) + e2r, v2)

        # 左端延长：P1(x) = -∫_x^{X[0]} f，P2(x) = -∫_x^{X[0]} P1
        xl = X[0]
        left = x < X[0]
        i1 = self.ml * (xl * xl - x * x) / 2.0 + self.ql * (xl - x)
        i2 = (self.ml / 2.0 * (xl * xl * (xl - x) - (xl ** 3 - x ** 3) / 3.0)
              + self.ql * (xl - x) ** 2 / 2.0)
        v1 = np.where(left, -i1, v1)
        v2 = np.where(left, i2, v2)

        if order == 1:
            return v1
        if order == 2:
            return v2
        raise ValueError("只实现 1/2 阶原函数")


# ------------------------------------------------------------
# 常见非线性
# ------------------------------------------------------------

def hard_clip(threshold: float = 1.0) -> Pwl:
    """硬削波 f(x) = clip(x, ±T)（[19] 式 (15)），2 个折点。"""
    return Pwl([-threshold, threshold], [1.0], [0.0])


def wavefold(tau: float = 0.7) -> Pwl:
    """[DAFx25] 式 (15) 的非对称波折叠：|x|≤τ 时 f=x，之外折回（f = ±2τ - x）。"""
    return Pwl([-tau, tau], [1.0], [0.0],
               left=(-1.0, -2.0 * tau), right=(-1.0, 2.0 * tau))


def algebraic_approx(k: int = 16) -> Pwl:
    """仓库 `AlgebraicWaveshaper` 的 f(x)=x/sqrt(1+x²) 的 k 段线性近似（便于对照）。"""
    return fit_from_callable(lambda v: v / np.sqrt(1.0 + v * v), -8.0, 8.0, k)


def fit_from_callable(fn, lo: float, hi: float, k: int, clamp: bool = True) -> Pwl:
    """把任意 f 用 k 段等距线性插值拟合成 Pwl（波表 / 查表路线）。"""
    X = np.linspace(lo, hi, k + 1)
    Y = fn(X)
    m = np.diff(Y) / np.diff(X)
    q = Y[:-1] - m * X[:-1]
    return Pwl(X, m, q, clamp=clamp)


# ============================================================
# AA-IIR：区间参数化的加权积分
# ============================================================

def mean_integral_pwl(x, beta, f: Pwl):
    """``I_n = ∫_0^1 f(x_n+uΔ)e^{β(1-u)}du``（闭式，Pwl），返回长度 len(x)-1。"""
    x = np.asarray(x, dtype=float)
    if len(x) < 2:
        return np.zeros(0, dtype=complex)
    x0, x1 = x[:-1], x[1:]
    delta = x1 - x0
    static = np.abs(delta) < 1e-300
    safe = np.where(static, 1.0, delta)
    lo = np.minimum(x0, x1)
    hi = np.maximum(x0, x1)
    i_lo = f.seg_index(lo)
    n_slot = (f.seg_index(hi) - i_lo).astype(np.int64) + 1
    acc = np.zeros(len(x0), dtype=complex)
    for k in range(int(n_slot.max())):
        mask = (k < n_slot) & ~static
        if not mask.any():
            break
        i = np.minimum(i_lo + k, f.n_seg - 1)
        a = np.maximum(lo, f.xe[i])
        b = np.minimum(hi, f.xe[i + 1])
        ua = (a - x0) / safe
        ub = (b - x0) / safe
        u_hi = np.maximum(ua, ub)
        u_lo = np.minimum(ua, ub)
        w_a = 1.0 - u_hi                       # 槽位起点到区间末端的采样距离
        h = u_hi - u_lo                        # 槽位在采样尺度上的长度 ∈ [0,1]
        # 被 mask 掉的槽位数值无意义（可能极大/负），先钳到安全值再算，避免 inf/nan
        w_a = np.where(mask, w_a, 0.0)
        h = np.where(mask, h, 0.0)
        xi_start = x0 + u_hi * delta           # 槽位起点的 ξ（沿路径）
        gamma0 = f.me[i] * xi_start + f.qe[i]
        alpha = -f.me[i] * delta
        E = np.exp(beta * h)
        integ = gamma0 * (E - 1.0) / beta + alpha * ((beta * h - 1.0) * E + 1.0) / (beta * beta)
        acc += np.where(mask, integ * np.exp(beta * w_a), 0.0)
    # Δ = 0：f 在区间上是常数，取极限 I = f(x)·(e^β-1)/β
    if static.any():
        acc = np.where(static, f.value(x0) * (np.exp(beta) - 1.0) / beta, acc)
    return acc


def mean_integral_numeric(x, beta, fn, order: int = 8):
    """任意 f + Gauss-Legendre 求积（自适应于 Δ，无除零问题）。"""
    x = np.asarray(x, dtype=float)
    if len(x) < 2:
        return np.zeros(0, dtype=complex)
    gx, gw = np.polynomial.legendre.leggauss(order)
    u = 0.5 * (gx + 1.0)
    w = 0.5 * gw
    x0, x1 = x[:-1], x[1:]
    delta = (x1 - x0)[:, None]
    xi = x0[:, None] + u[None, :] * delta
    weight = np.exp(beta * (1.0 - u))[None, :]
    return (fn(xi) * weight * w[None, :]).sum(axis=1)


def aa_iir(x, filters, integral):
    """AA-IIR waveshaper。``filters`` = (pairs, reals, direct)；``integral`` 算 I_n。"""
    pairs, reals, _direct = filters
    x = np.asarray(x, dtype=float)
    y = np.zeros(len(x))
    if len(x) < 2:
        return y
    for B, beta in pairs:
        I = integral(x, beta)
        v = np.zeros(len(x), dtype=complex)
        v[1:] = 2.0 * B * I
        y += np.real(ss.lfilter([1.0], [1.0, -np.exp(beta)], v))
    for A_, alpha in reals:
        I = integral(x, alpha + 0j)
        v = np.zeros(len(x))
        v[1:] = A_ * I.real
        y += ss.lfilter([1.0], [1.0, -np.exp(alpha)], v)
    return y


def aa_iir_pwl(x, filters, f: Pwl):
    """AA-IIR（分段线性非线性，闭式积分）。"""
    return aa_iir(x, filters, lambda xx, beta: mean_integral_pwl(xx, beta, f))


def aa_iir_fn(x, filters, fn, order: int = 8):
    """AA-IIR（任意 f，数值求积）。"""
    return aa_iir(x, filters, lambda xx, beta: mean_integral_numeric(xx, beta, fn, order))


def linear_transfer(filters, w):
    """线性情形 f(x)=x 时 AA-IIR 的**精确**传递函数 ``H(e^{jω})``（ω = w [rad/sample]）。

    对 f(x)=x：``I_n = A₀·x_n + A₁·(x_{n+1}-x_n)``，其中
    ``A₀ = (e^β-1)/β``、``A₁ = (e^β-β-1)/β²``，于是每个极点贡献一个一阶节

        ``Hʹ_β(z) = 2B·[A₁ + (A₀-A₁)z^{-1}] / (1 - e^β z^{-1})``

    （实极点把 ``2B`` 换成 ``A``）。注意 ``Hʹ`` 是**复状态** ŷ 的响应，实输出取
    ``y = Re(ŷ)``，故真实系统对实数输入的响应是

        ``H(e^{jω}) = [ Hʹ(e^{jω}) + conj(Hʹ(e^{-jω})) ] / 2``

    （直流处这才给出增益 1；直接用 ``Hʹ`` 会得到 ``1-j``、模 √2。）
    ``check_ws.py`` 第 8 节把它与时域实现逐点比对。
    """
    w = np.asarray(w, dtype=float)

    def raw(z):
        h = np.zeros(len(z), dtype=complex)
        for B, beta in pairs_:
            a0 = (np.exp(beta) - 1.0) / beta
            a1 = (np.exp(beta) - beta - 1.0) / (beta * beta)
            h += 2.0 * B * (a1 + (a0 - a1) / z) / (1.0 - np.exp(beta) / z)
        for A_, alpha in reals_:
            a0 = (np.exp(alpha) - 1.0) / alpha
            a1 = (np.exp(alpha) - alpha - 1.0) / (alpha * alpha)
            h += A_ * (a1 + (a0 - a1) / z) / (1.0 - np.exp(alpha) / z)
        return h

    pairs_, reals_, _direct = filters
    return 0.5 * (raw(np.exp(1j * w)) + np.conj(raw(np.exp(-1j * w))))


# ============================================================
# 对照方法
# ============================================================

def trivial(x, f):
    return f.value(x) if isinstance(f, Pwl) else f(x)


def oversampled(x, factor: int, f, zero_phase: bool = True):
    """OVS-r：线性插值上采样 → 过非线性 → 8 阶 Chebyshev I 抽取（同 [19] 描述）。

    ``zero_phase=True`` 用 ``sosfiltfilt``（零相位，幅度响应平方 → 阻带衰减翻倍）；
    ``False`` 用单程因果 ``sosfilt``（更接近「一路滤波+下采样」的朴素实现）。
    [19] 未说明抽取实现，两者的 SNR 差别足以改变 OVS-2 与 AA-FIR-1 的排序（见 README）。
    """
    n = len(x)
    idx = np.arange(n, dtype=float)
    fine = np.arange(0, n - 1, 1.0 / factor)
    up = np.interp(fine, idx, x)
    y = f.value(up) if isinstance(f, Pwl) else f(up)
    sos = ss.cheby1(8, 0.05, 0.8 / factor, output="sos")
    yf = ss.sosfiltfilt(sos, y) if zero_phase else ss.sosfilt(sos, y)
    return yf[::factor][:n]


def aa_fir1(x, f):
    """一阶 ADAA（矩形核）：``(F1(x_n)-F1(x_{n-1}))/(x_n-x_{n-1})``，退化点取 f(中点)。"""
    F1 = f.antideriv(x, 1)
    dx = np.diff(x)
    y = np.zeros(len(x))
    ok = np.abs(dx) > 1e-12
    q = np.where(ok, np.diff(F1) / np.where(ok, dx, 1.0), 0.0)
    mid = 0.5 * (x[1:] + x[:-1])
    y[1:] = np.where(ok, q, trivial(mid, f))
    y[0] = y[1] if len(y) > 1 else 0.0
    return y


def aa_fir2(x, f):
    """二阶 ADAA（三角核，等价于两个矩形核的复合），任意（含变步长）输入：

        ``s_i = (F₂(x_i) - F₂(x_{i-1})) / (x_i - x_{i-1})``
        ``y_n = 2·(s_n - s_{n-1}) / (x_n - x_{n-2})``

    等步长 Δ 时化为 ``(F₂(x_n) - 2F₂(x_{n-1}) + F₂(x_{n-2}))/Δ²``。

    ⚠ 用的是**二阶**原函数 F₂ = ∫F₁（不是 F₁）。用 F₁ 代入会得到一个「求导型」算子
    （对 f(x)=x 输出常数 1 而非 x），与三角核卷积完全不符 —— `check_ws.py` 第 7 节
    对等步长情形做了逐点求积对照。
    """
    F2 = f.antideriv(x, 2)
    n = len(x)
    y = np.zeros(n)
    if n < 3:
        return aa_fir1(x, f)
    dx1 = x[1:-1] - x[:-2]
    dx2 = x[2:] - x[1:-1]
    dxt = x[2:] - x[:-2]
    ok1 = np.abs(dx1) > 1e-12
    ok2 = np.abs(dx2) > 1e-12
    okt = np.abs(dxt) > 1e-12
    s1 = np.where(ok1, (F2[1:-1] - F2[:-2]) / np.where(ok1, dx1, 1.0), 0.0)
    s2 = np.where(ok2, (F2[2:] - F2[1:-1]) / np.where(ok2, dx2, 1.0), 0.0)
    val = np.where(okt, 2.0 * (s2 - s1) / np.where(okt, dxt, 1.0), 0.0)
    # 退化（连续样本相等）：退回一阶结果
    fallback = aa_fir1(x, f)
    y[2:] = np.where(ok1 & ok2 & okt, val, fallback[2:])
    y[:2] = y[2]
    return y


# ============================================================
# 滤波器：[19] 的参考代码
# ============================================================

def aa_iir_1():
    """[19] AA-IIR-1：2 阶 Butterworth，Fc = 0.45·Fs。"""
    return A.aa_filter_butter(2, 2.0 * np.pi * 0.45)


def aa_iir_2():
    """[19] AA-IIR-2：10 阶 Chebyshev II，rs = 60 dB，Fbs = 0.61·Fs。

    ⚠ [19] 的 `aaiir_sweep.m` 用 ``2*pi*0.61``（阻带边沿 0.61·Fs，在 Nyquist 之上），
    而振荡器篇 demo 用 ``pi*0.61``（= 0.61·Nyquist）。这里按 [19]。
    """
    return A.aa_filter_cheby2(10, 60.0, 2.0 * np.pi * 0.61)


# ============================================================
# 误差度量与解析锚点
# ============================================================

def clip_fourier(amplitude: float, threshold: float, kmax: int):
    """硬削波正弦的解析傅里叶系数 A_k（奇次谐波，与输入同相）。

    用**等间隔采样的 FFT** 求系数：被削波正弦是分段光滑的周期函数，M 点均匀采样 + DFT
    对 k < M/2 的系数是精确的（折回的尾巴 ~1/M 量级，可忽略）。⚠ 不能用固定点数的
    Gauss-Legendre：k 上千次振荡时高阶求积会退化成乱数（踩过）。
    """
    m = 1
    while m < max(4 * kmax, 1 << 15):
        m <<= 1
    t = np.arange(m) / m
    y = np.clip(amplitude * np.sin(2.0 * np.pi * t), -threshold, threshold)
    spec = np.fft.rfft(y)
    ks = np.arange(1, kmax + 1)
    ak = 2.0 * np.abs(spec[1:kmax + 1]) / m
    return ks, ak


def analytic_clip_snr(f0: float, fs: float, amplitude: float = 10.0, threshold: float = 1.0):
    """trivial 硬削波方法的解析 SNR：带内谐波功率 / 折叠分量功率（S 的独立锚点）。"""
    kmax = int(np.ceil(fs / 2.0 / f0)) + 4000
    ks, ak = clip_fourier(amplitude, threshold, kmax)
    p = ak ** 2 / 2.0
    inband = ks * f0 < fs / 2.0
    return 10.0 * np.log10(p[inband].sum() / p[~inband].sum())
