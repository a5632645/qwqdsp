"""AA-IIR 抗混叠任意波形振荡器 —— 论文算法的最小复现实现。

对应论文
--------
L. Gabrielli, S. D'Angelo, P. P. La Pastina, S. Squartini,
"Antiderivative Antialiasing for Arbitrary Waveform Generation",
IEEE/ACM Trans. Audio, Speech, Lang. Process., 30:2743-2753, 2022.
doi:10.1109/TASLP.2022.3198007
官方参考 MATLAB 代码: https://dangelo.audio/taslp-antialias-waveform

约定（与作者 MATLAB 参考代码一致）
--------------------------------
* 相位 ``x`` 以**周期**为单位，每采样前进 ``delta = f0 / fs``；周期 ``T = X[-1]``。
* 抗混叠滤波器 ``H(s)`` 用「每采样」归一化：极点 ``beta`` 本身就是每采样衰减率
  （``Re(beta) < 0``），递归因子 ``exp(beta)``，无需乘 ``1/fs``。
  这来自参考代码的极点归一化（例如 ``butter(2, 2*pi*0.45, 's')``，见
  :func:`aa_filter_butter` / :func:`aa_filter_cheby2` 的注释）。
* 分段线性周期函数 ``f`` 用 (X, m, q) 表示：段 ``j``（0-based）覆盖
  ``[X[j], X[j+1]]``，``f(x) = m[j]*x + q[j]``，``X`` 等距（间隔 ``T/k``）。
* 式号引用论文正文（Eq. (n)）。

实现说明
--------
论文式 (20)–(25)（简单实极点）与式 (26)（共轭极点对）是**逐采样递推**的。
本实现把它拆成两步：

1. 每个采样区间上的积分 ``I_n`` 用解析闭式、**向量化**算出来（式 (22)–(25) 的等价
   形式，见 :func:`aa_iir_forcing`）；
2. 递推 ``y_hat_{n+1} = exp(beta)*y_hat_n + 2B*I_n`` 用 ``scipy.signal.lfilter``
   一次算完（一阶线性递推，与逐点循环逐位等价）。

第一步的闭式（段 [a, b] 上 ``f(xi) = m*xi + c``，``gamma = beta/delta``）：

    ∫_a^b (m*xi + c) * exp(gamma*(b - xi)) dxi
      = (m*b + c)*(E - 1)/gamma - m*(E*(gamma*h - 1) + 1)/gamma^2,  E = exp(gamma*h)

它由式 (22) 的 ``F(m, q, ·)`` 差商 + 跨段求和折叠而来；等价性由
``check_impl.py`` 对作者参考代码逐点核对。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.signal as ss
import scipy.special as sp


# ============================================================
# 分段线性周期波形（表）
# ============================================================

@dataclass
class PwlTable:
    """分段线性周期函数 f(x)。

    段 ``j``（0-based，共 ``k`` 段）覆盖 ``[X[j], X[j+1]]``，``f(x) = m[j]*x + q[j]``；
    ``X`` 必须等距（间隔 ``T/k``，线性插值表天然满足）。
    ``wt`` 是均匀采样的整表（供 trivial / 过采样方法线性插值读取），与 (X, m, q)
    描述同一波形。
    """

    X: np.ndarray
    m: np.ndarray
    q: np.ndarray
    wt: np.ndarray | None = None

    def __post_init__(self) -> None:
        self.X = np.asarray(self.X, dtype=float)
        self.m = np.asarray(self.m, dtype=float)
        self.q = np.asarray(self.q, dtype=float)
        if self.wt is not None:
            self.wt = np.asarray(self.wt, dtype=float)
        assert self.X.ndim == 1 and self.m.shape == self.q.shape
        assert len(self.X) == len(self.m) + 1, "断点数应为段数 + 1"
        step = np.diff(self.X)
        assert np.allclose(step, step[0]), "只支持等距断点（线性插值表）"

    @property
    def T(self) -> float:
        """周期。"""
        return float(self.X[-1])

    @property
    def k(self) -> int:
        """每周期段数。"""
        return int(len(self.m))


# ------------------------------------------------------------
# 内置波形
# ------------------------------------------------------------

def saw_table() -> PwlTable:
    """理想锯齿波 f(x) = 2*x - 1, x ∈ [0, 1)（论文式 (1)，k = 1）。"""
    n = 2048
    phase = np.arange(n) / n
    return PwlTable(X=[0.0, 1.0], m=[2.0], q=[-1.0], wt=2.0 * phase - 1.0)


def escalation_ii_w3_table() -> PwlTable:
    """Massive "Escalation II" 第 3 张表（论文 §IV-C 用到的波形）。

    数值与作者 MATLAB 参考代码 ``generateEscalationII_w3.m`` 一致：
    ``X = [0..8]/8``，``m = [1 -1 1 -1 -1 1 -1 1]*2``，``q = [0 0 0 0 2 -2 2 -2]``，
    整表 2048 点。
    """
    X = np.arange(9) / 8.0
    m = np.array([1, -1, 1, -1, -1, 1, -1, 1], dtype=float) * 2.0
    q = np.array([0, 0, 0, 0, 2, -2, 2, -2], dtype=float)
    seg = 2048 // len(m)
    wt = np.concatenate(
        [m[j] * np.linspace(X[j], X[j + 1], seg) + q[j] for j in range(len(m))]
    )
    return PwlTable(X=X, m=m, q=q, wt=wt)


# ============================================================
# 分段线性表的取值 / 积分
# ============================================================

def _travel(x, tab: PwlTable):
    """把相位分解为 (p, r)：``x = p*T + r``，``r ∈ [0, T)``。"""
    x = np.asarray(x, dtype=float)
    p = np.floor(x / tab.T)
    r = x - p * tab.T
    return p.astype(np.int64), r


def pwl_value(x, tab: PwlTable):
    """按分段解析式取值（周期延拓），用于校验 / 直接项。"""
    _, r = _travel(x, tab)
    idx = np.clip((r * tab.k / tab.T).astype(np.int64), 0, tab.k - 1)
    return tab.m[idx] * r + tab.q[idx]


def lookup(x, tab: PwlTable):
    """均匀采样表的线性插值读取（周期延拓）——trivial / 过采样方法用的读取方式。"""
    assert tab.wt is not None, "该表没有整表采样"
    wt = tab.wt
    n = len(wt)
    pos = np.mod(np.asarray(x, dtype=float) / tab.T, 1.0) * n
    i0 = np.floor(pos).astype(np.int64) % n
    frac = pos - np.floor(pos)
    return (1.0 - frac) * wt[i0] + frac * wt[(i0 + 1) % n]


@dataclass
class _Antiderivatives:
    """周期分段线性 f 的 1/2 阶原函数的查表常数（F(0) = 0）。"""

    cs1: np.ndarray   # F1 在断点处的值
    cs2: np.ndarray   # F2 在断点处的值
    drift1: float     # F1 每周期增量 = ∫_0^T f
    drift2: float     # F2 每周期增量 = ∫_0^T F1


def _antiderivatives(tab: PwlTable) -> _Antiderivatives:
    X, m, q = tab.X, tab.m, tab.q
    seg = [0.0]
    seg1 = [0.0]
    for j in range(tab.k):
        a, b = X[j], X[j + 1]
        c1 = m[j] * (b * b - a * a) / 2.0 + q[j] * (b - a)
        seg.append(seg[-1] + c1)
        # ∫_a^b loc1(u) du，loc1(u) - loc1(a) = m[j]*(u²-a²)/2 + q[j]*(u-a)
        c2 = (
            seg[j] * (b - a)
            + m[j] * ((b**3 - a**3) / 3.0 - a * a * (b - a)) / 2.0
            + q[j] * ((b * b - a * a) / 2.0 - a * (b - a))
        )
        seg1.append(seg1[-1] + c2)
    cs1 = np.array(seg)
    cs2 = np.array(seg1)
    return _Antiderivatives(cs1=cs1, cs2=cs2, drift1=float(cs1[-1]), drift2=float(cs2[-1]))


def antiderivative(x, tab: PwlTable, order: int = 1, _pre: _Antiderivatives | None = None):
    """``f`` 的 ``order`` 阶原函数（``F(0) = 0``，连续）。

    ``order=1`` 返回 ``F1(x) = ∫_0^x f``；``order=2`` 返回 ``F2(x) = ∫_0^x F1``。
    用 ``x = p*T + r`` 的精确分解，避免大 ``x`` 下的抵消误差。
    """
    pre = _pre or _antiderivatives(tab)
    p, r = _travel(x, tab)
    X, m, q, k, T = tab.X, tab.m, tab.q, tab.k, tab.T
    idx = np.clip((r * k / T).astype(np.int64), 0, k - 1)
    a = X[idx]
    # loc1(r) = F1 在当前周期内的值
    loc1 = pre.cs1[idx] + m[idx] * (r * r - a * a) / 2.0 + q[idx] * (r - a)
    if order == 1:
        return p * pre.drift1 + loc1
    if order == 2:
        loc2 = (
            pre.cs2[idx]
            + pre.cs1[idx] * (r - a)
            + m[idx] * ((r**3 - a**3) / 3.0 - a * a * (r - a)) / 2.0
            + q[idx] * ((r * r - a * a) / 2.0 - a * (r - a))
        )
        return (
            p * pre.drift2
            + pre.drift1 * T * p * (p - 1) / 2.0
            + loc2
            + p * pre.drift1 * r
        )
    raise ValueError("只实现 1 阶与 2 阶原函数")


# ============================================================
# 抗混叠滤波器设计（对应论文 §IV-A 的两种滤波器）
# ============================================================

def _filter_from_zpk(z, p, gain):
    """零极点 → 部分分式 (pairs, reals, direct)。

    ``pairs``: [(B, beta)]，每个对应一对共轭极点（取上半平面那个），
    时域冲激响应贡献 ``2*Re(B*exp(beta*t))``；
    ``reals``: [(A, alpha)]，实极点，贡献 ``A*exp(alpha*t)``；
    ``direct``: 直接项 A0（见论文式 (8)）。
    """
    b, a = ss.zpk2tf(z, p, gain)
    r, poles, k = ss.residue(b, a)
    direct = float(np.atleast_1d(k)[0]) if np.size(k) else 0.0
    pairs, reals = [], []
    for pole, res in zip(poles, r):
        if pole.imag > 0:
            pairs.append((complex(res), complex(pole)))
        elif pole.imag == 0 and abs(pole.imag) < 1e-14 and abs(res.imag) < 1e-14:
            reals.append((float(res.real), float(pole.real)))
        # 下半平面极点：其共轭已在 pairs 里，跳过
    return pairs, reals, direct


def aa_filter_butter(order: int, wn_rad: float):
    """模拟 Butterworth 低通原型（``wn_rad`` 为「每采样」角频率，如 2*pi*0.45）。"""
    z, p, k = ss.butter(order, wn_rad, btype="low", analog=True, output="zpk")
    return _filter_from_zpk(z, p, k)


def aa_filter_cheby2(order: int, rs_db: float, wn_rad: float):
    """模拟 Chebyshev II 低通原型（``wn_rad`` 为阻带边沿「每采样」角频率）。"""
    z, p, k = ss.cheby2(order, rs_db, wn_rad, btype="low", analog=True, output="zpk")
    return _filter_from_zpk(z, p, k)


# 论文 §IV-A 的两种滤波器规格。
# ⚠ 论文正文写 "fc = 0.45 fs" / "f_bs = 0.61 fs"，但作者 MATLAB 参考代码里
#   cheby2 传的是 pi*0.61（不是 2*pi*0.61），即阻带边沿实际是 0.61*Nyquist。
#   这里照抄参考代码的数值以复现其 SNR。
AAIIR1_WN = 2.0 * np.pi * 0.45     # 2 阶 Butterworth
AAIIR2_RS = 60.0                   # 10 阶 Chebyshev II 阻带衰减 [dB]
AAIIR2_WN = np.pi * 0.61           # 10 阶 Chebyshev II 阻带边沿


def aa_iir_1():
    """论文 AA-IIR-1：2 阶 Butterworth，fc = 0.45 fs。"""
    return aa_filter_butter(2, AAIIR1_WN)


def aa_iir_2():
    """论文 AA-IIR-2：10 阶 Chebyshev II，rs = 60 dB，fbs(参考代码) = 0.61 Nyquist。"""
    return aa_filter_cheby2(10, AAIIR2_RS, AAIIR2_WN)


# ============================================================
# AA-IIR 振荡器
# ============================================================

def aa_iir_forcing(x, beta, tab: PwlTable):
    """式 (20)/(26) 的积分项 ``I_n``（不含 2B），n = 0..len(x)-2。

    ``I_n = (1/delta_n) ∫_{x_n}^{x_{n+1}} f(xi) exp(beta*(x_{n+1}-xi)/delta_n) dxi``

    ``delta_n = x_{n+1} - x_n``。闭式按「跨段」求和：把每个采样区间按表的断点
    切成若干子段（段内 f 是直线），每子段用解析积分。
    """
    x = np.asarray(x, dtype=float)
    if len(x) < 2:
        return np.zeros(0, dtype=complex)
    delta = np.diff(x)
    if np.any(delta <= 0):
        raise ValueError("相位必须严格递增")
    T, k = tab.T, tab.k
    t = x[:-1] * k / T      # 左端点的「段坐标」
    u = x[1:] * k / T       # 右端点
    base = np.floor(t)
    n_slot = (np.floor(u) - base).astype(np.int64) + 1
    gam = beta / delta      # 复；gamma * h ∈ [0, beta]
    acc = np.zeros_like(delta, dtype=complex)
    for s in range(int(n_slot.max())):
        mask = s < n_slot
        if not mask.any():
            break
        seg_a = np.clip(base + s, t, u)
        seg_b = np.clip(base + s + 1, t, u)
        i_seg = (base + s).astype(np.int64)
        j = np.mod(i_seg, k)                  # 段号（周期内）
        p = np.floor_divide(i_seg, k)         # 周期数
        # 局部坐标 = 段坐标（去整周期），避免大 x 下的抵消
        loc_b = (seg_b - p * k) * (T / k)
        h = (seg_b - seg_a) * (T / k)
        E = np.exp(gam * h)
        term = (
            (tab.m[j] * loc_b + tab.q[j]) * (E - 1.0) / gam
            - tab.m[j] * (E * (gam * h - 1.0) + 1.0) / (gam * gam)
        )
        # 指数权重以**整个采样区间**的右端 x_{n+1} 为基准（跨段时逐段衰减）
        term = term * np.exp(gam * (u - seg_b) * (T / k))
        acc += np.where(mask, term, 0.0)
    return acc / delta


def aa_iir_osc(x, filters, tab: PwlTable, include_direct: bool = False):
    """AA-IIR 振荡器（式 (20)/(26)）。

    ``filters`` 为 :func:`aa_filter_butter` / :func:`aa_filter_cheby2` 的返回值
    ``(pairs, reals, direct)``。返回与 ``x`` 等长的实信号。
    """
    pairs, reals, direct = filters
    x = np.asarray(x, dtype=float)
    y = np.zeros(len(x), dtype=float)
    if len(x) < 2:
        return y
    for B, beta in pairs:
        I = aa_iir_forcing(x, beta, tab)
        v = np.zeros(len(x), dtype=complex)
        v[1:] = 2.0 * B * I
        y += np.real(ss.lfilter([1.0], [1.0, -np.exp(beta)], v))
    for A, alpha in reals:
        I = aa_iir_forcing(x, alpha + 0j, tab)
        v = np.zeros(len(x))
        v[1:] = A * I.real
        y += ss.lfilter([1.0], [1.0, -np.exp(alpha)], v)
    if include_direct and direct != 0.0:
        y += direct * pwl_value(x, tab)
    return y


def phase_ramp(n_samples: int, f0: float, fs: float):
    """``x = [0, delta, 2*delta, ...)``，``delta = f0/fs``。"""
    return np.arange(n_samples, dtype=float) * (f0 / fs)


# ============================================================
# 经典方法（对照）
# ============================================================

def trivial(table: PwlTable, n_samples: int, f0: float, fs: float):
    """trivial：直接采样理想波形（锯齿用式 (1)，wavetable 用线性插值读取）。"""
    x = phase_ramp(n_samples, f0, fs)
    if table.wt is None:
        return pwl_value(x, table)
    return lookup(x, table)


def oversampled(table: PwlTable, n_samples: int, f0: float, fs: float, ratio: int):
    """OVS-r：高倍率生成 trivial 波形，再 ``scipy.signal.decimate`` 抽取。

    与 MATLAB ``decimate`` 同规格：8 阶 Chebyshev I（0.05 dB 纹波）、
    截止 0.8*(fs_in/2)/r、零相位滤波（filtfilt）。
    """
    hi = trivial(table, n_samples * ratio, f0, fs * ratio)
    return ss.decimate(hi, ratio, zero_phase=True)[:n_samples]


def aa_fir(table: PwlTable, n_samples: int, f0: float, fs: float, order: int = 1):
    """AA-FIR（论文 §III-B，矩形核 order=1 / 三角核 order=2），输入取时间。

    order=1: ``y_n = (F1(x_n) - F1(x_{n-1})) / delta``
    order=2: ``y_n = (F2(x_n) - 2*F2(x_{n-1}) + F2(x_{n-2})) / delta^2``
    """
    delta = f0 / fs
    x = phase_ramp(n_samples, f0, fs)
    y = np.zeros(n_samples)
    y[0] = np.nan   # 起点没有前值，最后用第 1 点填
    if order == 1:
        F1 = antiderivative(x, table, 1)
        y[1:] = (F1[1:] - F1[:-1]) / delta
    elif order == 2:
        F2 = antiderivative(x, table, 2)
        y[2:] = (F2[2:] - 2.0 * F2[1:-1] + F2[:-2]) / (delta * delta)
        y[1] = y[2]
    else:
        raise ValueError("AA-FIR 只实现 1/2 阶")
    y[0] = y[1]
    return y


def dpw(table: PwlTable, n_samples: int, f0: float, fs: float, order: int):
    """DPW-N（论文式 (11)，[Välimäki 2010]）：多项式波形 + N-1 阶差分。

    多项式取**分数相位** φ = frac(x) 的周期多项式（周期连续、无漂移）：

    * ``N=2``：``p(φ) = φ² - φ``，``p'' = 2``，``y = Δp/Δ``
    * ``N=3``：``p(φ) = φ³ - 3/2·φ² + φ/2``，``p'' = 3·(2φ-1)``，``y = Δ²p/(3Δ²)``

    ⚠ N=3 若用「不做周期连续性修正」的 ``s³/12``（``s = 2φ-1``），会在每个周期
    跳变处留下幅值 ``1/(12Δ²)`` 的尖峰，与 AA-FIR-2 不再相等（见 README 说明）。
    """
    if table.k != 1:
        raise ValueError("DPW 只对单段（锯齿）波形实现")
    delta = f0 / fs
    phi = np.mod(phase_ramp(n_samples, f0, fs), 1.0)
    y = np.zeros(n_samples)
    if order == 2:
        p = phi * phi - phi
        y[1:] = (p[1:] - p[:-1]) / delta
        y[0] = y[1]
    elif order == 3:
        p = phi**3 - 1.5 * phi * phi + 0.5 * phi
        y[2:] = (p[2:] - 2.0 * p[1:-1] + p[:-2]) / (3.0 * delta * delta)
        y[:2] = y[2]
    else:
        raise ValueError("DPW 只实现 2/3 阶")
    return y


# ============================================================
# SNR（论文 §IV：谐波功率 / 其余分量功率）
# ============================================================

def snr_db(y, fs: float, f0: float, guard_bins: int = 4, warmup: int = 4096,
           window: str = "blackmanharris"):
    """信噪比 [dB]：``10*log10(P_harm / P_rest)``。

    定义（本复现自行规定，论文未给实现细节；见 README「SNR 定义」）：
    * 去掉前 ``warmup`` 个采样（AA 递推的启动瞬态）；
    * ``window`` 加窗 + rFFT（抑制非整周期截断的泄漏）；
    * 信号 = 每个谐波 ``k*f0`` 附近 ``±guard_bins`` 个 bin 的功率之和（k = 1..Nyquist/f0）；
    * 噪声 = 其余 bin 功率之和；两端各 ``guard_bins`` 个 bin（DC 与 Nyquist 附近）
      在信号与噪声中都不计，保证两部分互斥。
    """
    y = np.asarray(y, dtype=float)[warmup:]
    n = len(y)
    w = ss.get_window(window, n, fftbins=True)
    spec = np.abs(np.fft.rfft(y * w)) ** 2
    n_bins = len(spec)
    harm = np.zeros(n_bins, dtype=bool)
    k = 1
    while k * f0 < fs / 2.0:
        c = int(round(k * f0 / fs * n))
        harm[max(c - guard_bins, 0):min(c + guard_bins + 1, n_bins)] = True
        k += 1
    used = np.zeros(n_bins, dtype=bool)
    used[guard_bins:n_bins - guard_bins] = True
    signal = float(spec[harm & used].sum())
    noise = float(spec[used & ~harm].sum())
    if noise <= 0.0:
        return float("inf")
    return 10.0 * np.log10(signal / noise)


def piano_frequencies():
    """88 键钢琴基频 [Hz]（MIDI 21..108，A4 = 440）。"""
    midi = np.arange(21, 109)
    return 440.0 * 2.0 ** ((midi - 69) / 12.0)


def analytic_saw_snr(f0: float, fs: float):
    """理想锯齿波点采样后的**解析** SNR：谐波功率 / 折叠分量功率。

    幅度为 ``2/(pi*k)`` 的理想锯齿波，第 k 次谐波超过 Nyquist 后折叠。
    功率与 ``1/k^2`` 成正比，故

        SNR = 10*log10( Σ_{k<=K} 1/k² / Σ_{k>K} 1/k² ),  K = floor(fs/2/f0)

    尾和用三伽玛函数精确计算（``Σ_{k>K} 1/k² = ψ'(K+1)``）。
    本函数给出「trivial 方法 SNR 应该等于多少」的独立锚点，用于检验 SNR 口径。
    """
    K = int(np.floor((fs / 2.0) / f0))
    tail = float(sp.polygamma(1, K + 1))
    return 10.0 * np.log10((np.pi ** 2 / 6.0 - tail) / tail)
