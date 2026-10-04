"""IIR 多相版（**并行全通和**）：解析（单边）上采样 + 复整形，8× 过采样，椭圆半带 rs=100 dB。

架构与 FIR 版（`polyphase_analytic_ovs.py`，同为 L=8）相同，只把三个滤波器换成「并行全通和」：

    x[n] @48k ──┬─ 实链路：半带全通和插值 ×2 **三级级联** → 实整形 P(x)=x+x²+x³
                └─ 复链路：**解析**（单边）上采样（插值三级 ∘ 解析全通）→ 复整形 H(z)=z+z²+z³ → Re
                ── 半带全通和抗混叠 + 抽取 ÷2 **三级级联** → y[n] @48k

半带结构本身只能 2×，8× 由 **3 级级联**半带插值给出：48k → 96k → 192k → 384k
（抽取端对称地三级级联）。**三级用同一组系数**（下面那组 N=19/rs=100 dB 节参数）：

* 限制通带的只有**第一级**——半带边沿 = 该级速率的 fs/4 = 24 kHz ≥ 需要的 22.7 kHz；
  后两级的半带交叉点（48 kHz、96 kHz）远高于信号带，只负责压掉各自补零产生、会落回信号带的镜像；
* 数值验证（脚本 [1b]）：三级复合响应通带 [0, 22.67k] 最大起伏 < 1e-6 dB，
  进入第一级阻带（≥ 25.33 kHz）后最坏 −100 dB（= 设计 rs），与单级半带一致。

解析转换（半带插值 → 解析全通 `1/2[A0(−ζ)+j·z⁻¹A1(−ζ)]`）**保持在最后一级之后、抽取之前**，
位置与 2× 版一致，不改动已验证的结构。每级的两条全通链都跑在该级**输入速率**上
（级 1 = 48k、级 2 = 96k、级 3 = 192k；解析级同样在 192k）。

为什么用全通和（依据 `labs/tf2ca/` 与 `include/qwqdsp/filter/{parallel_allpass,iir_hilbert}.hpp`）：

* 奇数阶低通 H = 1/2[A0(z) + A1(z)]（两条全通链，极点按 |p| 升序逐节交替分组）；
* 半带（Wp + Ws = 1）的极点成 ±jρ 对 → 两条链都是 **zeta = z^-2** 的函数，
  即 `P(z_L)`（z_L = z^-2），**分支跑在半速率** —— 这就是 IIR 也能多相的原因
  （直接型做不到：伴随矩阵 A^L 在高阶下病态，实测误差 4.5e-1）；
* 半带通带边沿 = fs_up/4 = **输入采样率/2**，正好是要求的通带。

因果约定（本文件的实测判据，全部以 scipy 的 `freqz` 真值为准）：

* 低速率全通链 ``A(zeta) = D#(zeta)/D(zeta)``，``D = prod(zeta - p_k)``，
  ``D# = prod(1 - p_k zeta)``（= ``hilbert_odd.zeta_chain_arg`` 的写法）。
  zeta 极点 ``p_k = -1/a_k``（**负实数，模 > 1**），节参数 ``a_k = ρ²``（`hilbert_odd.py:275`）；
  写成低速率滤波器（延迟 zeta^-1）就是 ``b = poly(p)``、``a = poly(1/p)·prod(-p)``，
  落在单位圆内 → 因果稳定。半带 ``H = 1/2[A0(zeta) + z^-1 A1(zeta)]``
  与 scipy 直给的椭圆半带最大偏差 8.6e-8（= 设计的固有残差：`ellip` 的 A 末位系数 ~1.9e-7 ≠ 0）；
* 频率旋转 ``z -> -jz`` 在 zeta 上就是 ``zeta -> -zeta``，即低速率系数**奇次项取负**：
  ``A(-zeta)`` 的系数 b[k] -> (-1)^k b[k]。解析滤波器
  ``H_an = 1/2[A0(-zeta) + j z^-1 A1(-zeta)] = H_hb(-j z)``（实测偏差 8.6e-8，负频 -100 dB）。
  注意：两条支路互为 Hilbert 对（相位差处处 90°），所以 ``H_an`` 是**平坦单边**滤波器
  （正频全通、负频 0），**不是**单边低通；
* 因此解析「插值」必须是**半带插值 ∘ 单边（Hilbert）滤波**：零插值产生 9k/39k 镜像，
  平坦单边的 H_an 会把 39k 也放过去，三次整形后 2·39+11 -> 7k 混进基带（实测 -32 dB）。
  级联后 = 单边低通（正频 0..fs_in/2、负频与镜像 0）。
  L=8 时全链路每输入采样 126 MAC（实）/ 198 MAC（复，多出解析级的 72），见 [5]。

设计闭式（`labs/tf2ca/hilbert_even.py`）：自由参数只有 (N, rs)。
  eps_s = sqrt(10^(rs/10)-1), eps_p = 1/eps_s → rp = 10log10(1+eps_p^2)
  q1 = exp(-pi*K(1-m1)/K(m1)), m1 = (eps_p^2)^2, q = q1^(1/N)
  k = theta2^2/theta3^2 (q 级数，theta2 从 m=0 起), Wp = (2/pi)*atan(sqrt(k)), Ws = 1-Wp

用法: python qwqdsp/labs/adaa_iir/polyphase_allpass_iir.py
输出: qwqdsp/labs/adaa_iir/output/polyphase_allpass_iir.png
"""
from __future__ import annotations

import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sg
from scipy.special import ellipk, ellipkm1

FS_IN = 48000.0
L = 8
NSTAGE = 3                       # 半带级联级数：2^3 = 8×
FS_UP = FS_IN * L
F1, F2 = 9000.0, 11000.0
A1 = A2 = 0.25
B1, B2, B3 = 1.0, 1.0, 1.0
RS_DB = 100.0
N_HB = 19
NX = 12000
NDISC, NANA = 4800, 2400

FUND = {F1, F2}
HARM = {2 * F1, 2 * F2, 3 * F1, 3 * F2}
SUM_IMD = {F1 + F2, 2 * F1 + F2, F1 + 2 * F2}
DIFF_IMD = {abs(F2 - F1), abs(2 * F1 - F2), abs(F1 - 2 * F2)}


def color_of(f):
    if f in DIFF_IMD:
        return "tab:red"
    if f in SUM_IMD:
        return "tab:orange"
    if f in HARM:
        return "tab:blue"
    return "0.35"


def db(x):
    return 20.0 * np.log10(np.maximum(np.abs(x), 1e-300))


# ------------------------------------------------------------
# 设计（闭式，公式与 labs/tf2ca/hilbert_even.py 一致）
# ------------------------------------------------------------

def k_from_q(q: float) -> float:
    n = 1
    while q ** (n * n) > 1e-18 and n < 500:
        n += 1
    m = np.arange(0, n + 1)
    s1 = np.sum(q ** (m * (m + 1)))
    s2 = np.sum(q ** (np.arange(1, n + 1) ** 2))
    return float(4 * np.sqrt(q) * s1 ** 2 / (1 + 2 * s2) ** 2)


def ripple_rp(rs_db: float) -> float:
    es = np.sqrt(10 ** (rs_db / 10) - 1)
    return float(10 * np.log10(1 + 1 / es ** 2))


def halfband_Wp(order: int, rs_db: float):
    es = np.sqrt(10 ** (rs_db / 10) - 1)
    k1 = (1.0 / es) ** 2
    m1 = k1 ** 2
    q1 = np.exp(-np.pi * ellipkm1(m1) / ellipk(m1))
    q = q1 ** (1.0 / order)
    k = k_from_q(q)
    return float(2 / np.pi * np.arctan(np.sqrt(k))), k


def design_halfband(order=N_HB, rs_db=RS_DB):
    """半带低通 + 极点按 |p| 升序逐节交替分成两条链（返回两条链的 z 域分母）。"""
    Wp, k = halfband_Wp(order, rs_db)
    rp = ripple_rp(rs_db)
    b, a = sg.ellip(order, rp, rs_db, Wp)
    poles = sg.ellip(order, rp, rs_db, Wp, output="zpk")[1]
    used = np.zeros(len(poles), dtype=bool)
    nodes = []
    for i in np.argsort(np.abs(poles)):          # 按 |p| 升序
        if used[i]:
            continue
        p = poles[i]
        if abs(p.imag) < 1e-12:
            nodes.append(np.array([p]))
            used[i] = True
        else:
            j = int(np.argmin(np.where(used, 1e9, np.abs(poles - np.conj(p)))))
            nodes.append(np.array([p, np.conj(p)]))
            used[i] = used[j] = True
    return dict(b=b, a=a, rp=rp, Wp=Wp, k=k, poles=poles, n_nodes=len(nodes))


# ------------------------------------------------------------
# 低速率全通链（变量 zeta = z^-2，延迟 zeta^-1）
# ------------------------------------------------------------

def chain_ba(alphas):
    """节参数 a=ρ² → 低速率全通链的 (b, a)（zeta^-1 的降幂，lfilter 约定）。

    链 ``A(zeta) = prod (1 - p zeta) / (zeta - p)``，zeta 极点 ``p = -1/a``
    （半带：`hilbert_odd.py:275`，负实数、|p|>1）。以 zeta^-1 为延迟、用 lfilter 递推：
        b = poly(p)                    （分子根 p，在单位圆外 → 零点）
        a = poly(1/p)·prod(-p)         （分母根 1/p，在单位圆内 → 极点，因果稳定）
    """
    a = np.asarray(alphas, dtype=float)
    p = -1.0 / a
    return np.poly(p), np.poly(1.0 / p) * np.prod(-p)


def rot_low(b, a):
    """频率旋转 ``zeta -> -zeta``：系数奇次项取负（A(-zeta)，即解析版链路）。

    低速率 H(zeta^-1) = Σ c[k] zeta^-k，代入 zeta -> -zeta 得 Σ c[k](-1)^k zeta^-k。
    """
    b = np.asarray(b, dtype=float).copy()
    a = np.asarray(a, dtype=float).copy()
    b[1::2] *= -1.0
    a[1::2] *= -1.0
    return b, a


def chain_resp(b, a, w):
    """低速率链在高速率频率 w 上的响应：在 zeta^-1 = exp(2jw) 处求值（= freqz 在 2w）。"""
    return sg.freqz(b, a, worN=2.0 * np.asarray(w, dtype=float))[1]


# ------------------------------------------------------------
# 直接高速率参考（把并联全通和展开成单个有理式，用于 [2] 的逐点对照）
# ------------------------------------------------------------

def _z2poly(c):
    """zeta^-1 的降幂系数 -> z^-1 的降幂系数（zeta^-1 = z^-2，插到偶次项）。"""
    c = np.asarray(c, dtype=float)
    n = len(c) - 1
    out = np.zeros(2 * n + 1)
    for k, v in enumerate(c):
        out[2 * k] = v
    return out


def _add(*terms):
    """按 lfilter 约定（下标 k ↔ z^-k）对齐相加：在高次端补零。"""
    n = max(len(t) for t in terms)
    return sum(np.pad(t, (0, n - len(t))) for t in terms)


def halfband_high_rate(bL, aL, bS, aS, gain):
    """半带结构的高速率有理式：``gain · [A0(z^-2) + z^-1 A1(z^-2)]``。

    gain=1 对应插值支路（通带增益 2，补零的 ×L），gain=0.5 对应抽取支路（单位增益）。
    """
    BL, AL, BS, AS = _z2poly(bL), _z2poly(aL), _z2poly(bS), _z2poly(aS)
    num = gain * _add(np.convolve(BL, AS), np.concatenate(([0.0], np.convolve(BS, AL))))
    den = np.convolve(AL, AS)
    return num, den


def analytic_only_high_rate(bLr, aLr, bSr, aSr):
    """单边（解析全通）高速率有理式：``A0(-z^-2) + j z^-1 A1(-z^-2)``（平坦单边、不选频）。"""
    BLr, ALr, BSr, ASr = _z2poly(bLr), _z2poly(aLr), _z2poly(bSr), _z2poly(aSr)
    ne = _add(np.convolve(BLr, ASr), 1j * np.concatenate(([0.0], np.convolve(BSr, ALr))))
    de = np.convolve(ALr, ASr)
    return ne, de


# ------------------------------------------------------------
# 多相结构（各级的两条链只跑该级输入速率）
# ------------------------------------------------------------

def hb_interp(x, bL, aL, bS, aS):
    """半带全通和插值 ×2（一级）：y[2m] = A0(x)[m]、y[2m+1] = A1(x)[m]。

    通带增益 2，零插值（×2 补零）后净增益 1；级联 NSTAGE 级给出 2^NSTAGE 倍。
    """
    y = np.zeros(2 * len(x))
    y[0::2] = sg.lfilter(bL, aL, x)
    y[1::2] = sg.lfilter(bS, aS, x)
    return y


def analytic_stage(u, bLr, aLr, bSr, aSr):
    """单边（Hilbert）滤波一级：``E'(z) = A0(-z^-2) + j z^-1 A1(-z^-2)``。

    逐点展开（y[2m] / y[2m+1]，v_e/v_o = u 的偶/奇子序列）：
        y[2m]   = A0r(v_e)[m] + j·A1r(v_o)[m-1]
        y[2m+1] = A0r(v_o)[m] + j·A1r(v_e)[m]
    （z^-1 把奇/偶支路解耦后给出交叉项；`E'` 是**平坦单边**、不选频 —— 抗镜像由前面的半带负责。）
    """
    ve = u[0::2]
    vo = u[1::2]
    a1ve = sg.lfilter(bSr, aSr, ve)
    a1vo = sg.lfilter(bSr, aSr, vo)
    y = np.zeros(2 * len(ve), dtype=complex)
    y[0::2] = sg.lfilter(bLr, aLr, ve) + 1j * np.concatenate(([0.0], a1vo[:-1]))
    y[1::2] = sg.lfilter(bLr, aLr, vo) + 1j * a1ve
    return y


def dec_stage(w, bL, aL, bS, aS):
    """半带结构抽取 ÷2（一级）：y[m] = 1/2·[A0(e)[m] + A1(o)[m-1]]（e/o = 偶/奇子序列）。

    相位取 ``filtered[2m]``（无额外延迟），因此三级级联后与「同结构有理式 + [::2] 三次」
    逐点一致（见 [2]）。
    """
    ce = sg.lfilter(bL, aL, w[0::2])
    co = sg.lfilter(bS, aS, w[1::2])
    m = len(ce)
    return 0.5 * (ce[:m] + np.concatenate(([0.0], co))[:m])


def upsample_real(x, bL, aL, bS, aS):
    """实链路 8× 插值：三级级联半带插值（同一组系数），分支分别跑在 48k/96k/192k。"""
    u = x
    for _ in range(NSTAGE):
        u = hb_interp(u, bL, aL, bS, aS)
    return u


def upsample_analytic(x, bL, aL, bS, aS, bLr, aLr, bSr, aSr):
    """复链路 8× 解析插值：三级级联半带插值 ∘ 末级之后的单边滤波（合起来 = 单边低通）。"""
    return analytic_stage(upsample_real(x, bL, aL, bS, aS), bLr, aLr, bSr, aSr)


def decimate(w, bL, aL, bS, aS):
    """8× 抽取：三级级联半带抽取（同一组系数），分支分别跑在 192k/96k/48k。"""
    v = w
    for _ in range(NSTAGE):
        v = dec_stage(v, bL, aL, bS, aS)
    return v


def line_levels(seg, fs, freqs):
    w = np.hanning(len(seg))
    sp = np.abs(np.fft.rfft(seg * w)) * 2.0 / w.sum()
    return {f: float(db(sp[max(0, int(round(f / (fs / len(seg)))) - 2):
                         int(round(f / (fs / len(seg)))) + 3].max())) for f in freqs}


def spec(seg, fs):
    w = np.hanning(len(seg))
    sp = np.abs(np.fft.rfft(seg * w)) * 2.0 / w.sum()
    return np.fft.rfftfreq(len(seg), 1.0 / fs), db(sp)


# ------------------------------------------------------------
# main
# ------------------------------------------------------------

def main() -> int:
    fails = 0
    d = design_halfband()
    print(f"半带：N={N_HB}（奇），rs={RS_DB:.0f} dB，由 rs 定死 rp={ripple_rp(RS_DB):.3e} dB，"
          f"Wp={d['Wp']:.7f}（crossover 0.5）")
    print(f"  第一级（fs=96k）：通带边沿 {d['Wp']*2*FS_IN/2000:.2f} kHz，"
          f"阻带边沿 {(1-d['Wp'])*2*FS_IN/2000:.2f} kHz；半带边沿 = fs/4 = {FS_IN/2000:.1f} kHz")
    print(f"  后两级（fs=192k/384k）半带交叉点 {FS_IN/1000:.1f}/{2*FS_IN/1000:.1f} kHz，"
          f"远高于信号带，只压各自补零产生的镜像")

    # 节参数 a = ρ²（直接取设计极点，比求根稳定）：按 |p| 升序、逐对交替分给两条链
    rho2 = np.sort(np.abs(d["poles"][d["poles"].imag > 0]) ** 2)
    alphas = [rho2[0::2], rho2[1::2]]
    for i, al in enumerate(alphas):
        print(f"  链{i}：{len(al)} 节，节参数 a=ρ² = {np.round(al, 6)}")
    i0, i1 = 0, 1                                    # 延迟 z^-1 落在短链 A1 上
    print(f"  延迟支路（虚部槽）= 链{i1}（节数 {len(alphas[i1])}）")

    # 低速率全通链（zeta 极点 p = -1/a），以及旋转（zeta -> -zeta）后的解析链
    bL, aL = chain_ba(alphas[i0])
    bS, aS = chain_ba(alphas[i1])
    bLr, aLr = rot_low(bL, aL)
    bSr, aSr = rot_low(bS, aS)

    # ---------- 1. 全通和复原半带 / 解析单边性 ----------
    w = np.linspace(-np.pi, np.pi, 32769)
    Hhb = 0.5 * (chain_resp(bL, aL, w) + np.exp(-1j * w) * chain_resp(bS, aS, w))
    Href = sg.freqz(d["b"], d["a"], worN=w)[1]
    err_hb = float(np.max(np.abs(Hhb - Href)))
    Han = 0.5 * (chain_resp(bLr, aLr, w) + 1j * np.exp(-1j * w) * chain_resp(bSr, aSr, w))
    Href_an = sg.freqz(d["b"], d["a"], worN=w - np.pi / 2)[1]     # H_hb(-j z)
    err_an = float(np.max(np.abs(Han - Href_an)))
    guard = np.pi * (1 - 2 * d["Wp"])
    core = (np.abs(w) > guard) & (np.abs(w) < np.pi - guard)
    pos, neg = core & (w > 0), core & (w < 0)
    neg_db = db(np.max(np.abs(Han[neg])))
    print(f"\n[1] 全通和复原半带 1/2[A0(zeta)+z^-1 A1(zeta)]：与 scipy 直给 H(z) 最大偏差 "
          f"{err_hb:.2e}（= ellip 的固有残差，A 末位系数 ~2e-7）")
    print(f"    解析滤波器 1/2[A0(-zeta)+j z^-1 A1(-zeta)] = H_hb(-jz)：最大偏差 {err_an:.2e}")
    print(f"    解析滤波器（正频 core）：|H| ∈ [{np.min(np.abs(Han[pos])):.6f}, "
          f"{np.max(np.abs(Han[pos])):.6f}]（理想 1，平坦单边）")
    print(f"    负频抑制 {neg_db:.1f} dB（设计 rs = −{RS_DB:.0f} dB）；"
          f"解析级（@384k）过渡带宽度 {(1-2*d['Wp'])*FS_UP/2000:.2f} kHz"
          f"（单边性在 |f| ≳ {(1-2*d['Wp'])*FS_UP/4000:.2f} kHz 外成立）")
    fails += err_hb > 1e-6 or err_an > 1e-6 or neg_db > -(RS_DB - 5)

    # ---------- 1b. 三级级联复合响应（数值验证「同一组系数可用」） ----------
    fk2 = np.linspace(1.0, FS_UP / 2 - 1, 400001)
    comp = np.ones_like(fk2, dtype=complex)
    for k in range(NSTAGE):
        rate = FS_UP / 2 ** k
        wk = 2 * np.pi * fk2 / rate
        comp *= 0.5 * (chain_resp(bL, aL, wk) + np.exp(-1j * wk) * chain_resp(bS, aS, wk))
    fpb = 0.45 * FS_IN                                   # 信号带上界 21.6 kHz（同 FIR 版）
    stop_lo = (1 - d["Wp"]) * 2 * FS_IN / 2              # 第一级阻带边沿 25.31 kHz
    pb_dev = db(np.max(np.abs(np.abs(comp[fk2 <= fpb]) - 1)))
    sb_worst = db(np.max(np.abs(comp[fk2 >= stop_lo])))
    trans = db(np.max(np.abs(comp[(fk2 > 2 * d["Wp"] * FS_IN / 2) & (fk2 < stop_lo)])))
    print(f"[1b] 三级级联复合响应（同一组系数）：通带 [0, {fpb/1000:.2f}k] 最大起伏 {pb_dev:.1f} dB；"
          f"阻带 [≥ {stop_lo/1000:.2f}k] 最坏 {sb_worst:.1f} dB")
    print(f"     过渡带 [{2*d['Wp']*FS_IN/2000:.2f}k, {stop_lo/1000:.2f}k] 最坏 {trans:.1f} dB"
          f"（半带固有过渡带，与单级相同；信号带 21.6k 以内不受影响）")
    fails += pb_dev > -80.0 or sb_worst > -95.0

    # ---------- 2. 多相 vs 直接（逐级有理式） ----------
    rng = np.random.default_rng(1)
    xt = rng.standard_normal(4000)
    nh, dh = halfband_high_rate(bL, aL, bS, aS, gain=1.0)
    nh1, dh1 = halfband_high_rate(bL, aL, bS, aS, gain=0.5)
    # 解析级的有理式（在 384k 上分支仍跑半速率）
    ne, de = analytic_only_high_rate(bLr, aLr, bSr, aSr)

    def naive_up_real(z):
        v = z
        for _ in range(NSTAGE):
            zv = np.zeros(2 * len(v))
            zv[::2] = v
            v = sg.lfilter(nh, dh, zv)
        return v

    def naive_dec(v):
        for _ in range(NSTAGE):
            v = sg.lfilter(nh1, dh1, v)[::2]
        return v

    ur = upsample_real(xt, bL, aL, bS, aS)
    e_up = float(np.max(np.abs(ur - naive_up_real(xt))))
    ua = analytic_stage(ur, bLr, aLr, bSr, aSr)
    e_an = float(np.max(np.abs(ua - sg.lfilter(ne, de, ur))))
    wu = rng.standard_normal(8 * len(xt))
    dd, dn = decimate(wu, bL, aL, bS, aS), naive_dec(wu)
    m = min(len(dd), len(dn))
    e_dec = float(np.max(np.abs(dd[:m] - dn[:m])))
    print(f"\n[2] 三级级联 vs 逐级「补零 + 同结构有理式」")
    print(f"    实插值最大偏差 {e_up:.2e}；解析级 {e_an:.2e}；抽取 {e_dec:.2e}")
    print(f"    （各级补零把通带降到 1/2，故插值增益 1 → 通带 2；抽取增益 0.5）")
    fails += e_up > 1e-9 or e_an > 1e-9 or e_dec > 1e-9

    # ---------- 3. 链路 ----------
    n = np.arange(NX)
    x = A1 * np.cos(2 * np.pi * F1 * n / FS_IN) + A2 * np.cos(2 * np.pi * F2 * n / FS_IN)
    zu = upsample_analytic(x, bL, aL, bS, aS, bLr, aLr, bSr, aSr)
    xr = upsample_real(x, bL, aL, bS, aS)
    y_c = decimate(np.real(B1 * zu + B2 * zu ** 2 + B3 * zu ** 3), bL, aL, bS, aS)
    y_r = decimate(B1 * xr + B2 * xr ** 2 + B3 * xr ** 3, bL, aL, bS, aS)

    seg = lambda y: y[NDISC:NDISC + NANA]
    freqs = sorted(f for f in (FUND | HARM | SUM_IMD | DIFF_IMD) if f < FS_IN / 2)
    lv_c, lv_r = line_levels(seg(y_c), FS_IN, freqs), line_levels(seg(y_r), FS_IN, freqs)
    ref_db = db(A1)

    def kind(f):
        return ("基波" if f in FUND else "自身谐波" if f in HARM
                else "求和互调" if f in SUM_IMD else "差频互调")

    print(f"\n[3] 输出谱（相对输入音幅度 dB；'—' = 未出现）")
    print(f"    {'频率':>7} {'归属':<9} {'实链路':>10} {'复链路':>10}")
    for f in freqs:
        rr, rc = lv_r[f] - ref_db, lv_c[f] - ref_db
        print(f"    {f:7.0f} {kind(f):<9} " +
              " ".join(f"{v:>10.1f}" if v > -140 else f"{'—':>10}" for v in (rr, rc)))

    diff_f = sorted(DIFF_IMD)
    d_r = [lv_r[f] - ref_db for f in diff_f]
    d_c = max(lv_c[f] - ref_db for f in DIFF_IMD)
    fund_r = max(lv_r[f] - ref_db for f in FUND)
    fund_c = max(lv_c[f] - ref_db for f in FUND)
    print(f"\n[4] 基波：实链路 {fund_r:+.1f} dB（基波 + 三阶贡献），复链路 {fund_c:+.1f} dB（单边 ⇒ 无三阶贡献）")
    print(f"    差频互调 2k/7k/13k：实链路 " +
          "/".join(f"{v:+.1f}" for v in d_r) + f" dB；复链路最强 {d_c:+.1f} dB（椭圆泄漏地板）")
    fails += (abs(fund_r - 1.14) > 1.0 or abs(fund_c) > 0.6
              or abs(d_r[0] + 12.0) > 1.5 or abs(d_r[1] + 26.6) > 1.5 or abs(d_r[2] + 26.6) > 1.5
              or d_c > -95.0)

    # ---------- 4. 代价 ----------
    n_sec = len(alphas[i0]) + len(alphas[i1])            # 每条链合计 9 个一阶节 = 一趟
    ticks = [2 ** k for k in range(NSTAGE)]              # 各级每输入采样 tick 数 1/2/4
    mac_up_r = n_sec * sum(ticks)                        # 三级级联插值
    mac_an_extra = 2 * n_sec * ticks[-1]                 # 解析级：2 条链 × 4 tick
    mac_dec = n_sec * sum(ticks)                         # 三级级联抽取
    mac_r = mac_up_r + mac_dec
    mac_an = mac_up_r + mac_an_extra + mac_dec
    t0 = time.perf_counter()
    for _ in range(5):
        upsample_analytic(x, bL, aL, bS, aS, bLr, aLr, bSr, aSr)
    t_up = (time.perf_counter() - t0) / 5 / len(x) * 1e6
    print(f"\n[5] 每输入采样 MAC（{n_sec} 个一阶节/趟；三级 tick 数 1/2/4）")
    print(f"    实链路：插值 {n_sec}×(1+2+4) = {mac_up_r} ＋抽取 {mac_dec} = {mac_r}")
    print(f"    复链路：插值 {mac_up_r}（三级半带）＋解析级 {mac_an_extra}"
          f"（2 链×4 tick）＋抽取 {mac_dec} = {mac_an}")
    print(f"    对照 FIR 多相 1536 → 省 {1536/mac_an:.1f}×；直接型椭圆 280 → 省 {280/mac_an:.2f}×；"
          f"实测复链路插值 {t_up:.3f} µs/输入采样")

    # ---------- 图 ----------
    fig = plt.figure(figsize=(13.2, 8.8))
    gs = fig.add_gridspec(2, 3, hspace=0.46, wspace=0.3)

    ax = fig.add_subplot(gs[0, 0])
    fk = np.linspace(-48e3 + 1, 48e3 - 1, 8193)
    comp_plot = np.ones_like(fk, dtype=complex)
    for k in range(NSTAGE):
        rate = FS_UP / 2 ** k
        wkk = 2 * np.pi * fk / rate
        comp_plot *= 0.5 * (chain_resp(bL, aL, wkk) + np.exp(-1j * wkk) * chain_resp(bS, aS, wkk))
    w384 = 2 * np.pi * fk / FS_UP
    Han_plot = 0.5 * (chain_resp(bLr, aLr, w384) + 1j * np.exp(-1j * w384) * chain_resp(bSr, aSr, w384))
    ax.plot(fk / 1000, db(np.abs(comp_plot)), lw=1.2, color="0.5", label="三级复合实响应 |H|")
    ax.plot(fk / 1000, db(np.abs(comp_plot * Han_plot)), lw=1.4, color="tab:blue",
            label="三级复合解析 |H·H_an|（单边）")
    for s in (-1, 1):
        ax.axvline(s * d["Wp"] * 2 * FS_IN / 2000, color="tab:red", ls=":", lw=1.2)
    ax.set_xlim(-48, 48)
    ax.set_ylim(-130, 8)
    ax.set_xlabel("频率 [kHz] @ 384 kHz")
    ax.set_ylabel("幅度 [dB]")
    ax.set_title(f"(a) 三级级联复合响应（N={N_HB}/rs={RS_DB:.0f} dB，同组系数）：解析负频 −100 dB", fontsize=9.5)
    ax.legend(fontsize=8, loc="lower center")
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 1])
    fz = np.linspace(0.02, 0.98, 2000) * np.pi
    ax.plot(fz / np.pi, np.abs(chain_resp(bLr, aLr, fz)), lw=1.1, label=f"链{i0}（实部槽）")
    ax.plot(fz / np.pi, np.abs(chain_resp(bSr, aSr, fz)), lw=1.1, label=f"链{i1}（虚部槽）")
    ax.set_ylim(0.9, 1.1)
    ax.set_xlabel("归一化频率 (×π)")
    ax.set_ylabel("|A|")
    ax.set_title("(b) 两条链各自都是全通（|A|≡1）：解析性由相位差给出", fontsize=9.5)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 2])
    dphi = np.unwrap(np.angle(chain_resp(bLr, aLr, fz)) -
                     np.angle(np.exp(-1j * fz) * chain_resp(bSr, aSr, fz)))
    ax.plot(fz / np.pi, np.abs(np.abs(dphi) * 180 / np.pi - 90), lw=1.2)
    ax.set_xlabel("归一化频率 (×π)")
    ax.set_ylabel("偏离 90° [deg]")
    ax.set_title("(c) 两链相位差：core 段近 90°（两端守护带由 crossover 决定）", fontsize=9.5)
    ax.grid(alpha=0.3)

    def plot_out(ax, y, lv, title):
        f_, s_ = spec(seg(y), FS_IN)
        m = f_ < 25000
        ax.plot(f_[m] / 1000, s_[m] - ref_db, lw=0.7, color="0.6")
        for f in freqs:
            v = lv[f] - ref_db
            if v > -140:
                ax.vlines(f / 1000, -140, v, color=color_of(f), lw=1.5)
        ax.set_xlim(0, 25)
        ax.set_ylim(-140, 6)
        ax.set_xlabel("频率 [kHz] @ 48 kHz")
        ax.set_ylabel("相对输入音 [dB]")
        ax.set_title(title, fontsize=9.5)
        ax.grid(alpha=0.3)

    plot_out(fig.add_subplot(gs[1, 0]), y_r, lv_r, "(d) 实链路（半带全通和）：差频互调 2k/7k/13k")
    plot_out(fig.add_subplot(gs[1, 1]), y_c, lv_c, "(e) 复链路（解析全通和 + 复整形）：无差频互调")

    ax = fig.add_subplot(gs[1, 2])
    names = ["FIR 多相", "直接型椭圆", "全通和多相"]
    vals = [1536, 280, mac_an]
    ax.bar(names, vals, color=["0.5", "0.7", "tab:blue"])
    for i, v in enumerate(vals):
        ax.text(i, v * 1.03, str(v), ha="center", fontsize=9)
    ax.set_yscale("log")
    ax.set_ylabel("MAC / 输入采样")
    ax.set_title(f"(f) 复链路代价：{vals[0]/vals[2]:.0f}× 少于 FIR、{vals[1]/vals[2]:.1f}× 少于直接型", fontsize=9.5)
    ax.grid(alpha=0.3, axis="y")

    fig.suptitle(f"IIR 多相（并行全通和，椭圆半带 N={N_HB}/rs={RS_DB:.0f} dB，L={L}）：解析上采样 + 复整形",
                 fontsize=11)
    out = Path(__file__).parent / "output" / "polyphase_allpass_iir.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=140)
    print(f"\n已写出 {out}")
    print(f"断言失败项：{fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
