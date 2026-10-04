"""IIR 版：解析（单边）上采样 + 复整形，全部用 **−100 dB 阻带椭圆滤波器**。

与 `polyphase_analytic_ovs.py`（FIR 版）同架构、同频率计划、同测试信号，只换滤波器技术：

    解析上采样：椭圆低通原型（通带 [0, f_p/2]、阻带自 fs_in−1.5·f_p，r_s = 100 dB）
                → **系数按 e^{jθk} 旋转**（θ = 2π·(f_p/2)/fs_up）→ 单边滤波器
                → 乘 2L（补零插值要 L 倍通带增益；解析侧再 ×2 使 Re(z) 对上带限插值）
    实数链路  ：椭圆低通（通带 [0, f_p]、阻带自 fs_in−f_p、r_s = 100 dB）×L
    抽取滤波器：椭圆低通（同规格、增益 1），递推跑在过采样率上再 ÷L

关于「IIR 多相」（本脚本给出实测结论，见 [1]）：

* **代数上**，IIR 可以做「状态空间块形式」：每 L 点更新一次状态
  `x[k+L] = A^L·x[k] + Σ_m A^m B·u[kL+L−1−m]`，`y[kL] = C·x[kL] + D·u[kL]`。
  它只对**抽取端**有意义（丢掉 L−1 个输出）；上采样端每个输出都要，直接递推已是最优。
* **数值上**，直接用高阶直接型的伴随矩阵 A^L 是灾难：12 阶椭圆的 a 系数最大 574、
  `‖A^8‖₂ = 1.2e7`，误差放大到 1e22（见 [1] 的实测）。所以工程上 IIR 超采样
  **就是**在过采样率上跑 SOS 递推（`lfilter`/sosfilt），省算力靠的是低阶（280 vs 1536 MAC/输入），
  不是靠多相分解。

用法: python qwqdsp/labs/adaa_iir/polyphase_analytic_iir.py
输出: qwqdsp/labs/adaa_iir/output/polyphase_analytic_iir.png
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
import scipy.signal as ss

# ------------------------------------------------------------
# 参数（与 FIR 版完全一致）
# ------------------------------------------------------------
FS_IN = 48000.0
L = 8
FS_UP = FS_IN * L
F1, F2 = 9000.0, 11000.0
A1 = A2 = 0.25
B1, B2, B3 = 1.0, 1.0, 1.0
FP = 0.45 * FS_IN
TAPS = 512
RIPPLE_DB = 0.01
STOP_DB = 100.0
NX = 12000
NDISC, NANA = 4800, 2400

FUND = {F1, F2}
HARM = {2 * F1, 2 * F2, 3 * F1, 3 * F2}
SUM_IMD = {F1 + F2, 2 * F1 + F2, F1 + 2 * F2}
DIFF_IMD = {abs(F2 - F1), abs(2 * F1 - F2), abs(F1 - 2 * F2)}
ALIAS = {FS_IN - 3 * F1, FS_IN - (2 * F1 + F2), FS_IN - (F1 + 2 * F2), FS_IN - 3 * F2}


def color_of(f):
    if f in DIFF_IMD:
        return "tab:red"
    if f in ALIAS:
        return "tab:purple"
    if f in SUM_IMD:
        return "tab:orange"
    if f in HARM:
        return "tab:blue"
    return "0.35"


def db(x):
    return 20.0 * np.log10(np.maximum(np.abs(x), 1e-300))


def freqresp(b, a, f, fs):
    """任意（含负）频率处的频率响应：H = B(e^{jω})/A(e^{jω})。"""
    f = np.atleast_1d(np.asarray(f, dtype=float))
    z = np.exp(-2j * np.pi * np.outer(f, np.arange(max(len(a), len(b)))) / fs)
    return (z[:, :len(b)] @ np.asarray(b)) / (z[:, :len(a)] @ np.asarray(a))


def group_delay(b, a, f, fs):
    """数值群延迟 [采样]（等间距 f 上对 unwrap 后的相位求导）。"""
    f = np.asarray(f, dtype=float)
    ph = np.unwrap(np.angle(freqresp(b, a, f, fs)))
    return -np.gradient(ph, 2 * np.pi * (f[1] - f[0]) / fs)


# ------------------------------------------------------------
# 设计
# ------------------------------------------------------------

def design_iir():
    n_a = int(ss.ellipord(FP / 2, FS_IN - 1.5 * FP, RIPPLE_DB, STOP_DB, fs=FS_UP)[0])
    b, a = ss.ellip(n_a, RIPPLE_DB, STOP_DB, FP / 2, btype="low", fs=FS_UP)
    th = 2 * np.pi * (FP / 2) / FS_UP
    k = np.arange(len(b))
    ba, aa = 2.0 * L * b * np.exp(1j * th * k), a * np.exp(1j * th * k)   # 单边化 + 2L
    n_r = int(ss.ellipord(FP, FS_IN - FP, RIPPLE_DB, STOP_DB, fs=FS_UP)[0])
    br, ar = ss.ellip(n_r, RIPPLE_DB, STOP_DB, FP, btype="low", fs=FS_UP)
    # 插值支路要 L 倍通带增益（补零），抽取支路保持单位增益 —— 两个不同的滤波器，别复用
    br_i = L * br
    return n_a, (ba, aa), (b, a), n_r, (br_i, ar), (br, ar)


def design_fir():
    h_lp = ss.remez(TAPS, [0, FP / 2, FS_IN - 1.5 * FP, FS_UP / 2], [1, 0], fs=FS_UP)
    k = np.arange(TAPS)
    tau = (TAPS - 1) / 2.0
    ha = 2.0 * L * h_lp * np.exp(1j * 2 * np.pi * (FP / 2) * (k - tau) / FS_UP)
    h_r = L * ss.remez(TAPS, [0, FP, FS_IN - FP, FS_UP / 2], [1, 0], fs=FS_UP)
    h_d = ss.remez(TAPS, [0, FP, FS_IN - FP, FS_UP / 2], [1, 0], fs=FS_UP)
    return ha, h_r, h_d


def polyphase_upsample(x, h, l):
    assert len(h) % l == 0
    return np.stack([np.convolve(x, h[p::l]) for p in range(l)], axis=1).ravel()


def decimate_fir(y, h, l):
    return np.convolve(y, h)[::l]


def decimate_iir(y, b, a, l):
    return ss.lfilter(b, a, y)[::l]


def ss_state_space(a, b):
    """可控标准型：x⁺ = A x + B u，y = C x + D u。"""
    n = len(a) - 1
    A = np.zeros((n, n))
    A[0, :] = -a[1:]
    if n > 1:
        A[1:, :-1] = np.eye(n - 1)
    B = np.zeros(n)
    B[0] = 1.0
    return A, B, b[1:] - b[0] * a[1:], b[0]


def ss_decimate(u, b, a, l):
    """状态空间块形式（每 L 点一次状态更新）——**仅用于本脚本 [1] 的等价性/病态演示**。"""
    A, B, C, D = ss_state_space(a, b)
    n = len(B)
    nb = len(u) // l
    U = np.asarray(u[:nb * l]).reshape(nb, l)
    G = np.empty((n, l))
    v = B.copy()
    for m in range(l):
        G[:, m] = v
        v = A @ v
    V = np.einsum("nm,km->kn", G, U[:, ::-1])
    AL = np.linalg.matrix_power(A, l)
    X = np.zeros(n)
    y = np.empty(nb)
    for k in range(nb):
        y[k] = C @ X + D * U[k, 0]
        X = AL @ X + V[k]
    return y


def line_levels(seg, fs, freqs):
    w = np.hanning(len(seg))
    sp = np.abs(np.fft.rfft(seg * w)) * 2.0 / w.sum()
    return {f: float(db(sp[max(0, int(round(f / (fs / len(seg)))) - 2):
                         int(round(f / (fs / len(seg)))) + 3].max())) for f in freqs}


def spec(seg, fs):
    w = np.hanning(len(seg))
    sp = np.abs(np.fft.rfft(seg * w)) * 2.0 / w.sum()
    return np.fft.rfftfreq(len(seg), 1.0 / fs), db(sp)


def main() -> int:
    fails = 0
    n_a, (ba, aa), (bp, ap), n_r, (br_i, ar), (br_d, ar_d) = design_iir()
    ha_f, hr_f, hd_f = design_fir()

    print(f"椭圆设计（fs_up = {FS_UP/1000:.0f} kHz，通带纹波 {RIPPLE_DB} dB，阻带 {STOP_DB:.0f} dB）")
    print(f"  解析原型：{n_a} 阶（ellipord 最小阶），极点最大模 {np.abs(np.roots(aa)).max():.4f}")
    print(f"  实数插值 / 抽取：{n_r} 阶；FIR 对照 {TAPS} taps")
    n_small = int(ss.ellipord(FP / 2, FS_IN - 1.5 * FP, 0.001, STOP_DB, fs=FS_UP)[0])
    print(f"  （通带纹波降到 0.001 dB 只需 {n_small} 阶：椭圆阶数对纹波几乎不敏感）")

    # ---------- 1. 规格 ----------
    fchk = np.concatenate([np.linspace(-FS_UP / 2, FS_UP / 2, 1 << 16)])
    Hi, Hr, Hp = (freqresp(ba, aa, fchk, FS_UP), freqresp(br_i, ar, fchk, FS_UP),
                  freqresp(bp, ap, fchk, FS_UP))
    pb = (fchk >= 0) & (fchk < FP)
    img = (fchk > FS_IN - FP) & (fchk < FS_IN + FP)
    neg = fchk < -2 * FP
    print(f"\n[0] 实测：椭圆原型 通带 {db(np.abs(Hp[np.abs(fchk) < FP/2]).min()):.3f}…0 dB，"
          f"阻带 {db(np.max(np.abs(Hp[np.abs(fchk) >= FS_IN - 1.5 * FP]))):.1f} dB")
    print(f"    解析滤波器 正通带偏差 {db(np.max(np.abs(np.abs(Hi[pb]) / (2 * L) - 1))):.1f} dB，"
          f"镜像带 {db(np.max(np.abs(Hi[img])) / (2 * L)):.1f} dB，"
          f"负频率侧 {db(np.max(np.abs(Hi[neg])) / (2 * L)):.1f} dB")
    fails += db(np.max(np.abs(Hi[img])) / (2 * L)) > -95.0
    fails += db(np.max(np.abs(Hi[neg])) / (2 * L)) > -95.0

    # ---------- 2. 「状态空间块形式」的等价性与病态 ----------
    rng = np.random.default_rng(0)
    u = rng.standard_normal(8000)
    b2, a2 = ss.butter(2, 0.2)
    e_small = float(np.max(np.abs(ss_decimate(u, b2, a2, L) - ss.lfilter(b2, a2, u)[::L])))
    A, B, C, D = ss_state_space(ap, bp)
    nrm = float(np.linalg.norm(np.linalg.matrix_power(A, L), 2))
    e_big = float(np.max(np.abs(ss_decimate(u, bp, ap, L) - ss.lfilter(bp, ap, u)[::L])))
    print(f"\n[1] IIR「多相」（状态空间块形式，每 L 点一次状态更新）")
    print(f"    代数等价：2 阶 Butterworth 上最大偏差 {e_small:.2e} ✓（公式没错）")
    print(f"    数值可行性：{n_r} 阶椭圆 a 系数最大 {np.abs(ap).max():.0f}、‖A^{L}‖₂ = {nrm:.1e}"
          f" → 同输入下偏差 {e_big:.1e} ✗（输出幅值量级 {np.max(np.abs(ss.lfilter(bp, ap, u))):.2f}）")
    print(f"    结论：块形式只对**抽取端**有意义，且高阶直接型伴随矩阵病态 → 工程上用 SOS 在过采样率上递推")
    fails += e_small > 1e-12

    # ---------- 3. 两条链路 ----------
    n = np.arange(NX)
    x = A1 * np.cos(2 * np.pi * F1 * n / FS_IN) + A2 * np.cos(2 * np.pi * F2 * n / FS_IN)
    xu = np.zeros(NX * L)
    xu[::L] = x

    z_iir = ss.lfilter(ba, aa, xu)
    z_fir = polyphase_upsample(x, ha_f, L)
    z_ri = ss.lfilter(br_i, ar, xu)
    z_rf = polyphase_upsample(x, hr_f, L)

    shape = lambda z: (B1 * z + B2 * z ** 2 + B3 * z ** 3).real
    y_ci = decimate_iir(shape(z_iir), br_d, ar_d, L)
    y_ri = decimate_iir(shape(z_ri), br_d, ar_d, L)
    y_cf = decimate_fir(shape(z_fir), hd_f, L)
    y_rf = decimate_fir(shape(z_rf), hd_f, L)

    # ---------- 4. 解析上采样质量 ----------
    t = np.arange(len(z_iir)) / FS_UP
    sl = slice(4000, 20000)
    best = min(((dd, float(np.max(np.abs(z_iir[sl] - (
        A1 * np.exp(2j * np.pi * F1 * (t[sl] - dd / FS_UP))
        + A2 * np.exp(2j * np.pi * F2 * (t[sl] - dd / FS_UP)))))))
        for dd in np.arange(-20, 80, 0.25)), key=lambda p: p[1])
    gd_i = group_delay(ba, aa, np.array([F1, F2]), FS_UP)
    gd_r = group_delay(br_i, ar, np.array([F1, F2]), FS_UP)
    dev_i = db(np.abs(freqresp(ba, aa, np.array([F1, F2]), FS_UP)) / (2 * L))
    print(f"\n[2] 解析上采样（IIR）质量")
    print(f"    群延迟 @9k/11k：解析支路 {gd_i[0]:.2f} / {gd_i[1]:.2f} 采样，"
          f"实数支路 {gd_r[0]:.2f} / {gd_r[1]:.2f}（FIR 版两支路都恒为 255.5，严格对齐）")
    print(f"    两音处幅度偏差 {dev_i[0]:+.2f} / {dev_i[1]:+.2f} dB")
    print(f"    与理想解析信号的最佳对齐误差 {db(best[1] / A1):.1f} dB（延迟 {best[0]:.2f}）"
          f"—— 椭圆非线性相位 + 频移 → Re(z) 不等于带限插值，但**只影响波形**")
    per = int(round(FS_UP / 1000.0))
    Zs = np.fft.fftshift(np.fft.fft(z_iir[4000:4000 + per * 32]))
    fzs = np.fft.fftshift(np.fft.fftfreq(per * 32, 1.0 / FS_UP))
    supp = db(np.max(np.abs(Zs[fzs < -2500])) / np.max(np.abs(Zs)))
    print(f"    z 的负频率侧最强分量 {supp:.1f} dB（整周期段无泄漏；FIR 版 −105.4 dB）")
    fails += supp > -95.0

    # ---------- 5. 输出谱 ----------
    seg = lambda y: y[NDISC:NDISC + NANA]
    freqs = sorted(f for f in (FUND | HARM | SUM_IMD | DIFF_IMD | ALIAS) if f < FS_IN / 2)
    chains = (("实IIR", y_ri), ("复IIR", y_ci), ("实FIR", y_rf), ("复FIR", y_cf))
    lv = {tag: line_levels(seg(y), FS_IN, freqs) for tag, y in chains}
    ref_db = db(A1)

    def kind(f):
        return ("基波" if f in FUND else "自身谐波" if f in HARM
                else "求和互调" if f in SUM_IMD else "差频互调" if f in DIFF_IMD else "混叠")

    print(f"\n[3] 输出谱（相对输入音幅度 dB；'—' = < -140）")
    print(f"    {'频率':>7} {'归属':<9} {'实IIR':>8} {'复IIR':>8} {'实FIR':>8} {'复FIR':>8}")
    for f in freqs:
        row = [f"{v:+.1f}" if (v := lv[tag][f] - ref_db) > -140 else "—" for tag, _ in chains]
        print(f"    {f:7.0f} {kind(f):<9} " + " ".join(f"{v:>8}" for v in row))

    d_ri = max(lv["实IIR"][f] - ref_db for f in DIFF_IMD)
    d_ci = max(lv["复IIR"][f] - ref_db for f in DIFF_IMD)
    a_ci = max(lv["复IIR"][f] - ref_db for f in ALIAS)
    print(f"\n[4] 差频互调最强：实 IIR {d_ri:+.1f} dB vs 复 IIR {d_ci:+.1f} dB（抑制 {d_ri - d_ci:.1f} dB）")
    print(f"    复 IIR 混叠最强 {a_ci:+.1f} dB（FIR 版 −116.6 dB）")
    fails += d_ci > -80.0 or a_ci > -80.0

    # ---------- 6. 代价 ----------
    macs_fir = 2 * TAPS + TAPS
    macs_iir = 2 * (n_a + 1) * L + (n_r + 1) * L
    nrep = 5
    t0 = time.perf_counter()
    for _ in range(nrep):
        polyphase_upsample(x, ha_f, L)
    t_fir_up = (time.perf_counter() - t0) / nrep / len(x) * 1e6
    t0 = time.perf_counter()
    for _ in range(nrep):
        ss.lfilter(ba, aa, xu)
    t_iir_up = (time.perf_counter() - t0) / nrep / len(x) * 1e6
    print(f"\n[5] 每输入采样 MAC：FIR 复链路 {macs_fir}  vs  IIR 复链路 {macs_iir}"
          f"（{macs_fir / macs_iir:.1f}× 更省：2 支路 × {n_a+1} × L 插值 + {n_r+1} × L 抽取）")
    print(f"    实测上采样 {t_fir_up:.2f} µs vs {t_iir_up:.2f} µs / 输入采样"
          f"（多相用的是向量化 convolve、IIR 是 lfilter 的标量 C 循环，别当严格公平比较）")

    # ---------- 图 ----------
    fig = plt.figure(figsize=(13.2, 12.2))
    gs = fig.add_gridspec(3, 2, hspace=0.45, wspace=0.22)

    ax = fig.add_subplot(gs[0, 0])
    fplot = np.linspace(-FS_UP / 2, FS_UP / 2, 4097)
    ax.plot(fplot / 1000, db(freqresp(bp, ap, fplot, FS_UP)), lw=0.9, alpha=0.5, ls="--",
            color="tab:green", label="椭圆原型 |H_lp|（双边）")
    ax.plot(fplot / 1000, db(freqresp(ha_f, [1.0], fplot, FS_UP)) - db(2 * L), lw=1.0,
            alpha=0.5, color="0.5", label="FIR 解析 |H_a|/(2L)")
    ax.plot(fplot / 1000, db(freqresp(ba, aa, fplot, FS_UP)) - db(2 * L), lw=1.8,
            color="tab:blue", label=f"IIR 解析 |H_a|/(2L)（单边）")
    for s in (-1, 1):
        ax.axvline(s * FP / 1000, color="tab:red", ls=":", lw=1.2)
    ax.set_xlim(-FS_UP / 2000, FS_UP / 2000)
    ax.set_ylim(-130, 8)
    ax.set_xlabel("频率 [kHz] @ 384 kHz")
    ax.set_ylabel("幅度 [dB]")
    ax.set_title(f"(a) 椭圆 {n_a} 阶（阻带 {STOP_DB:.0f} dB）vs FIR {TAPS} taps：过渡带更陡、阶数低",
                 fontsize=9.5)
    ax.legend(fontsize=7.5, loc="lower center", framealpha=0.95)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 1])
    fg = np.linspace(200, FS_UP / 2 - 100, 4096)
    ax.plot(fg / 1000, group_delay(ba, aa, fg, FS_UP), lw=1.2, label=f"IIR 解析（{n_a} 阶）")
    ax.plot(fg / 1000, group_delay(br_i, ar, fg, FS_UP), lw=1.0, alpha=0.75, label=f"IIR 实数（{n_r} 阶）")
    ax.axhline((TAPS - 1) / 2, color="0.4", ls="--", lw=1.0, label=f"FIR（恒为 {(TAPS-1)/2:.1f}）")
    ax.axvline(FP / 1000, color="tab:red", ls=":", lw=1.2)
    ax.set_xlim(0, FS_UP / 2000)
    ax.set_ylim(0, 60)
    ax.set_xlabel("频率 [kHz] @ 384 kHz")
    ax.set_ylabel("群延迟 [采样]")
    ax.set_title("(b) 群延迟：IIR 非线性 + 通带纹波 ⇒ 两支路对不齐（FIR 严格线性相位）", fontsize=9.5)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    tagof = {id(y_ri): "实IIR", id(y_ci): "复IIR", id(y_rf): "实FIR", id(y_cf): "复FIR"}

    def plot_out(ax, y, title):
        f_, s_ = spec(seg(y), FS_IN)
        m = f_ < 25000
        ax.plot(f_[m] / 1000, s_[m] - ref_db, lw=0.7, color="0.6")
        for f in freqs:
            v = lv[tagof[id(y)]][f] - ref_db
            if v > -140:
                ax.vlines(f / 1000, -140, v, color=color_of(f), lw=1.5)
        ax.axvspan(FS_IN / 2000, 25, color="tab:red", alpha=0.06)
        ax.set_xlim(0, 25)
        ax.set_ylim(-140, 6)
        ax.set_xlabel("频率 [kHz] @ 48 kHz")
        ax.set_ylabel("相对输入音 [dB]")
        ax.set_title(title, fontsize=9.5)
        ax.grid(alpha=0.3)

    plot_out(fig.add_subplot(gs[1, 0]), y_ri, "(c) 实链路 IIR @L=8：差频互调 2k/7k/13k 照旧存在")
    plot_out(fig.add_subplot(gs[1, 1]), y_ci, "(d) 复链路 IIR @L=8：无混叠、无差频互调")

    ax = fig.add_subplot(gs[2, 0])
    for y, lab, c in ((y_ci, "复链路 IIR（椭圆 −100 dB）", "tab:green"),
                      (y_cf, "复链路 FIR（512 taps）", "0.45")):
        f_, s_ = spec(seg(y), FS_IN)
        m = f_ < 25000
        ax.plot(f_[m] / 1000, s_[m] - ref_db, lw=0.8, label=lab, color=c, alpha=0.9)
    for f in sorted(DIFF_IMD):
        ax.axvline(f / 1000, color="tab:red", ls=":", lw=1.0)
    ax.set_xlim(0, 25)
    ax.set_ylim(-120, 6)
    ax.set_xlabel("频率 [kHz] @ 48 kHz")
    ax.set_ylabel("相对输入音 [dB]")
    ax.set_title("(e) 复链路 IIR vs FIR：输出线几乎重合（红虚线 = 差频互调位置）", fontsize=9.5)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[2, 1])
    names = ["FIR 复链路", "IIR 复链路"]
    vals = [macs_fir, macs_iir]
    ax.bar(names, vals, color=["0.5", "tab:blue"])
    for i, v in enumerate(vals):
        ax.text(i, v * 1.03, str(v), ha="center", fontsize=10)
    ax.set_ylabel("MAC / 输入采样")
    ax.set_title(f"(f) 代价：IIR 省 {vals[0]/vals[1]:.1f}×（含 I/Q 两路与抽取）", fontsize=9.5)
    ax.grid(alpha=0.3, axis="y")

    fig.suptitle(f"IIR 版（椭圆，阻带 {STOP_DB:.0f} dB）：解析上采样 + 复整形，与 FIR 同架构同信号",
                 fontsize=11)
    out = Path(__file__).parent / "output" / "polyphase_analytic_iir.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=140)
    print(f"\n已写出 {out}")
    print(f"断言失败项：{fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
