"""多相 FIR 超采样：**解析（单边）上采样** + 复整形 vs 传统实数链路。

架构（两条链路共用同一套 remez 规格、同一抽取滤波器）：

    x[n] @ fs_in（实数，两个纯音）
      ├── 链路 R：多相 FIR 插值（实数低通，通带 [0, 0.45·fs_in]、阻带封住镜像）
      │            → 实整形器 P(x) = x + b2·x² + b3·x³
      └── 链路 C：多相 FIR 插值（**解析滤波器**，单边通带 [0, 0.45·fs_in]）
                   → 复整形器 H(z) = z + b2·z² + b3·z³ → 取实部
      → 同一实数 remez 抗混叠低通 → 抽取 ÷L → y[n] @ fs_in

解析上采样滤波器：**一次 remez + 复调制**（不需要复数 remez）：

    h_lp = remez(通带 [0, f_p/2]、阻带 [fs_in − 1.5·f_p, fs_up/2])      （实低通原型）
    h_a  = 2 · h_lp[n] · exp(j·π·f_p·n/fs_up)                        （频移 f_p/2）

于是 D(f) = 2·H_lp(f − f_p/2)：
    f ∈ [0, f_p]            → 2（单边通带；×2 使 Re 还原成带限插值本身，不是它的一半）
    f ∈ [fs_in − f_p, …]    → 0（镜像，阻带电平）
    f ≤ −Δf                 → 0（负频率；Δf = 原型过渡带宽度）
    f ∈ [−Δf, 0]            → 过渡（"守护带"：解析性在 |f| < Δf 内不成立，这是本方法的固有代价）

用法: python qwqdsp/labs/adaa_iir/polyphase_analytic_ovs.py
输出: qwqdsp/labs/adaa_iir/output/polyphase_analytic_ovs.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as ss

# ------------------------------------------------------------
# 参数
# ------------------------------------------------------------
FS_IN = 48000.0
L = 8
FS_UP = FS_IN * L
F1, F2 = 9000.0, 11000.0        # 两个纯音（含高阶产物 27~33 kHz，需要超采样）
A1 = A2 = 0.25
B1, B2, B3 = 1.0, 1.0, 1.0      # b1·z + b2·z² + b3·z³
FP = 0.45 * FS_IN               # 21600：单边通带边（字面 [0, fs_in/2] 要求零宽过渡带，取 0.45）
TAPS = 512                      # = L·64，三条滤波器同一长度
NX = 6000                       # 0.125 s 输入（栅格 8 Hz，9k/11k 都在栅格上）
NDISC, NANA = 2400, 2400        # 丢掉瞬态 2400，分析 2400（栅格 20 Hz）

FUND = {F1, F2}
HARM = {2 * F1, 2 * F2, 3 * F1, 3 * F2}
SUM_IMD = {F1 + F2, 2 * F1 + F2, F1 + 2 * F2}
DIFF_IMD = {abs(F2 - F1), abs(2 * F1 - F2), abs(F1 - 2 * F2)}
ALIAS = {FS_IN - (3 * F1), FS_IN - (2 * F1 + F2), FS_IN - (F1 + 2 * F2), FS_IN - 3 * F2}


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


# ------------------------------------------------------------
# remez 设计
# ------------------------------------------------------------

def design():
    """三个滤波器：解析上采样（复）/ 实数上采样 / 实数抽取。

    补零插值把通带幅度降到 1/L，所以两个上采样滤波器都要在通带乘 L；
    解析滤波器再乘 2（使 Re(z) 还原成带限插值本身，而不是它的一半）。
    """
    # 解析上采样：实低通原型（通带 [0, FP/2]）+ 复频移 FP/2
    h_lp = ss.remez(TAPS, [0, FP / 2, FS_IN - 1.5 * FP, FS_UP / 2], [1, 0], fs=FS_UP)
    n = np.arange(TAPS)
    # 复频移的相位基准取在滤波器中心 n=tau，否则会额外引入一个常数相位
    # （群延迟对、但 Re(z) 会变成旋转过的信号，而不是带限插值本身）。
    tau = (TAPS - 1) / 2.0
    h_a = 2.0 * L * h_lp * np.exp(1j * 2 * np.pi * (FP / 2) * (n - tau) / FS_UP)
    # 实数插值 / 抽取：通带 [0, FP]，阻带自 fs_in−FP 起（抽取时折叠进通带的最近频点是 fs_in−FP）
    h_r = L * ss.remez(TAPS, [0, FP, FS_IN - FP, FS_UP / 2], [1, 0], fs=FS_UP)
    return h_lp, h_a, h_r, ss.remez(TAPS, [0, FP, FS_IN - FP, FS_UP / 2], [1, 0], fs=FS_UP)


# ------------------------------------------------------------
# 多相实现
# ------------------------------------------------------------

def polyphase_upsample(x, h, l):
    """多相插值：y[qL+p] = (h[p::L] * x)[q]。"""
    assert len(h) % l == 0
    outs = [np.convolve(x, h[p::l]) for p in range(l)]
    return np.stack(outs, axis=1).ravel()


def direct_upsample(x, h, l):
    """朴素实现（补零 + 卷积），仅用于验证多相分解。"""
    xu = np.zeros(len(x) * l, dtype=complex if np.iscomplexobj(h) else float)
    xu[::l] = x
    return np.convolve(xu, h)[:len(x) * l]


def aa_decimate(w, h, l):
    """实数抗混叠低通 + 抽取 ÷L。"""
    return np.convolve(w, h)[::l]


# ------------------------------------------------------------
# 测量工具
# ------------------------------------------------------------

def db(x):
    return 20.0 * np.log10(np.maximum(np.abs(x), 1e-300))


def line_levels(seg, fs, freqs):
    """Hann 窗 + 精确 bin 处取峰，归一化到「单位幅度余弦 = 0 dB」。"""
    w = np.hanning(len(seg))
    sp = np.abs(np.fft.rfft(seg * w)) * 2.0 / w.sum()
    out = {}
    for f in freqs:
        i = int(round(f / (fs / len(seg))))
        out[f] = float(db(sp[max(0, i - 2):i + 3].max()))
    return out


def spec(seg, fs):
    w = np.hanning(len(seg))
    sp = np.abs(np.fft.rfft(seg * w)) * 2.0 / w.sum()
    return np.fft.rfftfreq(len(seg), 1.0 / fs), db(sp)


def main() -> int:
    fails = 0
    h_lp, h_a, h_r, h_d = design()

    n = np.arange(NX)
    x = A1 * np.cos(2 * np.pi * F1 * n / FS_IN) + A2 * np.cos(2 * np.pi * F2 * n / FS_IN)

    # ---------- 1. 多相 vs 朴素 ----------
    err_poly = 0.0
    for h, tag in ((h_r, "h_r"), (h_a, "h_a")):
        a, b = polyphase_upsample(x, h, L)[:NX * L], direct_upsample(x, h, L)
        err_poly = max(err_poly, float(np.max(np.abs(a - b))))
    print(f"[1] 多相分解 vs 补零卷积：最大偏差 {err_poly:.2e}（阈值 1e-12）")
    fails += err_poly > 1e-12

    # ---------- 2. 滤波器规格 ----------
    NFFT = 1 << 16
    FA = np.fft.fftshift(np.fft.fftfreq(NFFT, 1.0 / FS_UP))          # 双边频率轴（带符号）
    HA = np.fft.fftshift(np.fft.fft(h_a, NFFT))
    HR = np.fft.fftshift(np.fft.fft(h_r, NFFT))
    HD = np.fft.fftshift(np.fft.fft(h_d, NFFT))
    HL = np.fft.fftshift(np.fft.fft(h_lp, NFFT))
    pb_lp = np.abs(FA) < FP / 2
    pb = np.abs(FA) < FP
    img = (FA > FS_IN - FP) & (FA < FS_IN + FP)
    neg = FA < -2 * FP
    print(f"[2] remez 规格（{TAPS} taps，过渡带 {FS_IN - 2 * FP:.0f} Hz）")
    print(f"    实数插值 h_r（通带增益 L={L}）：通带最大偏差 {db(np.max(np.abs(np.abs(HR[pb]) / L - 1))):.1f} dB，"
          f"阻带（|f| ≥ fs_in−f_p）{db(np.max(np.abs(HR[np.abs(FA) >= FS_IN - FP])) / L):.1f} dB")
    print(f"    实数抽取 h_d（通带增益 1）    ：通带最大偏差 {db(np.max(np.abs(np.abs(HD[pb]) - 1))):.1f} dB")
    print(f"    解析原型 h_lp    ：通带最大偏差 {db(np.max(np.abs(np.abs(HL[pb_lp]) - 1))):.1f} dB，"
          f"阻带 {db(np.max(np.abs(HL[np.abs(FA) >= FS_IN - 1.5 * FP]))):.1f} dB")
    print(f"    解析滤波器 h_a   ：单边通带 [0, f_p] 最大偏差 "
          f"{db(np.max(np.abs(np.abs(HA[(FA >= 0) & (FA < FP)]) / (2 * L) - 1))):.1f} dB，"
          f"镜像带 {db(np.max(np.abs(HA[img])) / (2 * L)):.1f} dB，"
          f"负频率侧 {db(np.max(np.abs(HA[neg])) / (2 * L)):.1f} dB")

    # ---------- 3. 解析上采样正确性 ----------
    z = polyphase_upsample(x, h_a, L)
    t = np.arange(len(z)) / FS_UP
    d = (TAPS - 1) / 2.0                       # 线性相位滤波器的群延迟（采样 @ fs_up）
    ref = A1 * np.exp(2j * np.pi * F1 * (t - d / FS_UP)) + A2 * np.exp(2j * np.pi * F2 * (t - d / FS_UP))
    sl = slice(2000, 20000)                    # 稳态段（瞬态 512 采样已过）
    tol = max(A1, A2)
    err = z[sl] - ref[sl]
    e_rms = db(np.sqrt(np.mean(np.abs(err) ** 2)) / tol)
    e_pk = db(np.max(np.abs(err)) / tol)
    print(f"[3] 解析上采样（群延迟 {d:.1f} 采样 @ fs_up，已补偿）")
    print(f"    z vs 理想解析信号：rms 误差 {e_rms:.1f} dB，峰值 {e_pk:.1f} dB"
          f"（残差全部落在镜像频率 39/57/87… kHz，即滤波器阻带泄漏）")
    # 负频率抑制：取两音公共周期的整数倍数（gcd(9000,11000)=1000 Hz → 384 采样），矩形窗零泄漏
    per = int(round(FS_UP / 1000.0))             # 384
    seg_len = per * 32                           # 12288 采样，已过瞬态
    N0 = 2000
    Zs = np.fft.fftshift(np.fft.fft(z[N0:N0 + seg_len]))
    fzs = np.fft.fftshift(np.fft.fftfreq(seg_len, 1.0 / FS_UP))
    supp = db(np.max(np.abs(Zs[fzs < -2500])) / np.max(np.abs(Zs)))
    print(f"    z 的负频率侧最强分量 {supp:.1f} dB（整周期段，无泄漏；相对通带峰）")
    fails += e_rms > -80.0 or supp > -80.0

    # ---------- 4. 两条链路 ----------
    x_up_r = polyphase_upsample(x, h_r, L)
    w_r = B1 * x_up_r + B2 * x_up_r ** 2 + B3 * x_up_r ** 3
    y_r = aa_decimate(w_r, h_d, L)

    w_c = (B1 * z + B2 * z ** 2 + B3 * z ** 3).real
    y_c = aa_decimate(w_c, h_d, L)

    w_n = B1 * x + B2 * x ** 2 + B3 * x ** 3          # 无超采样（L=1）参照
    y_n = w_n

    seg = lambda y: y[NDISC:NDISC + NANA]
    # 只列 fs_in/2 以内的频率：27~33 kHz 那些在 48 kHz 下根本不存在（被抽取滤波器/折叠处理）
    freqs = sorted(f for f in (FUND | HARM | SUM_IMD | DIFF_IMD | ALIAS) if f < FS_IN / 2)
    lv_r = line_levels(seg(y_r), FS_IN, freqs)
    lv_c = line_levels(seg(y_c), FS_IN, freqs)
    lv_n = line_levels(seg(y_n), FS_IN, freqs)
    ref_db = db(A1)                       # 参考 = 输入音幅度（不用某条链路自己的基波）

    def kind(f):
        return ("基波" if f in FUND else "自身谐波" if f in HARM
                else "求和互调" if f in SUM_IMD else "差频互调" if f in DIFF_IMD else "混叠")

    print(f"\n[4] 输出谱（相对**输入音幅度**，dB；'—' = 未出现（< -140 dB））")
    print(f"    {'频率':>7} {'归属':<9} {'实链路 L=8':>11} {'复链路 L=8':>11} {'实链路 L=1':>11}")
    for f in freqs:
        row = []
        for lv in (lv_r, lv_c, lv_n):
            v = lv[f] - ref_db
            row.append("—" if v < -140 else f"{v:+.1f}")
        print(f"    {f:7.0f} {kind(f):<9} {row[0]:>11} {row[1]:>11} {row[2]:>11}")

    # 判定：L=8 两条链路都不该有混叠；复链路不该有差频互调
    alias_r = max(lv_r[f] - ref_db for f in ALIAS)
    alias_c = max(lv_c[f] - ref_db for f in ALIAS)
    diff_c = max(lv_c[f] - ref_db for f in DIFF_IMD)
    print(f"\n[5] L=8 混叠最强：实链路 {alias_r:+.1f} dB，复链路 {alias_c:+.1f} dB（阈值 -80）")
    print(f"    复链路差频互调最强：{diff_c:+.1f} dB（阈值 -80；实链路同位置 "
          f"{max(lv_r[f] - ref_db for f in DIFF_IMD):+.1f} dB）")
    fails += alias_r > -80.0 or alias_c > -80.0 or diff_c > -80.0

    # ---------- 5. 代价 ----------
    macs_r = TAPS + TAPS              # 插值（实数）+ 抽取
    macs_c = 2 * TAPS + TAPS          # 插值（I/Q 两个实滤波）+ 抽取
    print(f"\n[6] 每输入采样 MAC 数（多相）：实链路 {macs_r}，复链路 {macs_c}（1.5×）"
          f"；每相 {TAPS//L} 抽头")

    # ---------- 图 ----------
    fig = plt.figure(figsize=(13.2, 12.4))
    gs = fig.add_gridspec(3, 2, hspace=0.42, wspace=0.22)

    ax = fig.add_subplot(gs[0, 0])                       # (a) 滤波器响应
    NA = 4096                      # 粗一点，等波纹才看得见（否则被压成色块）
    f2s = np.fft.fftshift(np.fft.fftfreq(NA, 1.0 / FS_UP))
    Hd = np.fft.fftshift(np.fft.fft(np.pad(h_a, (0, NA - len(h_a)))))
    Hr = np.fft.fftshift(np.fft.fft(np.pad(h_r, (0, NA - len(h_r)))))
    Hl = np.fft.fftshift(np.fft.fft(np.pad(h_lp, (0, NA - len(h_lp)))))
    # 先画对称的两条，再把**单边**的解析曲线画到最上层（否则被实数曲线盖住）
    ax.plot(f2s / 1000, db(Hr) - db(L), lw=1.0, alpha=0.45, color="tab:orange",
            label="实数上采样 |H_r|/L（双边）")
    ax.plot(f2s / 1000, db(Hl), lw=0.8, alpha=0.45, ls="--", color="tab:green",
            label="解析原型 |H_lp|（双边）")
    ax.plot(f2s / 1000, db(Hd) - db(2 * L), lw=2.0, color="tab:blue",
            label="解析上采样 |H_a|/(2L)（单边）")
    ax.axvline(-FP / 1000, color="tab:red", ls=":", lw=1.2)
    ax.axvline(FP / 1000, color="tab:red", ls=":", lw=1.2)
    ax.axhline(-100, color="0.7", lw=0.7, ls=":")
    ax.set_xlim(-FS_UP / 2000, FS_UP / 2000)
    ax.set_ylim(-150, 8)
    ax.set_xlabel("频率 [kHz] @ 384 kHz")
    ax.set_ylabel("幅度 [dB]")
    ax.set_title(f"(a) remez 设计：解析滤波器只放正频率（0~{FP/1000:.1f} kHz），负频率压到阻带", fontsize=9.5)
    ax.annotate("", xy=(30, -30), xytext=(-30, -30),
                arrowprops=dict(arrowstyle="-|>", color="tab:blue", lw=1.4))
    ax.text(0, -24, "只此一侧有信号", ha="center", fontsize=8, color="tab:blue")
    ax.legend(fontsize=7.5, loc="lower center", framealpha=0.95)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 1])                       # (b) 时域：z vs 理想解析信号
    t0 = 3000
    idx = np.arange(t0, t0 + 4 * L * 8)                  # 4 个低速率周期
    tt = (idx - d) / FS_UP
    ide = A1 * np.exp(2j * np.pi * F1 * tt) + A2 * np.exp(2j * np.pi * F2 * tt)
    ax.plot(idx - t0, z[idx].real, lw=1.1, color="tab:blue", label="Re(z)")
    ax.plot(idx - t0, ide.real, lw=0.9, ls="--", color="0.35", label="Re(理想解析信号)")
    ax.plot(idx - t0, z[idx].imag, lw=1.1, color="tab:orange", label="Im(z) = Hilbert(s)")
    ax.plot(idx - t0, ide.imag, lw=0.9, ls="--", color="0.6")
    ax.set_xlabel("采样 @ 384 kHz")
    ax.set_ylabel("幅度")
    ax.set_title(f"(b) 解析上采样：Re(z)=带限插值、Im(z)=其 Hilbert（误差 {e_rms:.0f} dB）", fontsize=9.5)
    ax.legend(fontsize=7.5, loc="upper right", ncol=2)
    ax.grid(alpha=0.3)

    def plot_out(ax, y, lv, ref, title):
        f_, s_ = spec(seg(y), FS_IN)
        m = f_ < 25000
        ax.plot(f_[m] / 1000, s_[m] - ref, lw=0.7, color="0.6")
        for f in freqs:
            if f > 24000 or lv[f] - ref < -140:
                continue
            ax.vlines(f / 1000, -140, lv[f] - ref, color=color_of(f), lw=1.5)
        ax.axvspan(FS_IN / 2000, 25000 / 1000, color="tab:red", alpha=0.06)
        ax.set_xlim(0, 25000 / 1000)
        ax.set_ylim(-140, 6)
        ax.set_xlabel("频率 [kHz] @ 48 kHz")
        ax.set_ylabel("相对输入音 [dB]")
        ax.set_title(title, fontsize=9.5)
        ax.grid(alpha=0.3)

    plot_out(fig.add_subplot(gs[1, 0]), y_r, lv_r, ref_db,
             f"(c) 实链路 @L={L}：无混叠，但有差频互调 2k/7k/13k")
    plot_out(fig.add_subplot(gs[1, 1]), y_c, lv_c, ref_db,
             f"(d) 复链路 @L={L}：无混叠、**无差频互调**")
    plot_out(fig.add_subplot(gs[2, 0]), y_n, lv_n, ref_db,
             "(e) 实链路 @L=1（无超采样）：27~33 kHz 折回 15~21 kHz（紫）")

    ax = fig.add_subplot(gs[2, 1])                       # (f) 叠加对照
    for y, ref, lab, c in ((y_n, ref_db, "实链路 L=1", "tab:purple"),
                           (y_r, ref_db, f"实链路 L={L}", "0.45"),
                           (y_c, ref_db, f"复链路 L={L}", "tab:green")):
        f_, s_ = spec(seg(y), FS_IN)
        m = f_ < 25000
        ax.plot(f_[m] / 1000, s_[m] - ref, lw=0.7, label=lab, color=c, alpha=0.85)
    for f in sorted(DIFF_IMD):
        ax.axvline(f / 1000, color="tab:red", ls=":", lw=1.0)
    ax.set_xlim(0, 25)
    ax.set_ylim(-120, 6)
    ax.set_xlabel("频率 [kHz] @ 48 kHz")
    ax.set_ylabel("相对输入音 [dB]")
    ax.set_title("(f) 三条输出谱叠加（红虚线 = 差频互调位置）", fontsize=9.5)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)

    fig.suptitle(f"多相 FIR 超采样（L={L}，8/16 kHz 双音，remez 设计）：解析上采样 + 复整形 vs 实数链路", fontsize=11)
    out = Path(__file__).parent / "output" / "polyphase_analytic_ovs.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=140)
    print(f"\n已写出 {out}")
    print(f"断言失败项：{fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
