"""直接项研究 —— 偶数阶椭圆低通的分解 / 类连续频谱。

第 1 步：偶数阶椭圆（阻带 -80 dB）做部分分式分解
    H(s) = Σ_m r_m/(s - p_m) + d
比较三种响应：
    H_orig   : 原始 b/a 形式（双真，|H(j∞)| = d ≠ 0）
    H_rebuild: 由部分分式重建 Σ r_m/(s-p_m) + d  —— 应与 H_orig 逐点相同
    H_proper : 丢掉直接项 Σ r_m/(s-p_m)           —— 严格真，阻带谷底被抬到 d

第 2 步：极低采样率纯音 → 高倍上采样，看"类连续"频谱
    (1) 零插值（= 冲激串）→ 频谱是采样谱以 f_in 为周期重复（镜像）
    (2) 频域砖墙      → 理想带限重建，只剩基音
    (3) 过 H / 过 H_p → 镜像被压到 d 量级

用法: python qwqdsp/labs/holters_parker/direct_term_study.py [1|2|all]
输出: qwqdsp/labs/holters_parker/output/even_elliptic_decompose.png
      qwqdsp/labs/holters_parker/output/step2_quasi_continuous.png
"""
import sys
from fractions import Fraction
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as ss

OUT_DIR = Path(__file__).parent / "output"

ORDER = 6          # 偶数阶（3 对共轭极点）
RIPPLE_DB = 0.1    # 通带纹波
STOP_DB = 80.0     # 阻带涟漪（-80 dB）


def partial_fraction(b, a):
    """部分分式展开: H(s) = Σ r_m/(s - p_m) + d"""
    r, p, d = ss.residue(b, a)
    return r, p, (d[0] if len(d) else 0.0)


def eval_partial(w, r, p, d):
    """复数响应 Σ r_m/(jw - p_m) + d"""
    s = 1j * np.asarray(w, dtype=float)
    return np.sum(r[:, None] / (s[None, :] - p[:, None]), axis=0) + d


def db(x, floor=1e-20):
    return 20 * np.log10(np.maximum(np.abs(x), floor))


# ------------------------------------------------------------
# 第 1 步
# ------------------------------------------------------------
def step1():
    b, a = ss.ellip(ORDER, RIPPLE_DB, STOP_DB, 1.0, analog=True)
    r, p, d = partial_fraction(b, a)
    print(f"偶数阶椭圆: N={ORDER}, 通带 {RIPPLE_DB} dB, 阻带 {STOP_DB} dB")
    print(f"  阶数(分子/分母) = {len(b)-1}/{len(a)-1}   直接项 d = b0/a0 = {b[0]/a[0]:.12g}")
    print(f"  residue() 给出的直接项 = {d:.12g}")
    print("  极点/留数（上半平面）:")
    for pole, res in zip(p, r):
        if pole.imag > 0:
            print(f"    p = {pole:+.6f}   r = {res:+.6f}")

    w = np.logspace(-3, 3, 4001)
    h_orig = ss.freqs(b, a, w)[1]
    h_rebuild = eval_partial(w, r, p, d)
    h_proper = eval_partial(w, r, p, 0.0)

    err = np.abs(h_rebuild - h_orig)
    print(f"\n  重建与原始的最大差: max|H_rebuild - H_orig| = {err.max():.3e}")
    print(f"  在 w->∞: |H_orig| = {abs(ss.freqs(b, a, np.array([1e9]))[1][0]):.6e},"
          f"  |H_proper| = {abs(eval_partial(np.array([1e9]), r, p, 0.0)[0]):.3e}")

    fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.4))

    # (a) 三条响应
    ax[0].plot(w, db(h_orig), lw=2.4, color="0.55", label="H_orig  (b/a, 双真)")
    ax[0].plot(w, db(h_rebuild), lw=1.2, ls="--", color="C0", label="H_rebuild = Σr/(s−p) + d")
    ax[0].plot(w, db(h_proper), lw=1.2, ls=":", color="C3", label="H_proper = Σr/(s−p)  (丢掉 d)")
    ax[0].axhline(db(np.array([d]))[0], color="C3", lw=0.8, alpha=0.5)
    ax[0].text(1.4e-3, db(np.array([d]))[0] + 3, f"d = {d:.1e}  ({20*np.log10(d):.1f} dB)",
               color="C3", fontsize=9)
    ax[0].set_xscale("log")
    ax[0].set_ylim(-110, 5)
    ax[0].set_xlabel("ω (rad/s)")
    ax[0].set_ylabel("|H| (dB)")
    ax[0].set_title("(a) 分解前 vs 重建后")
    ax[0].legend(fontsize=8, loc="lower left")
    ax[0].grid(alpha=0.3, which="both")

    # (b) 重建误差（应为机器精度）
    ax[1].semilogx(w, np.maximum(err / np.maximum(np.abs(h_orig), 1e-20), 1e-18), lw=1.0, color="C2")
    ax[1].set_ylim(1e-18, 1e-12)
    ax[1].set_xlabel("ω (rad/s)")
    ax[1].set_ylabel("|H_rebuild − H_orig| / |H_orig|")
    ax[1].set_title(f"(b) 重建误差 = 机器精度 (max {err.max():.1e})")
    ax[1].grid(alpha=0.3, which="both")

    # (c) 线性放大：传输零点（谷=0）与 H_proper 在那里抬到 d
    wz = np.logspace(0.25, 1.6, 4001)          # ≈1.8 … 40 rad/s，覆盖 2.19/2.92/7.76 三个零点
    hz_orig = np.abs(ss.freqs(b, a, wz)[1])
    hz_prop = np.abs(eval_partial(wz, r, p, 0.0))
    ax[2].plot(wz, hz_orig * 1e4, lw=1.2, color="0.55", label="H_orig  (×1e4)")
    ax[2].plot(wz, hz_prop * 1e4, lw=1.2, ls=":", color="C3", label="H_proper = H_orig − d (×1e4)")
    ax[2].axhline(d * 1e4, color="C0", lw=0.8, ls="--")
    ax[2].text(wz[1], d * 1e4 * 1.15, "d (×1e4)", color="C0", fontsize=9)
    ax[2].set_xscale("log")
    ax[2].set_xlabel("ω (rad/s)")
    ax[2].set_ylabel("|H| ×1e4  (0 = −∞ dB)")
    ax[2].set_title("(c) 阻带放大: 涟漪峰 ≈ d，谷 = 0；H_proper 把谷抬到 d")
    ax[2].legend(fontsize=8, loc="upper right")
    ax[2].grid(alpha=0.3, which="both")

    fig.suptitle(f"偶数阶椭圆 N={ORDER}, rp={RIPPLE_DB} dB, rs={STOP_DB} dB — "
                 f"d = |H(j∞)| = {d:.3g}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    OUT_DIR.mkdir(exist_ok=True)
    out = OUT_DIR / "even_elliptic_decompose.png"
    fig.savefig(out, dpi=130)
    print(f"\n图已保存: {out}")


# ------------------------------------------------------------
# 第 2 步
# ------------------------------------------------------------
def build(fs_in=20.0, tone=2.5, n=512, upsample=8192):
    """构造类连续域里的各条信号（step2 / step3 共用）。"""
    k = np.arange(n)
    x = np.sin(2.0 * np.pi * tone * k / fs_in)
    fs_fine = fs_in * upsample
    m = n * upsample
    f = np.fft.rfftfreq(m, 1.0 / fs_fine)

    # 零插值。两个用途要用**不同**的幅度约定：
    #   xz_filter : 喂给滤波器的冲激串表示 —— 连续模型里狄拉克权重是 T·x[k]，
    #               放到间距 h=T/L 的细网格上，样本值就是 (T/h)·x[k] = L·x[k]
    #   xz_plain  : 「插零后的高采样信号」本身 —— 值是 x[k]、别处为 0（幅度就是 1）
    #               直接项 d 乘的是它，所以栅格上的尖峰是 d·x[k]（= 1·d）
    xz_filter = np.zeros(m)
    xz_filter[::upsample] = x * upsample
    xz_plain = np.zeros(m)
    xz_plain[::upsample] = x
    Xz = np.fft.rfft(xz_filter)

    # 偶数阶椭圆，截止摆到 fs_in/2
    b, a = ss.ellip(ORDER, RIPPLE_DB, STOP_DB, 1.0, analog=True)
    r, p, d = partial_fraction(b, a)
    wc = 2.0 * np.pi * (fs_in / 2)
    H = ss.freqs(b, a, 2.0 * np.pi * f / wc)[1]

    t = np.arange(m) / fs_fine
    x_buf = np.sin(2.0 * np.pi * tone * t)          # 纯音的带限重构（解析已知）
    y_nod = np.fft.irfft(Xz * (H - d), n=m)         # 无 d：H_p 通路
    y_phys = y_nod + d * x_buf                      # 参照：d 乘在插值后的信号上
    y_spike = y_nod + d * xz_plain                  # ★ 直接项 = d × 插零信号（值 x[k]）
    return {"fs_in": fs_in, "tone": tone, "n": n, "upsample": upsample, "fs_fine": fs_fine, "m": m,
            "f": f, "x": x, "xz": xz_plain, "xz_filter": xz_filter, "Xz": Xz, "H": H, "d": d, "t": t,
            "x_buf": x_buf, "y_nod": y_nod, "y_phys": y_phys, "y_spike": y_spike,
            "X_buf": np.fft.rfft(d * x_buf)}


def step2(fs_in=20.0, tone=2.5, n=512, upsample=8192):
    """极低采样率纯音 → 高倍上采样，看"类连续"频谱。

    tone 取 fs_in/8：n 点里正好 n/8 个整周期，且落在细网格的 bin 上（无泄漏）。
    零插值乘 upsample：这是 L 倍插值的标准约定（否则通带幅度只有 1/L）。
    upsample 越大，细率 fs_fine 越高、full band 越宽，而镜像间距恒为 fs_in
    → 同一个横轴上的竖线数 ∝ upsample，最终糊成一片面积。
    """
    g = build(fs_in, tone, n, upsample)
    f, m, fs_fine, d = g["f"], g["m"], g["fs_fine"], g["d"]
    Xz, H, xz, x_buf, X_buf = g["Xz"], g["H"], g["xz"], g["x_buf"], g["X_buf"]
    t, y_nod, y_phys, y_spike = g["t"], g["y_nod"], g["y_phys"], g["y_spike"]
    y_ref = np.fft.irfft(Xz * H, n=m)
    print(f"  修复前后的差: max|y_spike(旧, 用 L·x[k]) − y_spike(新, 用 x[k])| = "
          f"{np.abs(y_spike - y_ref).max():.3e}  = d·(L−1) ≈ {d*(upsample-1):.3e}")
    print(f"  直接项量级: d·x_buf 振幅 = {d:.3e}（= 1·d）；"
          f"d×插零信号在栅格上的尖峰 = d·x[k] = {d:.3e}（同样是 1·d）✓")
    print(f"  频谱上的差别: d·冲激串 在镜像处也是 {20*np.log10(d):.1f} dB；"
          f"d·x_buf 在镜像处 = {20*np.log10(max(np.abs(X_buf[f>5*fs_in]).max(),1e-30)/np.max(np.abs(X_buf))):.0f} dB（无镜像）")

    span = 10 * fs_in                      # 图 A 用原来的范围（10 个采样率周期）
    idx = f <= span
    ref = np.max(np.abs(Xz))
    Xd = Xz * H
    Xp = Xz * (H - d)

    fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.4))

    # (a) 类连续频谱（无砖墙，叠加施加的滤波器幅度响应）
    ax[0].plot(f[idx], db(Xz[idx] / ref), lw=0.9, color="0.6", label=f"冲激串（零插值 ×{upsample}）")
    ax[0].plot(f[idx], db(H[idx]), lw=1.4, color="k", label="施加的滤波器 |H(jω)|（已缩放到截止 = fs_in/2）")
    ax[0].plot(f[idx], db(Xp[idx] / ref), lw=1.0, ls="--", color="C0", label="过 H_p（无 d）")
    ax[0].plot(f[idx], db(Xd[idx] / ref), lw=1.0, color="C3", label=f"= H_p + d·冲激串（含 d={d:.0e}）")
    ax[0].plot(f[idx], db(d * Xz[idx] / ref), lw=1.0, ls=":", color="C2", label="d·冲激串（差项）")
    ax[0].axhline(20 * np.log10(d), color="C3", lw=0.8, alpha=0.5)
    ax[0].axvline(fs_in / 2, color="0.3", ls=":", lw=1.0)
    ax[0].text(fs_in / 2 * 1.05, -8, f"fs_in/2 = {fs_in/2:g} Hz", fontsize=8, color="0.3")
    for j in (1, 2, 3):
        ax[0].axvline(j * fs_in - tone, color="0.88", lw=0.6)
        ax[0].axvline(j * fs_in + tone, color="0.88", lw=0.6)
    ax[0].set_xlim(0, span)
    ax[0].set_ylim(-110, 6)
    ax[0].set_xlabel("频率 (Hz)")
    ax[0].set_ylabel("dB（相对基音）")
    ax[0].set_title(f"(a) “类连续”频谱: fs_in={fs_in:g} Hz, 纯音 {tone:g} Hz")
    ax[0].legend(fontsize=8, loc="lower right")
    ax[0].grid(alpha=0.3)

    # (b) 时域: 无 d 版本 vs HP 模型的直接项（两者在图上重合）
    t = np.arange(m) / fs_fine
    sel = t <= 1.5
    ax[1].plot(t[sel], y_nod[sel], lw=1.6, color="C0", label="无 d 版本（H_p 通路）")
    ax[1].plot(t[sel], y_spike[sel], lw=0.8, ls="--", color="C3",
               label=f"+ d × 插零信号（振幅 {d:.0e}，图上重合）")
    ax[1].set_xlim(0, 1.5)
    ax[1].set_xlabel("时间 (s)")
    ax[1].set_ylabel("幅度")
    ax[1].set_title(f"(b) 时域: 直接项振幅 1·d，且只占 {100/upsample:.2f}% 的样本 → 看不出")
    ax[1].legend(fontsize=8, loc="upper right")
    ax[1].grid(alpha=0.3)

    # (c) 直接项本身放大 1e4 倍: HP 模型（只在栅格点）vs 带限模型（光滑）
    t0 = 0.5
    win = 0.025
    zs = (t > t0 - win) & (t < t0 + win)
    ax[2].plot(t[zs], d * x_buf[zs] * 1e4, lw=1.4, ls="--", color="0.55",
               label="参照：d × 带限插值（另一种模型，处处光滑）")
    nz = zs & (xz != 0)
    ax[2].plot(t[nz], d * xz[nz] * 1e4, "o", ms=5, color="C3",
               label="d × 插零信号（HP 模型，只在栅格点，振幅 = 1·d）")
    ax[2].axhline(0, color="0.7", lw=0.6)
    ax[2].axvline(t0, color="0.6", ls="--", lw=0.8)
    ax[2].set_xlim(t0 - win, t0 + win)
    ax[2].set_ylim(-1.25, 1.25)
    ax[2].set_xlabel("时间 (s)")
    ax[2].set_ylabel("幅度 ×1e4")
    ax[2].set_title("(c) 直接项（放大 1e4）: HP 模型只在栅格点有值")
    ax[2].legend(fontsize=8, loc="upper right")
    ax[2].grid(alpha=0.3)

    fig.suptitle(f"图 A：类连续频谱（10 个采样率周期）与直接项 "
                 f"(fs_in={fs_in:g} Hz, tone={tone:g} Hz, 上采样 ×{upsample})", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    OUT_DIR.mkdir(exist_ok=True)
    out = OUT_DIR / "step2_quasi_continuous.png"
    fig.savefig(out, dpi=130)

    # ---- 图 B：整个上采样带宽（0 … fs_fine/2），横轴对数 ----
    fig2, axb = plt.subplots(1, 2, figsize=(15.0, 4.4))
    f0 = f[1]                                   # 最低非零 bin（对数轴不能含 0）
    axb[0].plot(f, db(Xz / ref), lw=0.35, color="0.6", alpha=0.55,
                label=f"冲激串（×{upsample}）")
    axb[0].plot(f, db(H), lw=1.6, color="k", label="施加的滤波器 |H(jω)|")
    axb[0].plot(f, db(Xp / ref), lw=0.9, ls="--", color="C0", label="过 H_p（无 d）")
    axb[0].plot(f, db(Xd / ref), lw=0.9, color="C3", label=f"= H_p + d·冲激串 (d={d:.0e})")
    axb[0].plot(f, db(X_buf / ref), lw=1.2, ls="-.", color="C2",
                label="d·x_buf（正确的直接项：只剩基带，无镜像）")
    axb[0].axhline(20 * np.log10(d), color="C3", lw=0.8, alpha=0.5)
    axb[0].axvline(fs_in / 2, color="0.3", ls=":", lw=1.0)
    axb[0].text(fs_in / 2 * 1.15, -100, f"fs_in/2={fs_in/2:g}", fontsize=8, color="0.3")
    axb[0].set_xscale("log")
    axb[0].set_xlim(f0, fs_fine / 2)
    axb[0].set_ylim(-120, 6)
    axb[0].set_xlabel("频率 (Hz, 对数)")
    axb[0].set_ylabel("dB（相对基音）")
    axb[0].set_title(f"(a) 整个上采样带宽 对数轴: {f0:.3g} Hz … {fs_fine/2/1e3:.2f} kHz")
    axb[0].legend(fontsize=8, loc="lower left")
    axb[0].grid(alpha=0.3, which="both")

    axb[1].plot(f, np.abs(Xz) / ref, lw=0.35, color="0.6", alpha=0.55, label="冲激串")
    axb[1].plot(f, np.abs(H), lw=1.6, color="k", label="|H(jω)|")
    axb[1].plot(f, np.abs(Xp) / ref, lw=0.9, ls="--", color="C0", label="过 H_p（无 d）")
    axb[1].plot(f, np.abs(Xd) / ref, lw=0.9, color="C3", label="= H_p + d·冲激串")
    axb[1].axhline(d, color="C3", lw=0.8, alpha=0.5)
    axb[1].set_xscale("log")
    axb[1].set_xlim(f0, fs_fine / 2)
    axb[1].set_ylim(0, 1.05)
    axb[1].set_xlabel("频率 (Hz, 对数)")
    axb[1].set_ylabel("线性幅度（相对基音）")
    axb[1].set_title(f"(b) 同一频段线性幅度: 基音 1.0，其余 ≤ {d:.0e}")
    axb[1].legend(fontsize=8, loc="upper right")
    axb[1].grid(alpha=0.3, which="both")

    fig2.suptitle(f"图 B：整个上采样带宽（对数轴 {f0:.3g} Hz … {fs_fine/2/1e3:.2f} kHz）", fontsize=12)
    fig2.tight_layout(rect=(0, 0, 1, 0.93))
    out2 = OUT_DIR / "step2_full_band.png"
    fig2.savefig(out2, dpi=130)
    print(f"  图 A（10 个周期）: {out}")
    print(f"  图 B（整个带宽 0…{fs_fine/2:g} Hz）: {out2}")
    print(f"  镜像位置 f = k·{fs_in:g} ± {tone:g} Hz: "
          + ", ".join(f"{j*fs_in-tone:g}/{j*fs_in+tone:g}" for j in range(1, 4)) + " …")
    print(f"  冲激串的镜像与基音**等幅**（所以滤波器必须把它们全部压掉）")
    print(f"  过 H 后: 远镜像 ≈ d = {20*np.log10(d):.1f} dB，但第一个镜像 {fs_in-tone:g} Hz "
          f"落在过渡带(阻带边沿 {fs_in/2*2.126:.1f} Hz) → 只有约 -46 dB")
    print(f"  过 H_p（丢 d）与过 H 在此处差别很小 —— 镜像抗性主要由 rs/过渡带决定，不是由 d 决定")


def step3(fs_in=20.0, tone=2.5, n=512, upsample=8192, fs_out=8.0):
    """把第 2 步构造的类连续信号抽取到低速率，再看输出频谱。

    第 2 步那个类连续信号就是

        y_spike = 「无 d 滤波器的输出」 + 「d × 插零后的高采样信号」

    这里直接在细网格上抽取（每 D 个取一个）= 在低速率的时刻采样这个连续信号，
    再与「无 d 版本 y_nod」对比：**差值就是「d × 插零信号」这一项在输出率上的贡献**。
    """
    g = build(fs_in, tone, n, upsample)
    fs_fine = g["fs_fine"]
    D = int(round(fs_fine / fs_out))
    if abs(D * fs_out - fs_fine) > 1e-6:
        raise ValueError(f"fs_out={fs_out} 必须整除 fs_fine={fs_fine}")

    y_nod_out = g["y_nod"][::D]
    y_spike_out = g["y_spike"][::D]

    # 记录长度取整数个基音周期，避免 FFT 泄漏
    period = next(kk for kk in range(1, 100000)
                  if abs(kk * tone / fs_out - round(kk * tone / fs_out)) < 1e-9)
    tot = (len(y_nod_out) // period) * period
    y_nod_out, y_spike_out = y_nod_out[:tot], y_spike_out[:tot]

    fo = np.fft.rfftfreq(tot, 1.0 / fs_out)
    t_out = np.arange(tot) / fs_out
    Yn, Ys = np.fft.rfft(y_nod_out), np.fft.rfft(y_spike_out)
    Ydiff = Ys - Yn
    ref = np.max(np.abs(Yn))
    di = int(np.argmin(np.abs(fo - tone)))
    hits = (np.arange(tot) * D) % upsample == 0

    print(f"  抽取 {fs_in:g} → {fs_out:g} Hz（细率 {fs_fine:g} Hz，D = {D}）")
    print(f"  输出点数 {tot}（{tot/fs_out:.1f} s = {period} 的整数倍）；"
          f"其中 {100*hits.mean():.1f}% 落在输入栅格上")
    print(f"  基音 {tone:g} Hz: 无 d {20*np.log10(abs(Yn[di])/ref):+.2f} dB, "
          f"有 d·插零 {20*np.log10(abs(Ys[di])/ref):+.2f} dB, "
          f"差值 {20*np.log10(abs(Ydiff[di])/ref):+.2f} dB   （d 本身 = {20*np.log10(g['d']):.1f} dB）")
    print("  差值谱（d×插零 这一项在输出率上的全部内容）:")
    for i in np.argsort(np.abs(Ydiff))[-6:][::-1]:
        print(f"    {fo[i]:6.2f} Hz: {20*np.log10(max(abs(Ydiff[i]),1e-30)/ref):7.1f} dB")

    fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.4))

    # (a) 输出频谱
    ax[0].plot(fo, db(Yn / ref), lw=1.0, color="C0", label="无 d（H_p 输出）")
    ax[0].plot(fo, db(Ys / ref), lw=1.4, color="C3", label="无 d + d×插零信号")
    ax[0].axvline(tone, color="0.6", ls=":", lw=0.8)
    ax[0].set_xlim(0, fs_out / 2)
    ax[0].set_ylim(-120, 12)
    ax[0].set_xlabel(f"频率 (Hz) @ 输出率 {fs_out:g} Hz")
    ax[0].set_ylabel("dB（相对基音）")
    ax[0].set_title(f"(a) 抽取到 {fs_out:g} Hz 后的频谱")
    ax[0].legend(fontsize=8, loc="lower right")
    ax[0].grid(alpha=0.3)

    # (b) 差值谱 = d×插零这一项在输出率上的内容
    ax[1].plot(fo, db(Ydiff / ref), lw=1.2, color="C2")
    ax[1].axvline(tone, color="0.6", ls=":", lw=0.8)
    ax[1].text(tone, -114, f" {tone:g} Hz（基音）", fontsize=8, color="0.4", rotation=90)
    ax[1].set_xlim(0, fs_out / 2)
    ax[1].set_ylim(-120, 0)
    ax[1].set_xlabel(f"频率 (Hz) @ 输出率 {fs_out:g} Hz")
    ax[1].set_ylabel("dB（相对基音）")
    ax[1].set_title("(b) 差值谱: d×插零 在输出率上贡献了什么")
    ax[1].grid(alpha=0.3)

    # (c) 输出时域
    nt = min(24, tot)
    ax[2].stem(t_out[:nt], y_spike_out[:nt], linefmt="C3-", markerfmt="C3o", basefmt=" ",
               label="无 d + d×插零")
    ax[2].plot(t_out[:nt], y_nod_out[:nt], ".-", ms=4, lw=0.8, color="C0", label="无 d")
    ax[2].set_xlim(-0.05, nt / fs_out)
    ax[2].set_xlabel("时间 (s)")
    ax[2].set_ylabel("幅度")
    ax[2].set_title(f"(c) 输出时域: {100*hits.mean():.0f}% 的样本落在输入栅格上")
    ax[2].legend(fontsize=8, loc="upper right")
    ax[2].grid(alpha=0.3)

    fig.suptitle(f"第 3 步：把第 2 步的类连续信号抽取到 {fs_out:g} Hz "
                 f"（{fs_in:g}→{fs_out:g}，D={D}）", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    OUT_DIR.mkdir(exist_ok=True)
    out = OUT_DIR / "step3_resampled.png"
    fig.savefig(out, dpi=130)
    print(f"图已保存: {out}")


def step4(fs_in=20.0, fs_out=8.0, tone=2.5, n_in=512):
    """直接构建 Holters–Parker 重采样器（按用户指定的三步）。

    1) 用同一个偶数阶椭圆原型的**严格真部分** H_p = Σ r_m/(s-p_m)，用「一阶复极点滤波器
       状态迭代」求无 d 输出并在新采样率的时刻取值（即 algorithm.md §4 的循环）；
    2) fs_in : fs_out = p : q（既约）→ 输出第 q·n 个点恰好落在输入第 p·n 个点上
       （因为 t = q·n·(p/q) = p·n），于是直接叠加 d·x[p·n]；
    3) 先算无 d 的 y，再叠加直接项。
    """
    b, a = ss.ellip(ORDER, RIPPLE_DB, STOP_DB, 1.0, analog=True)
    r, p, d = partial_fraction(b, a)
    cutoff = min(fs_in, fs_out) / 2
    w = 2.0 * np.pi * cutoff / fs_in          # 通带边沿 rad/sample（原型 fpass 归一化在 1）
    lam = p * w                               # 单极点（rad/sample）
    c = r * w                                 # 对应留数
    pbar = np.exp(lam)

    x = np.sin(2.0 * np.pi * tone * np.arange(n_in) / fs_in)
    inc = fs_in / fs_out
    fr = Fraction(int(round(fs_in)), int(round(fs_out)))
    pn, qn = fr.numerator, fr.denominator

    # --- 1) 无 d：状态迭代 + 任意时刻取值 ---
    state = c * x[0]
    rpos, phase = 0, 0.0
    y = []
    while rpos < len(x) - 1:
        y.append(float(np.sum((state * np.exp(phase * lam)).real)))
        phase += inc
        new_rpos = min(rpos + int(np.floor(phase)), len(x) - 1)
        phase -= np.floor(phase)
        for kk in range(rpos + 1, new_rpos + 1):
            state = state * pbar + c * x[kk]
        rpos = new_rpos
    y_nod = np.array(y)

    # --- 2)+3) 叠加直接项：输出第 q·n 点 ← 输入第 p·n 点 ---
    y_fix = y_nod.copy()
    hit_out = [qn * nn for nn in range((len(y_nod) - 1) // qn + 1) if pn * nn < n_in]
    for i_out in hit_out:
        y_fix[i_out] += d * x[pn * (i_out // qn)]

    pn_, qn_ = pn, qn
    print(f"  HP 重采样 {fs_in:g} → {fs_out:g} Hz（fs_in/fs_out = {pn_}/{qn_}，截止 {cutoff:g} Hz）")
    print(f"  输出点数 {len(y_nod)}；命中点 {len(hit_out)} 个（每 {qn_} 个输出一个，命中率 "
          f"{100*len(hit_out)/len(y_nod):.1f}%）")

    # --- 评估：输出频谱 vs 解析理想（跳过起始暂态，避免泄漏假线）---
    period = next(kk for kk in range(1, 100000)
                  if abs(kk * tone / fs_out - round(kk * tone / fs_out)) < 1e-9)
    skip = period * 4                      # 跳过 4 个周期长度的暂态
    tot = ((len(y_nod) - skip) // period) * period
    yo_nod, yo_fix = y_nod[skip:skip + tot], y_fix[skip:skip + tot]
    fo = np.fft.rfftfreq(tot, 1.0 / fs_out)
    Yn, Yf = np.fft.rfft(yo_nod), np.fft.rfft(yo_fix)
    # 原型通带边沿归一化在 1 rad/s，而目标截止对应 ω=1 → 基音处 ω = tone/cutoff
    Ht = complex(ss.freqs(b, a, np.array([tone / cutoff]))[1][0])
    di = int(np.argmin(np.abs(fo - tone)))
    ref = abs(Yn[di])
    print(f"  分析段 {tot} 点（跳过前 {skip} 点暂态） 基音 {tone:g} Hz: "
          f"无 d {20*np.log10(abs(Yn[di])/ref):+.2f} dB, "
          f"+命中点补 d {20*np.log10(abs(Yf[di])/ref):+.2f} dB, "
          f"解析 |H| = {20*np.log10(abs(Ht)):+.2f} dB")
    print("  输出谱里最显著的几根线（相对基音）:")
    for i in np.argsort(np.abs(Yf))[-6:][::-1]:
        print(f"    {fo[i]:6.2f} Hz: {20*np.log10(max(abs(Yf[i]),1e-30)/ref):7.1f} dB")

    fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.4))
    t_out = np.arange(tot) / fs_out
    ax[0].plot(fo, db(Yn / ref), lw=1.0, color="C0", label="无 d（状态迭代输出）")
    ax[0].plot(fo, db(Yf / ref), lw=1.4, color="C3", label="+ d·x[p·n] 于第 q·n 点")
    ax[0].set_xlim(0, fs_out / 2)
    ax[0].set_ylim(-120, 12)
    ax[0].set_xlabel(f"频率 (Hz) @ {fs_out:g} Hz")
    ax[0].set_ylabel("dB（相对基音）")
    ax[0].set_title(f"(a) HP 重采样器输出谱（{fs_in:g}→{fs_out:g}）")
    ax[0].legend(fontsize=8, loc="lower right")
    ax[0].grid(alpha=0.3)

    ax[1].plot(fo, db((Yf - Yn) / ref), lw=1.2, color="C2")
    ax[1].set_xlim(0, fs_out / 2)
    ax[1].set_ylim(-120, 0)
    ax[1].set_xlabel(f"频率 (Hz) @ {fs_out:g} Hz")
    ax[1].set_ylabel("dB（相对基音）")
    ax[1].set_title("(b) 直接项（命中点补法）在输出率上的内容")
    ax[1].grid(alpha=0.3)

    nt = min(32, tot)
    ax[2].plot(t_out[:nt], yo_nod[:nt], ".-", ms=4, lw=0.8, color="C0", label="无 d")
    ax[2].plot(t_out[:nt], yo_fix[:nt], ".", ms=7, color="C3", label="+ d·x[p·n]")
    ax[2].set_xlim(-0.05, nt / fs_out)
    ax[2].set_xlabel("时间 (s)")
    ax[2].set_ylabel("幅度")
    ax[2].set_title(f"(c) 输出时域: 命中点每 {qn_} 个输出一个")
    ax[2].legend(fontsize=8, loc="upper right")
    ax[2].grid(alpha=0.3)

    fig.suptitle(f"第 4 步：HP 重采样器 + 精确命中点直接项（{fs_in:g}→{fs_out:g}，{pn_}/{qn_}）",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    OUT_DIR.mkdir(exist_ok=True)
    out = OUT_DIR / "step4_hp_resampler.png"
    fig.savefig(out, dpi=130)
    print(f"图已保存: {out}")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("1", "all"):
        print("=" * 70)
        print("第 1 步：偶数阶椭圆的部分分式分解与重建")
        print("=" * 70)
        step1()
    if which in ("2", "all"):
        print("\n" + "=" * 70)
        print("第 2 步：极低采样率纯音的“类连续”频谱")
        print("=" * 70)
        step2()
    if which in ("3", "all"):
        print("\n" + "=" * 70)
        print("第 3 步：重采样到其他采样率后的输出频谱")
        print("=" * 70)
        step3()
    if which in ("4", "all"):
        print("\n" + "=" * 70)
        print("第 4 步：直接构建 HP 重采样器（状态迭代 + 精确命中点补 d）")
        print("=" * 70)
        step4()


if __name__ == "__main__":
    main()
