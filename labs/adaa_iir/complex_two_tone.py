"""复数整形器的**双音**实验：直接按上采样速率合成解析信号（两个纯音），
配一个「只有 2 阶和 3 阶」的解析整形器 H(z) = b1·z + b2·z² + b3·z³，看输出长什么样。

背景（Vicanek, *Complex Waveshapers*, 2025）：H 的 Taylor 系数就是谐波幅度，
`F = Re H(z)/r`、`G = Im H(z)/r` 是 Hilbert 对。关键性质：**解析信号只有正频率，
解析函数只产生频率求和**——所以差频（经典 IMD：f2−f1、2f1−f2、f1−2f2、DC）**不生成**。

双音下 `z = r1·e^{iω1t} + r2·e^{iω2t}`，`z^n` 的多项式展开给出**多重和**
（k1·ω1 + k2·ω2, k1+k2 = n），于是：

  1 阶：f1、f2
  2 阶：2f1、2f2、**f1+f2**（交叉项 2·b2·r1r2）
  3 阶：3f1、3f2、**2f1+f2、f1+2f2**（交叉项 3·b3·r1²r2、3·b3·r1r2²）

三条数值结论：
  * 复链路**没有任何差频/直流**（红标一次都不出现），实链路同多项式则一大堆；
  * 复链路里**每个音自己的谐波序列不受另一个音影响**（2f1 幅度恒为 b2·r1²，与 r2 无关）；
  * 复链路的新增分量全部**高于两个载频**，实链路的差频**低于载频**（因此过任何输出滤波器都活着）。

用法: python qwqdsp/labs/adaa_iir/complex_two_tone.py
输出: qwqdsp/labs/adaa_iir/output/complex_two_tone.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np

FS = 1000.0 * 8192          # 「上采样速率」= 8.192 MHz（低速率 1 kHz × 8192）
F1, F2 = 300.0, 400.0       # 两个纯音（fs_low/2 = 500，「输出 Nyquist」）
R1, R2 = 0.5, 0.5           # 幅度
B1, B2, B3 = 1.0, 1.0, 1.0  # H(z) = b1 z + b2 z² + b3 z³
OUT_NYQ = 500.0             # 低速率 1 kHz 的 Nyquist（画出来看哪些分量能活过输出滤波）
N = 81920                   # 10 ms；f1/f2 各 3/4 个整周期；谱栅格 100 Hz，所有产物都落在栅格上


# ------------------------------------------------------------
# 理论：每条谱线的余弦幅度
# ------------------------------------------------------------

def h_complex(z):
    """解析整形器 H(z) = b1 z + b2 z² + b3 z³（只含 1/2/3 阶）。"""
    return B1 * z + B2 * z ** 2 + B3 * z ** 3


def h_real(x):
    """同一个多项式，但作用在实信号上（“实链路”对照）。"""
    return B1 * x + B2 * x ** 2 + B3 * x ** 3


def theory_complex():
    """复链路期望谱线。k1+k2 = n ≤ 3，幅度 = b_n·(多重组合数)·r1^k1 r2^k2。"""
    return {
        F1: B1 * R1,
        F2: B1 * R2,
        2 * F1: B2 * R1 ** 2,
        2 * F2: B2 * R2 ** 2,
        F1 + F2: 2 * B2 * R1 * R2,
        3 * F1: B3 * R1 ** 3,
        3 * F2: B3 * R2 ** 3,
        2 * F1 + F2: 3 * B3 * R1 ** 2 * R2,
        F1 + 2 * F2: 3 * B3 * R1 * R2 ** 2,
    }


def theory_real():
    """实链路期望谱线。x = r1cos1 + r2cos2，逐项展开：
    cos³θ = (3cosθ+cos3θ)/4，cos²θcosφ = cosφ/2 + (cos(2θ+φ)+cos(2θ−φ))/4。"""
    return {
        0.0: B2 * (R1 ** 2 + R2 ** 2) / 2,
        F1: B1 * R1 + B3 * (3 * R1 ** 3 / 4 + 3 * R1 * R2 ** 2 / 2),
        F2: B1 * R2 + B3 * (3 * R2 ** 3 / 4 + 3 * R1 ** 2 * R2 / 2),
        2 * F1: B2 * R1 ** 2 / 2,
        2 * F2: B2 * R2 ** 2 / 2,
        F1 + F2: B2 * R1 * R2,
        abs(F2 - F1): B2 * R1 * R2,
        3 * F1: B3 * R1 ** 3 / 4,
        3 * F2: B3 * R2 ** 3 / 4,
        2 * F1 + F2: 3 * B3 * R1 ** 2 * R2 / 4,
        abs(2 * F1 - F2): 3 * B3 * R1 ** 2 * R2 / 4,
        F1 + 2 * F2: 3 * B3 * R1 * R2 ** 2 / 4,
        abs(F1 - 2 * F2): 3 * B3 * R1 * R2 ** 2 / 4,
    }


SUM_IMD = {F1 + F2, 2 * F1 + F2, F1 + 2 * F2}                       # 求和互调（橙）
DIFF_IMD = {abs(F2 - F1), abs(2 * F1 - F2), abs(F1 - 2 * F2), 0.0}  # 差频/直流（红）
FUND = {F1, F2}
OWN_HARM = {2 * F1, 2 * F2, 3 * F1, 3 * F2}


def color_of(f):
    if f in DIFF_IMD:
        return "tab:red"
    if f in SUM_IMD:
        return "tab:orange"
    if f in OWN_HARM:
        return "tab:blue"
    return "0.35"


# ------------------------------------------------------------
# 测量
# ------------------------------------------------------------

def spectrum(y):
    """返回 (频率, 余弦幅度, dB re 1.0)。矩形窗 + 整周期 → 精确谱线。

    单边谱的 ×2 只适用于正频率；DC（与恰好落在 Nyquist 的 bin）不折半，不能乘 2。
    """
    y = np.asarray(y, dtype=float)
    sp = np.abs(np.fft.rfft(y)) * 2.0 / len(y)
    sp[0] *= 0.5
    fr = np.fft.rfftfreq(len(y), 1.0 / FS)
    return fr, sp, 20.0 * np.log10(np.maximum(sp, 1e-30))


def measure(y, theory):
    """在理论频率处取幅度，并找出理论之外的峰（> -80 dB）。返回 (rows, unexpected)。"""
    fr, sp, db = spectrum(y)
    rows = []
    for f, th in sorted(theory.items()):
        i = int(round(f / (FS / N)))
        rows.append((f, sp[i], 20 * np.log10(max(th, 1e-30)), 20 * np.log10(max(sp[i], 1e-30))))
    unexpected = []
    for i in range(len(fr)):
        if db[i] < -80:
            continue
        if min((abs(fr[i] - f) for f in theory), default=1e9) > FS / N * 1.5:
            unexpected.append((float(fr[i]), float(db[i])))
    return rows, unexpected


def main() -> int:
    t = np.arange(N) / FS
    z = R1 * np.exp(2j * np.pi * F1 * t) + R2 * np.exp(2j * np.pi * F2 * t)

    # 三种输出
    y_no_fund = (h_complex(z) - B1 * z).real          # 仅 2/3 阶：H = z² + z³
    y_cplx = h_complex(z).real                        # H = z + z² + z³
    y_real = h_real(z.real)                           # 实链路，同一多项式

    print(f"fs = {FS:.0f} Hz（低速率 1 kHz × 8192，直接合成、无插值）")
    print(f"两音 f1 = {F1:.0f} Hz、f2 = {F2:.0f} Hz，r1 = r2 = {R1}，"
          f"H(z) = {B1}·z + {B2}·z² + {B3}·z³\n")

    cases = [
        ("(b) 复链路  仅 2/3 阶：H = z² + z³", y_no_fund, theory_complex_wo_fund()),
        ("(c) 复链路  含基波：H = z + z² + z³", y_cplx, theory_complex()),
        ("(d) 实链路  同一多项式 P(x) = x + x² + x³", y_real, theory_real()),
    ]
    fails = 0
    for tag, y, th in cases:
        rows, unexp = measure(y, th)
        print(f"{tag}")
        print(f"   {'频率':>7} {'归属':<10} {'实测':>8} {'理论':>8}   判定")
        for f, amp, th_db, me_db in rows:
            kind = ("基波" if f in FUND else "自身谐波" if f in OWN_HARM
                    else "求和互调" if f in SUM_IMD else "差频/直流" if f in DIFF_IMD else "其他")
            ok = abs(me_db - th_db) < 0.01
            fails += not ok
            print(f"   {f:7.0f} {kind:<10} {me_db:7.1f} {th_db:7.1f}   "
                  + ("✓" if ok else f"✗ 差 {me_db - th_db:+.2f} dB"))
        fails += len(unexp)
        print(f"   理论之外的峰（> -80 dB）：{len(unexp)} 个"
              + ("" if not unexp else "  " + "，".join(f"{f:.0f} Hz {v:.1f} dB" for f, v in unexp[:6])))
        print()

    # ---- 断言：差频区必须数值为零；通道必须独立 ----
    fr, sp, db = spectrum(y_cplx)
    low = fr < 250.0
    worst = float(db[low].max())
    print(f"[检查 1] 复链路 0..250 Hz（差频区）最大分量 = {worst:.1f} dB（阈值 -200 dB）")
    fails += worst > -200.0
    dmax = max(db[int(round(f / (FS / N)))] for f in DIFF_IMD)
    print(f"[检查 2] 复链路 差频/直流处最大分量 = {dmax:.1f} dB（阈值 -200 dB）")
    fails += dmax > -200.0
    i2 = int(round(2 * F1 / (FS / N)))
    z_solo = R1 * np.exp(2j * np.pi * F1 * t)
    sp_solo = np.abs(np.fft.rfft(h_complex(z_solo).real)) * 2 / N
    dev = abs(sp[i2] - sp_solo[i2])
    print(f"[检查 3] 通道独立：2f1 幅度  双音 {sp[i2]:.6f} / 单音 {sp_solo[i2]:.6f}"
          f"  → 偏差 {dev:.2e}（阈值 1e-12）")
    fails += dev > 1e-12

    # ---- 图 ----
    fig, axes = plt.subplots(4, 1, figsize=(10.0, 10.4), sharex=True)
    panels = [
        ("(a) 解析双音输入 z = r1·e^{iω1t} + r2·e^{iω2t}（单边谱）", None, None),
        ("(b) 复链路 H = z² + z³（仅 2/3 阶，无基波）", y_no_fund, theory_complex_wo_fund()),
        ("(c) 复链路 H = z + z² + z³（基波 + 2/3 阶）", y_cplx, theory_complex()),
        ("(d) 实链路 同一多项式 P(x) = x + x² + x³ —— 差频/直流出现了", y_real, theory_real()),
    ]
    for ax, (title, y, th) in zip(axes, panels):
        if y is None:
            fr_i = np.array([F1, F2])
            db_i = np.array([20 * np.log10(R1), 20 * np.log10(R2)])
            ax.vlines(fr_i, -70, db_i, color="0.35", lw=1.6)
        else:
            fr_i, _, db_i = spectrum(y)
            loud = db_i > -70
            ax.vlines(fr_i[loud], -70, db_i[loud],
                      color=[color_of(f) for f in fr_i[loud]], lw=1.6)
            for f in th:
                ax.annotate(f"{f:.0f}", (f, 1.0), ha="center", va="bottom",
                            fontsize=6.5, color=color_of(f), rotation=90)
        ax.axvline(OUT_NYQ, color="tab:red", ls=":", lw=1.4)
        ax.axvspan(0, OUT_NYQ, color="tab:red", alpha=0.05)
        ax.set_xlim(0, 1350)
        ax.set_ylim(-70, 10)
        ax.set_ylabel("幅度 [dB]")
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].text(OUT_NYQ + 15, -62, "输出 Nyquist 500 Hz", color="tab:red", fontsize=8)
    handles = [plt.Line2D([], [], color=c, lw=2, label=lab) for c, lab in [
        ("0.35", "基波"), ("tab:blue", "各音自身谐波"), ("tab:orange", "求和互调（> 载频）"),
        ("tab:red", "差频/直流（< 载频，滤波器救不了）")]]
    axes[0].legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.9)
    fig.suptitle(f"复数整形器的双音输出：fs = {FS/1e6:.3f} MHz 直接合成，"
                 f"f1 = {F1:.0f} Hz，f2 = {F2:.0f} Hz，H(z) = {B1}z + {B2}z² + {B3}z³",
                 fontsize=10.5)
    fig.tight_layout()
    out = Path(__file__).parent / "output" / "complex_two_tone.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"\n已写出 {out}")
    print(f"断言失败项：{fails}")
    return 1 if fails else 0


def theory_complex_wo_fund():
    """H = z² + z³（去掉 b1 项）——「仅 2 阶和 3 阶谐波」的字面读法。"""
    th = theory_complex()
    th.pop(F1)
    th.pop(F2)
    return th


if __name__ == "__main__":
    raise SystemExit(main())
