"""实验：非线性的「查表 / 分段线性拟合（CPWL）」能否喂给 ADAA-IIR？代价与精度如何？

问题背景
--------
[DAFx25] "Simplifying Antiderivative Antialiasing with Lookup Table Integration"
证明 **FIR 版 ADAA（ADAA-LUT）** 可以用数值积分出来的表代替解析原函数。但 **IIR 版**
的积分核 ``I_n = ∫₀¹ f(x_n+uΔ)e^{β(1-u)}du``（``γ=β/Δ``）依赖逐样本差分 Δ，
不能预计算与 Δ 无关的表，所以「IIR 能不能 LUT 化」是个真问题。

``ws.py`` 给出的答案：只要 f **分段线性**，每段积分有闭式（``mean_integral_pwl``，
按「采样区间跨过的段数」循环求和），于是用 k 段 CPWL 拟合任意 f 即可 LUT 化 IIR。

本脚本量化三件事：
1. 拟合大小 K 扫描下的 SNR（AA-IIR-1/2 与 FIR 对照 AA-FIR-1/2）；
2. 无拟合误差的「解析/精确」参照 ``aa_iir_fn``（Gauss-Legendre 数值求积），
   据此指出 K 多大以后 SNR 不再提升（拟合误差不再主导）；
3. 代价：每采样平均跨段数 + 运行耗时。

自检：K=512 的 AA-IIR 波形与数值求积参照的最大偏差应 < 1e-3（退出码反映）。

用法: python qwqdsp/labs/adaa_iir/exp_ws_tanh.py
输出: qwqdsp/labs/adaa_iir/output/exp_ws_tanh.png
      qwqdsp/labs/adaa_iir/output/exp_ws_tanh.csv
      stdout: SNR 表 + 自检 + 代价
"""
from __future__ import annotations

import csv
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np

import aaiir as A
import ws as W

# ------------------------------------------------------------
# 实验配置（[DAFx25] 的 tanh 参数）
# ------------------------------------------------------------
FS = 44100.0            # 采样率 [Hz]
F0 = 2000.0             # 输入正弦频率 [Hz]
AMP = 3.0               # 输入幅度（过 tanh 后明显削波）
DUR = 1.5               # 时长 [s]（> 1 s）
N = int(round(FS * DUR))
ALPHA, BETA = 1.0, 0.3  # f(x) = alpha * tanh(x / beta)

FIT_LO, FIT_HI = -4.0, 4.0      # CPWL 拟合区间（覆盖输入 ±3 并留余量）
KS = [4, 8, 16, 32, 64, 128, 512]
QUAD_ORDER = 8                  # aa_iir_fn 的 Gauss-Legendre 阶数
QUAD_ORDER_CHK = 16             # 参照自身的收敛性检查
CHECK_K = 512                   # 自检用的 K
CHECK_TOL = 1e-3                # 自检阈：波形最大偏差
WARMUP = 4096                   # 自检/绘图的启动瞬态
REPEAT = 3                      # 计时重复次数（取最快）


def f_exact(x):
    """目标非线性 f(x) = alpha·tanh(x/beta)。"""
    return ALPHA * np.tanh(np.asarray(x, dtype=float) / BETA)


# ------------------------------------------------------------
# 每采样跨段数（CPWL 积分的代价）
# ------------------------------------------------------------

def segments_crossed(x, f: W.Pwl):
    """每个采样区间 [x_n, x_{n+1}] 跨过的拟合断点数（= 需循环求和的槽位数 - 1）。"""
    x = np.asarray(x, dtype=float)
    lo = np.minimum(x[:-1], x[1:])
    hi = np.maximum(x[:-1], x[1:])
    return (f.seg_index(hi) - f.seg_index(lo)).astype(np.int64)


def timeit(fn, repeat: int = REPEAT):
    """跑 repeat 次取最快耗时，返回 (输出, 秒)。"""
    best, out = float("inf"), None
    for _ in range(repeat):
        t0 = time.perf_counter()
        out = fn()
        best = min(best, time.perf_counter() - t0)
    return out, best


# ------------------------------------------------------------
# 主流程
# ------------------------------------------------------------

def main() -> int:
    t = np.arange(N) / FS
    x = AMP * np.sin(2.0 * np.pi * F0 * t)          # 输入正弦（相位对 2π 非整周期，加窗处理）

    iir1, iir2 = W.aa_iir_1(), W.aa_iir_2()

    print(f"配置：fs={FS:.0f} Hz, f0={F0:.0f} Hz, 幅度={AMP}, 时长={DUR} s（N={N}）")
    print(f"      非线性 f(x) = {ALPHA}·tanh(x/{BETA})，CPWL 拟合区间 [{FIT_LO}, {FIT_HI}]")
    print(f"      输入范围 [{x.min():.3f}, {x.max():.3f}]")

    # ---------- 精确参照：Gauss-Legendre 数值求积 ----------
    print("\n[参照] aa_iir_fn(..., fn=f, order=%d)（无数值拟合误差）" % QUAD_ORDER)
    ref = {}
    for name, filt in (("AA-IIR-1", iir1), ("AA-IIR-2", iir2)):
        y, dt = timeit(lambda flt=filt: W.aa_iir_fn(x, flt, f_exact, order=QUAD_ORDER))
        ref[name] = {"y": y, "snr": A.snr_db(y, FS, F0, warmup=WARMUP), "sec": dt}
        print(f"  {name}: SNR = {ref[name]['snr']:7.2f} dB   耗时 {dt * 1e3:7.2f} ms")

    # 参照自身收敛性：order 8 vs 16（确认「精确」名副其实）
    conv = {}
    for name, filt in (("AA-IIR-1", iir1), ("AA-IIR-2", iir2)):
        y16 = W.aa_iir_fn(x, filt, f_exact, order=QUAD_ORDER_CHK)
        conv[name] = float(np.max(np.abs(ref[name]["y"] - y16)))
    print(f"  收敛性检查 order {QUAD_ORDER} vs {QUAD_ORDER_CHK}：波形最大偏差 "
          f"AA-IIR-1 = {conv['AA-IIR-1']:.2e}, AA-IIR-2 = {conv['AA-IIR-2']:.2e}")

    # ---------- K 扫描 ----------
    rows = []
    sig = {name: [] for name in ("AA-IIR-1", "AA-IIR-2", "AA-FIR-1", "AA-FIR-2")}
    cross, secs = [], {name: [] for name in sig}
    fit_err = []

    for K in KS:
        pwl = W.fit_from_callable(f_exact, FIT_LO, FIT_HI, K, clamp=True)
        err = float(np.max(np.abs(pwl.value(x) - f_exact(x))))   # 拟合值本身的偏差
        fit_err.append(err)

        y11, t11 = timeit(lambda p=pwl: W.aa_iir_pwl(x, iir1, p))
        y12, t12 = timeit(lambda p=pwl: W.aa_iir_pwl(x, iir2, p))
        y21, t21 = timeit(lambda p=pwl: W.aa_fir1(x, p))
        y22, t22 = timeit(lambda p=pwl: W.aa_fir2(x, p))

        snr = {
            "AA-IIR-1": A.snr_db(y11, FS, F0, warmup=WARMUP),
            "AA-IIR-2": A.snr_db(y12, FS, F0, warmup=WARMUP),
            "AA-FIR-1": A.snr_db(y21, FS, F0, warmup=WARMUP),
            "AA-FIR-2": A.snr_db(y22, FS, F0, warmup=WARMUP),
        }
        for k in sig:
            sig[k].append(snr[k])
        secs["AA-IIR-1"].append(t11)
        secs["AA-IIR-2"].append(t12)
        secs["AA-FIR-1"].append(t21)
        secs["AA-FIR-2"].append(t22)

        c = segments_crossed(x, pwl)
        cross.append((float(c.mean()), int(c.max())))
        rows.append((K, snr, float(c.mean()), int(c.max()), err, t11, t12, t21, t22))

    # ---------- 输出 SNR 表 ----------
    print("\nSNR [dB] 表（Δ = 与精确参照之差；K 小时拟合出的 f 本身已明显偏离，"
          "其 SNR 不可与参照直接比较）")
    hdr = (f"{'K':>5} | {'AA-IIR-1':>9} {'AA-IIR-2':>9} {'AA-FIR-1':>9} {'AA-FIR-2':>9}"
           f" | {'ΔIIR-1':>7} {'ΔIIR-2':>7} | {'跨段/采样':>10} {'最大':>5} | {'拟合偏差':>9}")
    print(hdr)
    print("-" * len(hdr))
    for K, snr, cm, cx, err, *_ in rows:
        print(f"{K:>5} | {snr['AA-IIR-1']:>9.2f} {snr['AA-IIR-2']:>9.2f}"
              f" {snr['AA-FIR-1']:>9.2f} {snr['AA-FIR-2']:>9.2f}"
              f" | {snr['AA-IIR-1'] - ref['AA-IIR-1']['snr']:>+7.2f}"
              f" {snr['AA-IIR-2'] - ref['AA-IIR-2']['snr']:>+7.2f}"
              f" | {cm:>10.3f} {cx:>5d} | {err:>9.2e}")
    print(f"{'精确':>5} | {ref['AA-IIR-1']['snr']:>9.2f} {ref['AA-IIR-2']['snr']:>9.2f}"
          f" {'—':>9} {'—':>9} | {'0.00':>7} {'0.00':>7} | {'—':>10} {'—':>5} | {'—':>9}")
    print("注：K=4/8 处 SNR 反而高于精确参照，并非「更好」——拟合偏差高达 "
          "0.47~0.68，拟合出的 f 是远比真 tanh 平缓的另一条曲线（谐波更少、更少混叠）。"
          "有意义的比较只看 K ≥ 16 起的收敛段。")

    # ---------- K 饱和点 ----------
    # 判据不能用「首次 |Δ|<1dB」：K 太小时 CPWL 拟合出的 f 本身是另一条（远更平缓的）
    # 非线性，谐波更少、SNR 反而虚高（见上表拟合偏差 0.47~0.68）。故取「从此往后
    # 全部 K 都在参照 ±1 dB 内」的最小 K，即 SNR 从上方收敛到方法自身下限的拐点。
    print("\n拟合误差不再主导的 K（自该 K 起，所有 K' ≥ K 均有 |SNR - 精确参照| < 1 dB）")
    SAT_TOL = 1.0
    satK = {}
    for name in ("AA-IIR-1", "AA-IIR-2"):
        r = ref[name]["snr"]
        sat = next((K for i, K in enumerate(KS)
                    if all(abs(v - r) < SAT_TOL for v in sig[name][i:])), None)
        satK[name] = sat
        print(f"  {name}: 精确参照 {r:.2f} dB；K = {sat}"
              f"（K={CHECK_K} 时 {sig[name][-1]:.2f} dB，残差 {sig[name][-1] - r:+.2f} dB）")
    for name in ("AA-FIR-1", "AA-FIR-2"):
        print(f"  {name}: 无解析参照；K={KS[-2]}→{CHECK_K} 变化 "
              f"{sig[name][-1] - sig[name][-2]:+.2f} dB"
              f"（已进入平台，残留为该方法自身的混叠下限）")
    print("  => 说明：拟合误差只要求中等 K（几十段量级）即可压到 <1 dB；"
          "再加大 K 只是逼近方法自身下限，且代价线性增长。")

    # ---------- 自检 ----------
    print(f"\n自检：K={CHECK_K} 的 AA-IIR 波形 vs 数值求积参照（阈 {CHECK_TOL:g}）")
    ok = True
    dev = {}
    tail = slice(WARMUP, None)
    for name, filt in (("AA-IIR-1", iir1), ("AA-IIR-2", iir2)):
        pwl = W.fit_from_callable(f_exact, FIT_LO, FIT_HI, CHECK_K, clamp=True)
        y = W.aa_iir_pwl(x, filt, pwl)
        d_all = float(np.max(np.abs(y - ref[name]["y"])))
        d_tail = float(np.max(np.abs(y[tail] - ref[name]["y"][tail])))
        dev[name] = d_tail
        good = d_tail < CHECK_TOL
        ok &= good
        print(f"  {name}: 全段最大偏差 {d_all:.3e}，去瞬态后 {d_tail:.3e}"
              f"  -> {'通过' if good else '未通过'}")

    # ---------- 代价 ----------
    print("\n代价：耗时（N 采样，取最快一次）与每采样平均跨段数")
    print(f"{'K':>5} | {'IIR-1 ms':>9} {'IIR-2 ms':>9} {'FIR-1 ms':>9} {'FIR-2 ms':>9}"
          f" | {'跨段/采样':>10} {'跨段总计':>9}")
    for K, snr, cm, cx, err, t11, t12, t21, t22 in rows:
        print(f"{K:>5} | {t11 * 1e3:>9.2f} {t12 * 1e3:>9.2f} {t21 * 1e3:>9.2f} {t22 * 1e3:>9.2f}"
              f" | {cm:>10.3f} {cm * N:>9.0f}")
    print(f"参照（求积 order={QUAD_ORDER}）：IIR-1 {ref['AA-IIR-1']['sec'] * 1e3:.2f} ms, "
          f"IIR-2 {ref['AA-IIR-2']['sec'] * 1e3:.2f} ms（与 K 无关）")
    kbig = rows[-1]
    budget = 1e6 / FS
    print(f"  · IIR 的 CPWL 代价 ≈ O(采样数 × 每区间**最大**跨段数)："
          f"K={CHECK_K} 时平均跨 {kbig[2]:.1f} 段、最大跨 {kbig[3]} 段，"
          f"``mean_integral_pwl``（ws.py，未改）的 Python 循环按最大值转 {kbig[3]} 遍"
          f"全长数组，故实测耗时随 K 近似线性（IIR-1: {rows[0][5] * 1e3:.0f} → "
          f"{kbig[5] * 1e3:.0f} ms，即每遍 ≈ "
          f"{(kbig[5] - rows[0][5]) / (kbig[3] - rows[0][3]) * 1e3:.1f} ms）。"
          f"若改为只遍历实际用到的槽位（平均 {kbig[2]:.1f} 段），预计可降到约 "
          f"{kbig[5] * kbig[2] / kbig[3] * 1e3:.0f} ms。")
    print(f"  · 实时预算 1/fs = {budget:.2f} µs/采样：K={CHECK_K} 的 IIR-1 每采样 "
          f"{kbig[5] / N * 1e6:.2f} µs = 预算的 {kbig[5] / N * 1e6 / budget * 100:.0f}%"
          f"（单实例勉强实时）；K={KS[3]} 时 {rows[3][5] / N * 1e6:.2f} µs"
          f"（{rows[3][5] / N * 1e6 / budget * 100:.0f}%）。")
    print(f"  · 对照：Gauss-{QUAD_ORDER} 数值求积参照每采样 "
          f"{ref['AA-IIR-1']['sec'] / N * 1e6:.2f} µs（IIR-1）/ "
          f"{ref['AA-IIR-2']['sec'] / N * 1e6:.2f} µs（IIR-2），"
          f"与 K 无关且比 K={CHECK_K} 的 CPWL 快 "
          f"{kbig[5] / ref['AA-IIR-1']['sec']:.0f}×。"
          f"即：只有当 K 小（≲32，此时每采样 "
          f"{rows[3][5] / N * 1e6:.2f} µs）CPWL 路线才与求积/FIR 同量级。")
    print(f"  · FIR（ADAA-LUT）代价与 K 无关（{min(r[7] for r in rows) * 1e3:.0f}~"
          f"{max(r[7] for r in rows) * 1e3:.0f} ms），因为只需一次查表读 F₁/F₂。")

    # ---------- 结论 ----------
    print("\n结论")
    print(f"  1) 能 LUT 化：分段线性 f 让 AA-IIR 的区间积分有闭式，k 段 CPWL 的 SNR"
          f"从上方单调收敛到数值求积参照；K≥{satK['AA-IIR-1']}（IIR-1）/ "
          f"K≥{satK['AA-IIR-2']}（IIR-2）后拟合误差 < 1 dB，K={CHECK_K} 时残差 "
          f"{sig['AA-IIR-1'][-1] - ref['AA-IIR-1']['snr']:+.2f} / "
          f"{sig['AA-IIR-2'][-1] - ref['AA-IIR-2']['snr']:+.2f} dB，"
          f"波形偏差 {dev['AA-IIR-1']:.1e} / {dev['AA-IIR-2']:.1e}。")
    print(f"  2) 精度只需「几十段」：K=32 已把拟合误差压到 "
          f"{abs(sig['AA-IIR-1'][3] - ref['AA-IIR-1']['snr']):.2f}/"
          f"{abs(sig['AA-IIR-2'][3] - ref['AA-IIR-2']['snr']):.2f} dB；"
          f"再加 K 只是逼近方法自身混叠下限（K=128→512 仅 "
          f"{sig['AA-IIR-1'][-1] - sig['AA-IIR-1'][-2]:+.2f} dB）。")
    print(f"  3) 代价随 K 线性增长（ws.py 的槽位循环按最大跨段数 {rows[-1][3]} 转全长数组）："
          f"IIR-1 从 K=4 的 {rows[0][5] * 1e3:.0f} ms 涨到 K=512 的 "
          f"{rows[-1][5] * 1e3:.0f} ms（每采样 {rows[-1][5] / N * 1e6:.1f} µs，"
          f"占实时预算 {rows[-1][5] / N * 1e6 / (1e6 / FS) * 100:.0f}%）；"
          f"同样的 AA-IIR-1 用 Gauss-8 求积只需 "
          f"{ref['AA-IIR-1']['sec'] / N * 1e6:.2f} µs/采样。")
    print(f"  4) FIR 版（ADAA-LUT）代价与 K 无关（{min(r[7] for r in rows) * 1e3:.0f}~"
          f"{max(r[7] for r in rows) * 1e3:.0f} ms 全程），这是 LUT 化的天然优势；"
          f"但 FIR-2 的 SNR 上限（K=512 时 {sig['AA-FIR-2'][-1]:.2f} dB）明显低于 "
          f"IIR-2（{ref['AA-IIR-2']['snr']:.2f} dB）。")
    print("  5) 注意点：本实验用等距 CPWL 拟合到 [-4,4]，而 f 在 |x|≳1 已饱和，"
          "两端大量段位浪费；对这类「中心陡、两端平」的 f，非等距分段（按曲率布点）"
          "可用更少的段达到同样精度。")

    # ---------- 写 CSV / 图 ----------
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    with open(out / "exp_ws_tanh.csv", "w", newline="", encoding="utf-8") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["K", "AA-IIR-1", "AA-IIR-2", "AA-FIR-1", "AA-FIR-2",
                      "seg_per_sample", "seg_max", "fit_dev", "iir1_ms", "iir2_ms",
                      "fir1_ms", "fir2_ms"])
        for K, snr, cm, cx, err, t11, t12, t21, t22 in rows:
            wtr.writerow([K] + [f"{snr[k]:.3f}" for k in
                                ("AA-IIR-1", "AA-IIR-2", "AA-FIR-1", "AA-FIR-2")]
                         + [f"{cm:.4f}", cx, f"{err:.3e}",
                            f"{t11 * 1e3:.3f}", f"{t12 * 1e3:.3f}",
                            f"{t21 * 1e3:.3f}", f"{t22 * 1e3:.3f}"])

    style = {
        "AA-IIR-1": dict(color="tab:red", ls="-", marker="o", ms=4, lw=1.8),
        "AA-IIR-2": dict(color="tab:blue", ls="-", marker="s", ms=4, lw=1.8),
        "AA-FIR-1": dict(color="tab:red", ls="--", marker="^", ms=4, lw=1.4),
        "AA-FIR-2": dict(color="tab:blue", ls="--", marker="v", ms=4, lw=1.4),
    }
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for name in sig:
        ax.semilogx(KS, sig[name], label=name, **style[name])
    for name, c in (("AA-IIR-1", "tab:red"), ("AA-IIR-2", "tab:blue")):
        ax.axhline(ref[name]["snr"], color=c, ls=":", lw=1.2, alpha=0.7,
                   label=f"{name} 精确参照（数值求积）")
    ksat = max(v for v in satK.values() if v)
    ax.axvline(ksat, color="0.45", ls="-.", lw=1.2,
               label=f"拟合误差 <1 dB 门限（两档取较严者 K={ksat}）")
    ax.set_xlabel("CPWL 拟合段数 K（查表大小）")
    ax.set_ylabel("SNR [dB]")
    ax.set_title(f"ADAA-LUT 用于 tanh(β={BETA}, α={ALPHA})：SNR vs 查表大小 K\n"
                 f"fs={FS:.0f} Hz, f0={F0:.0f} Hz, 幅度 {AMP}, N={N}")
    ax.set_xticks(KS)
    ax.set_xticklabels([str(k) for k in KS])
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    png = out / "exp_ws_tanh.png"
    fig.savefig(png, dpi=150)

    print(f"\n已写出 {png}\n已写出 {out / 'exp_ws_tanh.csv'}")
    if not ok:
        print("\n自检未通过：K=%d 的 IIR 波形偏差超过 %g，拟合或调用有误。"
              % (CHECK_K, CHECK_TOL))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
