# -*- coding: utf-8 -*-
"""
run.py
======

把数字 IIR 低通滤波器 H(z)=B(z)/A(z) 分解为「两条全通链之和」的探索脚本：

    H = (A0 + A1)/2        （低通）
    Hc= (A0 - A1)/2        （功率互补高通，|H|^2 + |Hc|^2 = 1）

它会

1. 用 scipy 设计一个奇数阶低通（butter / cheby1 / cheby2 / ellip）；
2. 用 :mod:`tf2ca` 求出两条全通链的分母 (D0, D1)，打印每个一阶/二阶全通节点；
3. 数值校验：频响最大偏差、全通链平坦度；
4. 对照 MATLAB ``tf2ca`` 的发布样例（EMQF 文档中 ellip(7,2,40,0.3)）；
5. 画 6 张图：幅频对比、z 平面极点分配、两链相位、极点模长交替规则、
   冲激响应对比、LP/HP 互补幅频。

用法
----
    python run.py                       # 默认 ellip 7 阶, cutoff=0.3, rp=2, rs=40
    python run.py --kind butter --order 7 --cutoff 0.3
    python run.py --check               # 跑一遍四族 × 多阶数的批量校验表
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 无界面后端
import matplotlib.pyplot as plt
import numpy as np

from tf2ca import (
    allpass_flatness,
    allpass_response,
    bc_antisymmetric_sqrt,
    design,
    max_response_error,
    reconstruct,
    response,
    sections,
    split_by_radius,
    tf_response,
    to_biquads,
)

SCRIPT_DIR = Path(__file__).resolve().parent


# ------------------------------------------------------------
# 打印与校验
# ------------------------------------------------------------


def describe_chain(name: str, d: np.ndarray) -> None:
    """打印一条全通链的阶数与各一阶/二阶节点。"""
    print(f"  {name}(z): 阶数 = {d.size - 1}")
    for k, (order, poly) in enumerate(to_biquads(d)):
        num = poly[::-1]
        print(
            f"    # 节点 {k}: {order} 阶全通"
            f"  分子 = {np.array2string(num, precision=6, suppress_small=True)}"
            f"  分母 = {np.array2string(poly, precision=6, suppress_small=True)}"
        )


def report(b: np.ndarray, a: np.ndarray, d0: np.ndarray, d1: np.ndarray) -> None:
    """校验并打印分解结果。"""
    b_hat, a_hat, bc = reconstruct(d0, d1)
    err = max_response_error(b, a, d0, d1)
    flat = allpass_flatness(d0, d1)
    print(f"  频响最大偏差 |H_allpass - H_direct| = {err:.3e}")
    print(f"  全通链平坦度 max||A_i|-1|          = {flat:.3e}")
    print(f"  还原分母/分子最大偏差              = "
          f"{max(np.max(np.abs(a_hat / a_hat[0] - a / a[0])), np.max(np.abs(b_hat - b))):.3e}")
    print(f"  互补分子 B_c 系数（反镜像，c[k] = -c[n-k]）:")
    print(f"    {np.array2string(bc, precision=6, suppress_small=True)}")


def validate_reference() -> None:
    """对照 EMQF 文档里 MATLAB tf2ca 的发布样例 ellip(7,2,40,0.3)。

    参考值取自 vadkudr/EMQFfilters 的 EMQFdemo（四位小数）。
    """
    print("== 对照 MATLAB tf2ca 发布样例：ellip(7, 2, 40, 0.3) ==")
    b = np.array([0.0237, -0.0386, 0.0604, -0.0148, -0.0148, 0.0604, -0.0386, 0.0237])
    a = np.array([1.0, -4.4713, 10.0742, -14.0199, 12.9236, -7.8394, 2.9149, -0.5207])
    d0_ref = np.array([1.0, -2.5163, 3.3183, -2.2130, 0.7457])
    d1_ref = np.array([1.0, -1.9550, 1.8366, -0.6983])
    d0, d1 = split_by_radius(a)

    def rel(p, q):
        p, q = np.asarray(p, float), np.asarray(q, float)
        if p.size != q.size:
            return None
        return float(np.linalg.norm(p / np.linalg.norm(p) - q / np.linalg.norm(q)))

    pair = [rel(d0, d0_ref), rel(d1, d1_ref)]
    swap = [rel(d0, d1_ref), rel(d1, d0_ref)]
    same = (pair[0] is not None and pair[1] is not None) or (swap[0] is not None and swap[1] is not None)
    print(f"  我们的 D0 = {np.array2string(d0, precision=4)}, D1 = {np.array2string(d1, precision=4)}")
    print(f"  参考的 d0 = {d0_ref}, d1 = {d1_ref}")
    print(f"  与参考一致（允许分支互换，四位小数容差）：{same}")
    print()


def run_check() -> None:
    """批量校验：四族 × 多阶数 × 多截止频率的频响偏差。"""
    print("== 批量校验（split_by_radius）==")
    print(f"{'family':8s} {'order':>5s} {'maxErr(freq)':>14s} {'allpassFlat':>13s}")
    worst = 0.0
    for kind in ("butter", "cheby1", "cheby2", "ellip"):
        for order in (3, 5, 7, 9, 11):
            for cutoff in (0.1, 0.2, 0.35, 0.5, 0.7):
                b, a = design(kind, order, cutoff)
                d0, d1 = split_by_radius(a)
                err = max_response_error(b, a, d0, d1)
                flat = allpass_flatness(d0, d1)
                worst = max(worst, err)
                flag = "" if err < 1e-4 else "   <- 注意"
                print(f"{kind:8s} {order:5d} {err:14.3e} {flat:13.3e}{flag}")
    print(f"  最差频响偏差 = {worst:.3e}（窄带高阶受极点提取的条件数限制）")


# ------------------------------------------------------------
# 绘图
# ------------------------------------------------------------


def plot(b, a, d0, d1, out: Path) -> None:
    """画六张图说明分解的正确性与极点交替规则。"""
    w = np.linspace(1e-6, np.pi - 1e-6, 4000)
    ref = tf_response(b, a, w)
    lp, hp = response(d0, d1, w)
    a0 = allpass_response(d0, w)
    a1 = allpass_response(d1, w)
    bc = bc_antisymmetric_sqrt(b, a)

    fig, ax = plt.subplots(2, 3, figsize=(16, 8.5))

    # (0,0) 幅频：直接型 vs 并行全通
    ax[0, 0].plot(w / np.pi, 20 * np.log10(np.abs(ref) + 1e-18), "k", lw=2, label="direct B/A")
    ax[0, 0].plot(w / np.pi, 20 * np.log10(np.abs(lp) + 1e-18), "r--", lw=1, label="(A0+A1)/2")
    ax[0, 0].set(title="mag: direct-form vs parallel-allpass", xlabel="$\\omega/\\pi$",
                 ylabel="dB", ylim=(-80, 3))
    ax[0, 0].grid(alpha=0.3)
    ax[0, 0].legend()

    # (0,1) z 平面极点/零点
    theta = np.linspace(0, 2 * np.pi, 400)
    ax[0, 1].plot(np.cos(theta), np.sin(theta), "k:", lw=0.8)
    for d, color, lab in ((d0, "#1f77b4", "D0 poles"), (d1, "#d62728", "D1 poles")):
        p = np.roots(d)
        ax[0, 1].plot(p.real, p.imag, "x", color=color, ms=9, label=lab)
    z = np.roots(b)
    ax[0, 1].plot(z.real, z.imag, "o", mfc="none", mec="g", ms=7, label="B zeros")
    ax[0, 1].axhline(0, color="gray", lw=0.5)
    ax[0, 1].axvline(0, color="gray", lw=0.5)
    ax[0, 1].set(title="z-plane: pole assignment", xlabel="Re", ylabel="Im", aspect="equal")
    ax[0, 1].grid(alpha=0.3)
    ax[0, 1].legend(loc="upper left", fontsize=8)

    # (0,2) 两链相位 & |H| = |cos(dphi/2)|
    ax[0, 2].plot(w / np.pi, np.unwrap(np.angle(a0)), label="$\\arg A_0$")
    ax[0, 2].plot(w / np.pi, np.unwrap(np.angle(a1)), label="$\\arg A_1$")
    ax[0, 2].plot(w / np.pi, np.abs(np.cos((np.unwrap(np.angle(a0)) - np.unwrap(np.angle(a1))) / 2)),
                  "k--", label="$|\\cos(\\Delta\\varphi/2)|$")
    ax[0, 2].set(title="allpass phases", xlabel="$\\omega/\\pi$", ylabel="rad")
    ax[0, 2].grid(alpha=0.3)
    ax[0, 2].legend(fontsize=8)

    # (1,0) 极点模长排序 + 交替分配
    nodes = sections(a)
    nodes.sort(key=lambda s: abs(s["poles"][0]))
    radii = [abs(n["poles"][0]) for n in nodes]
    for k, r in enumerate(radii):
        ax[1, 0].stem([k], [r], linefmt=("#1f77b4" if k % 2 == 0 else "#d62728"))
        ax[1, 0].plot(k, r, "o", color=("#1f77b4" if k % 2 == 0 else "#d62728"))
    ax[1, 0].set(title="pole |p| sorted -> alternate branches", xlabel="rank (|p| asc)",
                 ylabel="|pole|")
    ax[1, 0].grid(alpha=0.3)

    # (1,1) 冲激响应对比
    n_imp = 200
    x = np.zeros(n_imp)
    x[0] = 1.0
    from scipy.signal import lfilter

    y_direct = lfilter(b, a, x)
    y_lp = lfilter(d0[::-1], d0, x) / 2 + lfilter(d1[::-1], d1, x) / 2
    ax[1, 1].plot(y_direct, "k", lw=2, label="direct")
    ax[1, 1].plot(y_lp, "r--", label="allpass sum")
    ax[1, 1].set(title=f"impulse response (max|Δ|={np.max(np.abs(y_direct - y_lp)):.1e})",
                 xlabel="n")
    ax[1, 1].grid(alpha=0.3)
    ax[1, 1].legend(fontsize=8)

    # (1,2) 互补 LP/HP
    ax[1, 2].plot(w / np.pi, 20 * np.log10(np.abs(lp) + 1e-18), label="LP=(A0+A1)/2")
    ax[1, 2].plot(w / np.pi, 20 * np.log10(np.abs(hp) + 1e-18), label="HP=(A0-A1)/2")
    ax[1, 2].plot(w / np.pi, 20 * np.log10(np.abs(lp) ** 2 + np.abs(hp) ** 2), "k:", lw=1,
                  label="$|H|^2+|H_c|^2$")
    ax[1, 2].set(title="power-complementary LP/HP", xlabel="$\\omega/\\pi$", ylabel="dB",
                 ylim=(-80, 3))
    ax[1, 2].grid(alpha=0.3)
    ax[1, 2].legend(fontsize=8)

    fig.suptitle("tf2ca: IIR lowpass = 1/2 (A0 + A1)", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"图已保存：{out}")


# ------------------------------------------------------------
# main
# ------------------------------------------------------------


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="IIR lowpass -> sum of two allpass chains (tf2ca)")
    p.add_argument("--kind", default="ellip", choices=["butter", "cheby1", "cheby2", "ellip"])
    p.add_argument("--order", type=int, default=7, help="奇数阶")
    p.add_argument("--cutoff", type=float, default=0.3, help="Nyquist 归一截止频率 (0..1)")
    p.add_argument("--ripple", type=float, default=2.0, help="通带波纹 dB (cheby1/ellip)")
    p.add_argument("--stop", type=float, default=40.0, help="阻带衰减 dB (cheby2/ellip)")
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--check", action="store_true", help="跑批量校验表")
    args = p.parse_args(argv)

    if args.order % 2 == 0:
        p.error("只能分解奇数阶：--order 必须为奇数")

    if args.check:
        run_check()
        return

    out = args.out or SCRIPT_DIR / "output" / f"{args.kind}_{args.order}_{args.cutoff:g}.png"

    print(f"== 设计 {args.kind} 低通：order={args.order}, cutoff={args.cutoff} ==")
    b, a = design(args.kind, args.order, args.cutoff, args.ripple, args.stop)
    print(f"  B = {np.array2string(b, precision=6, suppress_small=True)}")
    print(f"  A = {np.array2string(a, precision=6, suppress_small=True)}")
    print()

    print("== 极点分组（按 |p| 升序交替）==")
    d0, d1 = split_by_radius(a)
    describe_chain("D0", d0)
    describe_chain("D1", d1)
    print()

    print("== 校验 ==")
    report(b, a, d0, d1)
    print()

    validate_reference()
    plot(b, a, d0, d1, out)


if __name__ == "__main__":
    main()
