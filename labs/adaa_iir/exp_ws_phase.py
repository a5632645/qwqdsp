"""AA-IIR waveshaper 的**相位 / 等效延迟**（[19] §3.2 明确提醒的一点）。

因果 IIR 核不是线性相位 → AA-IIR 的等效延迟随频率变化；而 AA-FIR-1/2 的延迟是
**恒定** 0.5 / 1 采样（矩形/三角核的质心）。对「非线性在反馈环里」的用法（wave digital
结构、带非线性反馈的滤波器），这是要单独算进延迟预算的量。

做法：线性情形 f(x)=x 时 AA-IIR 严格等价于一个线性滤波器（论文附录 B 给了一阶低通的
闭式 H(z)，`check_ws.py` 第 4 节已逐点核对）。`ws.linear_transfer` 把这一闭式推广到
任意极点组合（推导见其 docstring），直接给出精确的 H(e^{jω}) → 群延迟 = -d∠H/dω。

用法: python qwqdsp/labs/adaa_iir/exp_ws_phase.py
输出: qwqdsp/labs/adaa_iir/output/ws_group_delay.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np

import ws as W

FS = 44100.0
PROBE = (100.0, 1000.0, 5000.0, 10000.0, 20000.0)


def fir_transfer(kind: int, w):
    """AA-FIR 在线性情形下的精确响应：一阶 = (1+z⁻¹)/2，二阶 = (1+z⁻¹)²/4。"""
    z = np.exp(-1j * np.asarray(w))
    return (1.0 + z) / 2.0 if kind == 1 else ((1.0 + z) / 2.0) ** 2


def main() -> int:
    w = np.linspace(1e-6, np.pi, 200001)
    freqs = w / (2.0 * np.pi) * FS
    curves = {
        "AA-FIR-1": fir_transfer(1, w),
        "AA-FIR-2": fir_transfer(2, w),
        "AA-IIR-1": W.linear_transfer(W.aa_iir_1(), w),
        "AA-IIR-2": W.linear_transfer(W.aa_iir_2(), w),
    }
    fig, ax = plt.subplots(figsize=(9, 4.5))
    print(f"{'方法':<10}" + "".join(f"{f / 1000:>9.1f}k" for f in PROBE) + f"{'波动':>9}")
    scale = w[1] / (2.0 * np.pi / FS)
    for m, H in curves.items():
        ph = np.unwrap(np.angle(H))
        gd = -np.gradient(ph, w)            # 群延迟 [采样]
        vals = [float(np.interp(f, freqs, gd)) for f in PROBE]
        band = (freqs > 20.0) & (freqs < 20000.0)
        print(f"{m:<10}" + "".join(f"{v:>10.2f}" for v in vals)
              + f"{gd[band].max() - gd[band].min():>9.2f}")
        ax.semilogx(freqs, gd, label=m)
    del scale
    ax.set_xlim(20, 20000)
    ax.set_xlabel("频率 [Hz]")
    ax.set_ylabel("群延迟 [采样]")
    ax.set_title("等效延迟：AA-FIR 恒定（0.5 / 1 采样），AA-IIR 随频率变化")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out = Path(__file__).parent / "output" / "ws_group_delay.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150)
    print("\n直流增益（应≈1；AA-IIR-2 因丢掉直接项 A0=1e-3 而略低）：")
    for m, H in curves.items():
        print(f"  {m:<10} {abs(H[0]):.6f}  ({20 * np.log10(abs(H[0])):+.4f} dB)")
    print(f"\n已写出 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
