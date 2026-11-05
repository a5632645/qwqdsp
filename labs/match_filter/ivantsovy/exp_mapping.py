"""频率映射：数字频率 -> 模拟频率（fs = 2 Hz，Nyquist = 1）。

这套映射是**精确**的：对论文因子规则构造出的任何原型，都有

    |Hd(e^{j theta})| = |Ha(j rho(theta))|,

其中 ``rho`` 只由设计的 ``(omega, sigma)`` 决定：

    rho(theta) = (omega/pi) t / sqrt(1 + (sigma/pi)^2 t^2),   t = tan(theta/2)

推导：一次因子 ``1 + a S`` 的数字幅度为 ``sqrt(cos^2 + kappa_a^2 sin^2)``，
``kappa_a = sqrt(omega^2 a^2 + sigma^2)/pi``；提出
``M0 = sqrt(cos^2(θ/2) + (sigma/pi)^2 sin^2(θ/2))`` 后正好是
``M0 * sqrt(1 + a^2 rho^2)``，即 ``M0`` 乘以该因子在 ``S = j rho`` 处的模拟幅度。
``M0`` 在分子分母上出现同样多次（亏欠因子补齐）而相消；原点零点
``(omega/pi) sin(θ/2)`` 也等于 ``M0 * rho``，同样吻合。

``sigma = 0`` 时退化为双线性 ``rho = (omega/pi) tan(theta/2)``（在 Nyquist 发散）；
``sigma > 0`` 时 ``rho(pi) = omega/sigma`` 有限。双点解 ``sigma = 2 omega / u`` 给出
``rho(pi) = u/2``，即 Nyquist 恰好映到 ``fs/2``。

映射被三点唯一钉死。归一化形式只有一个形状参数：

    y / (fs/2) = sin( arctan( A tan(theta/2) ) )

它自动通过 (0, 0) 与 (Nyquist, Nyquist)；Nyquist 条件 ``sigma = 2 omega / u`` 把
``A`` 写成 ``sigma/pi``，而"通过 (fc, fc)"再把它唯一确定：

    A = x_c / ( tan(pi x_c / 2) * sqrt(1 - x_c^2) ),    x_c = 2 fc / fs

（与 ``sigma/pi`` 逐点吻合到机器精度，见 :func:`near_nyquist` 的输出。）所以三点
同时匹配时映射没有剩余自由度；``x_c -> 1`` 时 ``A ≈ (pi/2) sqrt((1-x_c)/2) -> 0``，
低中频段被压平是唯一解，不是参数没调好。

``A = sigma/pi`` 的另一种写法：

    y / (fs/2) = sin( arctan( A tan(theta/2) ) ),   A = sigma/pi

低频档 ``sigma -> 2``，``A -> 2/pi = 0.637``（曲线固定、与 fc 无关）；``fc -> fs/2``
时 ``sigma -> 0``，``A -> 0``，映射被压平到 ``y ≈ A·tan(θ/2)·(fs/2)``——这就是
"高频档看起来严重畸变"的来源：不是增益错，是频率被压缩。

x / y 都取 Hz（fs = 2），因此两轴都以 Nyquist = 1 归一。

用法: python qwqdsp/labs/match_filter/ivantsovy/exp_mapping.py
产物: output/two_point_mapping.png、output/two_point_mapping_near_nyquist.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import ivantsovy as I
import two_point as TP

FS = 2.0
POLES = np.array([-1.0])
LOWPASS_ZEROS = np.array([])
OUT = Path(__file__).parent / "output"
REF_FS = 48000.0
"""截止频率按 48 kHz 下的数值给出，再按比例换算到 fs = 2 Hz。

低频段沿用 5 点测试里的 20 / 117.6 / 691.4 Hz，4 kHz 以上加密到 8 个点——
映射的形状基本都在这段里变化（u 从 12 掉到 2）。
"""

REF_CUTOFFS = ((20.0, 117.5904, 691.3754, 2000.0)
               + tuple(float(v) for v in np.geomspace(4065.0, 23900.0, 8)))
CASES = tuple((ref, FS * ref / REF_FS) for ref in REF_CUTOFFS)

NEAR_CUTOFFS = tuple(float(v) for v in np.geomspace(0.65, 0.999, 6))
"""贴 Nyquist 的 6 档截止频率（单位 Hz，fs = 2；1.0 本身无解，取 0.999）。"""

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False


def rho(theta, omega: float, sigma: float) -> np.ndarray:
    """解析映射 ``rho(theta)``。"""
    t = np.tan(np.asarray(theta, dtype=float) / 2.0)
    return (omega / np.pi) * t / np.sqrt(1.0 + (sigma / np.pi) ** 2 * t ** 2)


def numeric_mapping(fc: float, fs: float = FS, n: int = 20001):
    """用一阶低通的幅度反解出的映射（用于核对解析式）。"""
    design = TP.lowpass(fc, fs)
    x = np.linspace(0.0, fs / 2.0, n)
    hd = np.clip(I.digital_magnitude(design, x, fs), 1e-12, 1.0 - 1e-15)
    return x, fc * np.sqrt(np.maximum(1.0 / hd ** 2 - 1.0, 0.0))


def bilinear(fc: float, fs: float = FS, n: int = 20001):
    design = I.design(LOWPASS_ZEROS, POLES, fs, fc, 0.0)
    x = np.linspace(0.0, fs / 2.0, n)
    hd = np.clip(I.digital_magnitude(design, x, fs), 1e-12, 1.0 - 1e-15)
    return x, fc * np.sqrt(np.maximum(1.0 / hd ** 2 - 1.0, 0.0))


def near_nyquist() -> None:
    """fc ∈ [0.65, 1.0)、x ∈ [0.65, 1.0] 的放大图。"""
    n = len(NEAR_CUTOFFS)
    colors = [plt.get_cmap("viridis")(v) for v in np.linspace(0.08, 0.95, n)]
    listed = matplotlib.colors.ListedColormap(colors)
    norm = matplotlib.colors.BoundaryNorm(np.arange(n + 1), n)
    scalar = plt.cm.ScalarMappable(cmap=listed, norm=norm)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))
    x = np.linspace(0.65, FS / 2.0, 4001)

    ax = axes[0]
    ax.plot([0.65, 1.0], [0.65, 1.0], color="0.4", lw=1.0, ls="--", label="恒等映射 y = x")
    for index, fc in enumerate(NEAR_CUTOFFS):
        omega, sigma = TP.omega_sigma(fc, FS)
        ax.plot(x, fc * rho(np.pi * x / FS * 2.0, omega, sigma), color=colors[index],
                lw=1.9)
    xb, yb = bilinear(NEAR_CUTOFFS[0], FS)
    keep = yb <= 1.02
    ax.plot(xb[keep], yb[keep], color="0.55", lw=1.2, ls=":",
            label="双线性 σ=0（超出 1.02 已截断）")
    ax.set_xlim(0.65, 1.0)
    ax.set_ylim(0.0, 1.02)
    ax.set_xlabel("数字频率 [Hz]（fs = 2，Nyquist = 1）")
    ax.set_ylabel("等效模拟频率 [Hz]")
    ax.set_title("贴 Nyquist 的映射：y = fc·ρ(θ)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)
    labels = []
    for fc in NEAR_CUTOFFS:
        _, sigma = TP.omega_sigma(fc, FS)
        labels.append(f"{fc:.3f}  A={sigma / np.pi:.3f}")
    cb = fig.colorbar(scalar, ax=ax, ticks=np.arange(n) + 0.5, pad=0.015)
    cb.ax.set_yticklabels(labels, fontsize=7.5)
    cb.set_label("截止频率 [Hz] 与 A = σ/π")

    ax = axes[1]
    for index, fc in enumerate(NEAR_CUTOFFS):
        omega, sigma = TP.omega_sigma(fc, FS)
        theta = np.pi * x / FS * 2.0
        slope = np.gradient(fc * rho(theta, omega, sigma), x)
        ax.plot(x, np.minimum(slope, 6.0), color=colors[index], lw=1.7)
    ax.axhline(1.0, color="0.4", lw=1.0, ls="--", label="恒等映射 dy/dx = 1")
    ax.set_xlim(0.65, 1.0)
    ax.set_ylim(0.0, 6.0)
    ax.set_xlabel("数字频率 [Hz]（fs = 2，Nyquist = 1）")
    ax.set_ylabel("局部斜率 dy/dx")
    ax.set_title("局部斜率（0.65–1.0 段，1 = 恒等）")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)
    cb = fig.colorbar(scalar, ax=ax, ticks=np.arange(n) + 0.5, pad=0.015)
    cb.ax.set_yticklabels(labels, fontsize=7.5)
    cb.set_label("截止频率 [Hz] 与 A = σ/π")

    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    path = OUT / "two_point_mapping_near_nyquist.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"图：{path}")

    print("A = σ/π 与三点钉死公式 x_c/(tan(π x_c/2)·√(1-x_c²)) 的对比：")
    for fc in NEAR_CUTOFFS:
        _, sigma = TP.omega_sigma(fc, FS)
        x_c = 2.0 * fc / FS
        pin = x_c / (np.tan(np.pi * x_c / 2.0) * np.sqrt(1.0 - x_c ** 2))
        print(f"  fc={fc:.4f}: A={sigma / np.pi:.9f}  pin={pin:.9f}  "
              f"diff={sigma / np.pi - pin:+.1e}")
    print("贴 Nyquist：x = 0.65 / 0.75 / 0.85 / 0.95 / 1.00 处的 y（Hz，fs = 2）")
    for fc in NEAR_CUTOFFS:
        omega, sigma = TP.omega_sigma(fc, FS)
        y = fc * rho(np.pi * np.array([0.65, 0.75, 0.85, 0.95, 1.0]), omega, sigma)
        print(f"  fc={fc:.4f} (u={FS / fc:7.3f}, A={sigma / np.pi:.4f}): "
              + " ".join(f"{v:8.5f}" for v in y))


def main() -> int:
    # 12 条曲线按顺序均匀取色（截止频率在 4 k–24 k 段很密，用 fc 直接上色会挤在一起）
    n = len(CASES)
    colors = [plt.get_cmap("viridis")(v) for v in np.linspace(0.08, 0.95, n)]
    listed = matplotlib.colors.ListedColormap(colors)
    norm = matplotlib.colors.BoundaryNorm(np.arange(n + 1), n)
    scalar = plt.cm.ScalarMappable(cmap=listed, norm=norm)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))

    ax = axes[0]
    ax.plot([0.0, 1.0], [0.0, 1.0], color="0.4", lw=1.0, ls="--", label="恒等映射 y = x")
    worst = 0.0
    x = np.linspace(1e-9, FS / 2.0, 4001)
    for ref, fc in CASES:
        omega, sigma = TP.omega_sigma(fc, FS)
        color = colors[CASES.index((ref, fc))]
        ax.plot(x, fc * rho(np.pi * x / FS * 2.0, omega, sigma), color=color, lw=1.8)
        xn, yn = numeric_mapping(fc, FS)
        worst = max(worst, float(np.max(np.abs(yn - fc * rho(np.pi * xn / FS * 2.0,
                                                            omega, sigma)))))
    xb, yb = bilinear(CASES[0][1], FS)
    keep = yb <= 1.15
    ax.plot(xb[keep], yb[keep], color="0.55", lw=1.2, ls=":",
            label="双线性 σ=0（超出 1.15 已截断）")
    ax.text(0.03, 0.66, "σ=0 在 Nyquist 处发散（y → ∞）", color="0.45", fontsize=8.5)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.15)
    ax.set_xlabel("数字频率 [Hz]（fs = 2，Nyquist = 1）")
    ax.set_ylabel("等效模拟频率 [Hz]")
    ax.set_title("解析映射 y = fc·ρ(θ)（与数值反解的偏差 < 1e-7）")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)
    cb = fig.colorbar(scalar, ax=ax, ticks=np.arange(n) + 0.5, pad=0.015)
    cb.ax.set_yticklabels([f"{ref:g}" for ref, _ in CASES], fontsize=7)
    cb.set_label("截止频率 [Hz]")

    ax = axes[1]
    xs = np.linspace(1e-4, FS / 2.0 * 0.9999, 4001)
    for index, (ref, fc) in enumerate(CASES):
        omega, sigma = TP.omega_sigma(fc, FS)
        theta = np.pi * xs / FS * 2.0
        slope = np.gradient(fc * rho(theta, omega, sigma), xs)
        ax.plot(xs, np.minimum(slope, 8.0), color=colors[index], lw=1.6)
    ax.axhline(1.0, color="0.4", lw=1.0, ls="--", label="恒等映射 dy/dx = 1")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 8.0)
    ax.set_xlabel("数字频率 [Hz]（fs = 2，Nyquist = 1）")
    ax.set_ylabel("局部斜率 dy/dx")
    ax.set_title("局部斜率：<1 模拟轴收缩，>1 模拟轴拉伸")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)
    cb = fig.colorbar(scalar, ax=ax, ticks=np.arange(n) + 0.5, pad=0.015)
    cb.ax.set_yticklabels([f"{ref:g}" for ref, _ in CASES], fontsize=7)
    cb.set_label("截止频率 [Hz]")

    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    path = OUT / "two_point_mapping.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"图：{path}")
    print(f"解析式与数值反解的最大偏差：{worst:.2e}")

    print("x = 0.25 / 0.50 / 0.75 / 0.90 / 0.95 处的 y（Hz，fs = 2）：")
    for ref, fc in CASES:
        omega, sigma = TP.omega_sigma(fc, FS)
        y = fc * rho(np.pi * np.array([0.25, 0.5, 0.75, 0.9, 0.95]), omega, sigma)
        print(f"  {ref:8.1f} Hz (u={FS / fc:8.2f}): " + " ".join(f"{v:8.5f}" for v in y))
    near_nyquist()
    print("低 fc 极限（omega≈u, sigma≈2）时 y 与 fc 无关：")
    print("  y = (2/pi) t / sqrt(1 + (2t/pi)^2),  t = tan(theta/2)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
