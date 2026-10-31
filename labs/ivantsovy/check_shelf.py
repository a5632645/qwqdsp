"""High-shelf 原型与 §4 离散化的正确性自检。

用法: python qwqdsp/labs/ivantsovy/check_shelf.py
退出码 0 表示全部通过。
"""
from __future__ import annotations

import numpy as np

import ivantsovy as I
import shelf as S

FS = 48000.0
SIGMA = np.pi * np.sqrt(2.0 / 3.0)


def require(cond: bool, message: str) -> None:
    if not cond:
        raise AssertionError(message)


def to_db(magnitude: float) -> float:
    return float(20.0 * np.log10(max(abs(magnitude), 1e-300)))


def levels(zeros, poles, f_ref: float) -> tuple[float, float, float]:
    low = to_db(I.analog_magnitude(zeros, poles, np.array([1e-4 * f_ref]), f_ref)[0])
    mid = to_db(I.analog_magnitude(zeros, poles, np.array([f_ref]), f_ref)[0])
    high = to_db(I.analog_magnitude(zeros, poles, np.array([1e6 * f_ref]), f_ref)[0])
    return low, mid, high


def check_prototypes() -> None:
    f_ref = 1000.0
    for order in (4, 8):
        for boost in (12.0, -12.0, 24.0, -24.0):
            zeros, poles = S.butterworth_shelf(order // 2, boost)
            low, mid, high = levels(zeros, poles, f_ref)
            require(abs(low) < 1e-9, f"Butter shelf DC != 0 dB: {low}")
            require(abs(mid - boost / 2.0) < 1e-6,
                    f"Butter shelf midpoint mismatch: {mid} vs {boost / 2.0}")
            require(abs(high - boost) < 1e-6,
                    f"Butter shelf top mismatch: {high} vs {boost}")

    for order in (4, 8):
        for boost in (12.0, -12.0, 24.0, -24.0):
            for ripple in (0.5, 1.0, 3.0):
                zeros, poles = S.chebyshev2_shelf(order // 2, boost, ripple)
                low, mid, high = levels(zeros, poles, f_ref)
                sign = 1.0 if boost > 0.0 else -1.0
                expected_high = boost - sign * ripple
                require(abs(low) < 1e-9, f"Cheb2 shelf DC != 0 dB: {low}")
                require(abs(mid - boost / 2.0) < 0.01,
                        f"Cheb2 midpoint mismatch: {mid} vs {boost / 2.0}")
                require(abs(high - expected_high) < 0.02,
                        f"Cheb2 high end mismatch: {high} vs {expected_high}")

                f = np.geomspace(0.01 * f_ref, 1e3 * f_ref, 200000)
                h = I.analog_magnitude(zeros, poles, f, f_ref)
                db = 20.0 * np.log10(np.maximum(h, 1e-300))
                extreme = db.max() if boost > 0.0 else db.min()
                require(abs(extreme - boost) < 0.05,
                        f"Cheb2 shelf extreme {extreme} != {boost} dB")
    print("PASS  shelf 三电平：DC 0 dB、中点 boost/2、Butterworth 高频 boost、"
          "Chebyshev II 高频 boost-ripple 且棚带极值到 boost")


def check_digital() -> None:
    f = np.linspace(20.0, FS / 2.0, 100001)
    rows: list[tuple[str, int, float, float, float, float]] = []

    def run(family: str, order: int, boost: float, ripple: float, f_ref: float,
            limit: float) -> None:
        if family == "butter":
            zeros, poles = S.butterworth_shelf(order // 2, boost)
        else:
            zeros, poles = S.chebyshev2_shelf(order // 2, boost, ripple)
        design = I.design(zeros, poles, FS, f_ref, SIGMA, gain=1.0, prewarp=True)
        analog = I.analog_magnitude(zeros, poles, f, f_ref)
        digital = I.digital_magnitude(design, f, FS)
        err = 20.0 * np.log10(np.maximum(digital, 1e-300)) - 20.0 * np.log10(
            np.maximum(analog, 1e-300))
        max_db = float(np.max(np.abs(err)))
        _, poles_d = I.digital_zeros_poles(design)
        require(np.max(np.abs(poles_d)) < 1.0,
                f"unstable shelf: {family} N={order}, boost={boost}, f_ref={f_ref}")
        require(max_db < limit,
                f"shelf deviation too large: {family} N={order}, boost={boost}, "
                f"ripple={ripple}, f_ref={f_ref}: {max_db:.3f} dB")
        rows.append((family, order, boost, ripple, f_ref, max_db))

    for order in (4, 8):
        for boost in (12.0, -12.0, 24.0, -24.0):
            for f_ref in (1000.0, 5000.0):
                run("butter", order, boost, 0.0, f_ref, 0.05)

    for order in (4, 8):
        for boost in (12.0, -12.0):
            for ripple in (0.5, 1.0, 3.0):
                for f_ref in (1000.0, 5000.0):
                    run("cheb2", order, boost, ripple, f_ref, 1.5)

    run("cheb2", 16, 24.0, 1.0, 5000.0, 1.0)
    run("cheb2", 16, 24.0, 3.0, 5000.0, 3.0)

    print(f"{'family':<8}{'N':>4}{'boost':>8}{'ripple':>8}{'f_ref':>8}{'max dB':>10}")
    for family, order, boost, ripple, f_ref, max_db in rows:
        print(f"{family:<8}{order:>4}{boost:>8.0f}{ripple:>8.1f}{f_ref:>8.0f}{max_db:>10.3f}")
    worst = max(row[-1] for row in rows)
    print(f"PASS  {len(rows)} 组 shelf：稳定，且与模拟原型的最大 dB 偏差 {worst:.3f} dB")


def main() -> int:
    check_prototypes()
    check_digital()
    print("\n全部 shelf 检查通过。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
