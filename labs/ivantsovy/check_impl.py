"""正确性自检：论文公式、四种经典类型、预畸变、稳定性与零点位置。

用法: python qwqdsp/labs/ivantsovy/check_impl.py
退出码 0 表示全部通过。
"""
from __future__ import annotations

import numpy as np
from scipy import signal

import ivantsovy as I

FS = 48000.0
F_REF = 5000.0
ZETA = 1.0 / np.sqrt(2.0)
SIGMAS = (0.0, 1.0, np.pi * np.sqrt(2.0 / 3.0), 2.0)
ZOOM_SIGMA = np.pi * np.sqrt(2.0 / 3.0)


def require(cond: bool, message: str) -> None:
    if not cond:
        raise AssertionError(message)


def poly_from_factors(factors, gain: float = 1.0) -> np.ndarray:
    out = np.array([gain])
    for fac in factors:
        out = np.convolve(out, fac.gain * fac.poly / fac.norm)
    return out


def design_tf(zeros, poles, fs: float, f_ref: float, sigma: float,
              prewarp: bool = False, gain: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    design = I.design(zeros, poles, fs, f_ref, sigma, prewarp=prewarp, gain=gain)
    b = poly_from_factors(design.num, design.gain)
    a = poly_from_factors(design.den)
    return b / a[0], a / a[0]


def check_identity() -> None:
    for x in (0.0, 0.3, 1.0, 9.6, 40.0):
        for y in (0.0, 0.5, 1.0, 2.0):
            for sigma in (0.5, 2.0, ZOOM_SIGMA, 3.0):
                p1, p2 = I.phi12(x, y, sigma)
                x2, y2 = x * x, y * y
                nu = np.sqrt(x2 * x2 + 2 * x2 * sigma**2 * (2 * y2 - 1) + sigma**4)
                kappa = x2 * (2 * y2 - 1) + sigma**2
                d = np.pi**2 + np.pi * np.sqrt(2) * np.sqrt(max(nu + kappa, 0.0)) + nu
                require(abs(1 + p1 + p2 - 4 * np.pi**2 / d) < 1e-12,
                        f"identity failed at x={x}, y={y}, sigma={sigma}")
    print("PASS  1 + phi1 + phi2 = 4*pi^2/d  (论文结论)")


def check_first_order() -> None:
    for sigma in SIGMAS:
        for fs, fc in ((48000.0, 1000.0), (44100.0, 10000.0)):
            z, p, _ = signal.buttap(1)
            b, a = design_tf(z, p, fs, fc, sigma)
            omega = fs / fc
            alpha = I.phi(0.0, sigma)
            beta = I.phi(omega, sigma)
            g = 1.0 / (1.0 + alpha)
            b_ref = g * (1.0 + beta) * np.array([1.0, alpha])
            a_ref = np.array([1.0, beta])
            err = max(np.max(np.abs(b - b_ref)), np.max(np.abs(a - a_ref)))
            require(err < 1e-12,
                    f"first-order LP mismatch sigma={sigma}, fs={fs}, fc={fc}: {err:.3e}")

            b, a = design_tf(np.array([0.0]), np.array([-1.0]), fs, fc, sigma)
            g = omega / (2.0 * np.pi)
            b_ref = g * (1.0 + beta) * np.array([1.0, -1.0])
            err = max(np.max(np.abs(b - b_ref)), np.max(np.abs(a - a_ref)))
            require(err < 1e-12,
                    f"first-order HP mismatch sigma={sigma}, fs={fs}, fc={fc}: {err:.3e}")
    print("PASS  一阶 LP / HP = 论文 (2.8) / Table 2.1")


def second_order_reference(kind: str, sigma: float,
                           omega: float) -> tuple[np.ndarray, np.ndarray]:
    phi0 = I.phi(0.0, sigma)
    beta1, beta2 = I.phi12(omega, ZETA, sigma)
    if kind == "lp":
        alpha1, alpha2 = 2.0 * phi0, phi0 * phi0
        g = 1.0 / (1.0 + alpha1 + alpha2)
    elif kind == "hp":
        alpha1, alpha2 = -2.0, 1.0
        g = omega * omega / (4.0 * np.pi * np.pi)
    elif kind == "bp":
        alpha1, alpha2 = phi0 - 1.0, -phi0
        g = (omega / (2.0 * np.pi)) * 2.0 * ZETA / (2.0 + alpha1)
    elif kind == "bs":
        alpha1, alpha2 = I.phi12(omega, 0.0, sigma)
        g = 1.0 / (1.0 + alpha1 + alpha2)
    else:
        raise ValueError(kind)
    b_ref = g * (1.0 + beta1 + beta2) * np.array([1.0, alpha1, alpha2])
    a_ref = np.array([1.0, beta1, beta2])
    return b_ref, a_ref


def check_second_order() -> None:
    pole_pair = np.array([-ZETA + 1j * np.sqrt(1.0 - ZETA**2),
                          -ZETA - 1j * np.sqrt(1.0 - ZETA**2)])
    zeros = {
        "lp": np.array([]),
        "hp": np.array([0.0, 0.0]),
        "bp": np.array([0.0]),
        "bs": np.array([1j, -1j]),
    }
    gains = {"lp": 1.0, "hp": 1.0, "bp": 2.0 * ZETA, "bs": 1.0}

    for sigma in (ZOOM_SIGMA, 2.0):
        for fs, fc in ((48000.0, 1000.0), (44100.0, 10000.0)):
            omega = fs / fc
            for kind in ("lp", "hp", "bp", "bs"):
                b, a = design_tf(zeros[kind], pole_pair, fs, fc, sigma,
                                 gain=gains[kind])
                b_ref, a_ref = second_order_reference(kind, sigma, omega)
                err = max(np.max(np.abs(b - b_ref)), np.max(np.abs(a - a_ref)))
                require(err < 1e-12,
                        f"2nd-order {kind} mismatch sigma={sigma}, fs={fs}, fc={fc}: {err:.3e}")
    print("PASS  二阶 LP / HP / BP / BS = 论文 (3.8) / Table 3.1")


def check_four_pole_example() -> None:
    for k_gain in (1.0, 3.0):
        a1 = 1.0 / np.sqrt(k_gain**2 + np.sqrt(2.0) * k_gain + 1.0)
        a2 = 1.0 / np.sqrt(k_gain**2 - np.sqrt(2.0) * k_gain + 1.0)
        b1 = 2.0 + np.sqrt(2.0) * k_gain
        b2 = 2.0 - np.sqrt(2.0) * k_gain
        roots = np.concatenate([np.roots([an**2, an**2 * bn, 1.0])
                                for an, bn in ((a1, b1), (a2, b2))])
        for sigma in (0.0, 2.0, ZOOM_SIGMA):
            f_ref = 1.0 / (2.0 * np.pi)  # 论文例子里的 s 已按 1 rad/s 归一化
            omega = FS / f_ref
            design = I.design(np.array([]), roots, FS, f_ref, sigma)
            b = poly_from_factors(design.num, design.gain)
            a = poly_from_factors(design.den)

            c0 = I.phi(0.0, sigma)
            b_ref = np.array([1.0])
            for _ in range(4):
                b_ref = np.convolve(b_ref, np.array([1.0, c0]) / (1.0 + c0))
            a_ref = np.array([1.0])
            for an, bn in ((a1, b1), (a2, b2)):
                p1, p2 = I.phi12(omega * an, 0.5 * an * bn, sigma)
                if bn >= 0.0:
                    c1, c2 = p1, p2
                else:
                    c1, c2 = p1 / p2, 1.0 / p2
                a_ref = np.convolve(a_ref, np.array([1.0, c1, c2]) / (1.0 + c1 + c2))

            scale = max(np.max(np.abs(b_ref)), np.max(np.abs(a_ref)), 1.0)
            err = max(np.max(np.abs(b - b_ref)), np.max(np.abs(a - a_ref))) / scale
            require(err < 1e-12,
                    f"example 4.10 mismatch K={k_gain}, sigma={sigma}: {err:.3e}")
    print("PASS  四极点例子 = 论文 (4.10)（含 b<0 / 右半平面极点分支）")


def check_bilinear_limit() -> None:
    for order in (1, 2, 5, 8):
        fs, fc = 48000.0, 3000.0
        z, p, _ = signal.buttap(order)
        z = I.complete_conjugates(z)
        p = I.complete_conjugates(p)
        b_s = np.poly(z) if len(z) else np.array([1.0])
        a_s = np.poly(p)
        b_s = b_s / b_s[-1]
        a_s = a_s / a_s[-1]
        wc = 2.0 * np.pi * fc
        b_scaled = np.array([c / wc**i for i, c in enumerate(b_s[::-1])])[::-1]
        a_scaled = np.array([c / wc**i for i, c in enumerate(a_s[::-1])])[::-1]
        b_ref, a_ref = signal.bilinear(b_scaled, a_scaled, fs)

        b, a = design_tf(z, p, fs, fc, 0.0)
        f = np.linspace(0.0, fs / 2.0, 2000)
        _, h_here = signal.freqz(b, a, worN=2.0 * np.pi * f / fs)
        _, h_ref = signal.freqz(b_ref, a_ref, worN=2.0 * np.pi * f / fs)
        err = np.max(np.abs(np.abs(h_here) - np.abs(h_ref)))
        require(err < 1e-8,
                f"sigma=0 != bilinear for order {order}: {err:.3e}")
    print("PASS  sigma = 0 退化为标准双线性变换（N = 1,2,5,8）")


def check_prewarp() -> None:
    for sigma in (2.0, ZOOM_SIGMA):
        omega = I.prewarped_omega(F_REF, FS, sigma)
        theta = np.arccos((omega * omega - np.pi**2 - sigma**2) /
                          (omega * omega + np.pi**2 - sigma**2))
        f_mapped = FS / (2.0 * np.pi) * theta
        require(abs(f_mapped - F_REF) < 1e-9,
                f"F(omega) != f_ref for sigma={sigma}: {f_mapped}")

        f_match = F_REF * omega / sigma
        for order in (1, 2, 5, 8, 10):
            z, p, _ = signal.buttap(order)
            z = I.complete_conjugates(z)
            p = I.complete_conjugates(p)
            design = I.design(z, p, FS, F_REF, sigma, prewarp=True)
            hd_ref = I.digital_magnitude(design, np.array([F_REF]), FS)[0]
            require(abs(hd_ref - 1.0 / np.sqrt(2.0)) < 1e-12,
                    f"prewarp cutoff mismatch: order={order}, sigma={sigma}")

            hd_nyq = I.digital_magnitude(design, np.array([FS / 2.0]), FS)[0]
            ha_match = I.analog_magnitude(z, p, np.array([f_match]), F_REF)[0]
            require(abs(hd_nyq - ha_match) <= 1e-9 * ha_match,
                    f"prewarp Nyquist mismatch: order={order}, sigma={sigma}")

    z, p, _ = signal.cheb2ap(8, 80.0)
    z = I.complete_conjugates(z)
    p = I.complete_conjugates(p)
    design = I.design(z, p, FS, F_REF, 2.0, prewarp=True)
    hd_ref = I.digital_magnitude(design, np.array([F_REF]), FS)[0]
    ha_ref = I.analog_magnitude(z, p, np.array([F_REF]), F_REF)[0]
    require(abs(hd_ref - ha_ref) <= 1e-12 * ha_ref,
            "Chebyshev II prewarp reference mismatch")
    print("PASS  预畸变：F(omega)=f_ref；Butterworth 任意阶 Hd(f_ref)=1/sqrt(2)；"
          "Nyquist=A(f_ref*omega/sigma)；Chebyshev II 参考点精确对齐")


def check_properties() -> None:
    cases: list[tuple[str, np.ndarray, np.ndarray]] = []
    for order in range(1, 11):
        z, p, _ = signal.buttap(order)
        cases.append((f"Butterworth-{order}", z, p))
        for rs in (40.0, 60.0, 80.0, 100.0):
            z, p, _ = signal.cheb2ap(order, rs)
            cases.append((f"ChebyshevII-{order}-{rs:g}dB", z, p))

    for name, z, p in cases:
        z = I.complete_conjugates(z)
        p = I.complete_conjugates(p)
        for fc in (2000.0, 5000.0, 10000.0):
            design = I.design(z, p, FS, fc, ZOOM_SIGMA)
            dc = I.digital_magnitude(design, np.array([0.0]), FS)[0]
            require(abs(dc - 1.0) < 1e-12, f"DC gain != 1: {name}, fc={fc}")

            nyq = I.digital_magnitude(design, np.array([FS / 2.0]), FS)[0]
            ref = I.analog_magnitude(z, p, np.array([FS / ZOOM_SIGMA]), fc)[0]
            require(abs(nyq - ref) <= 1e-9 * max(ref, 1e-300),
                    f"Nyquist mismatch: {name}, fc={fc}: {nyq} vs {ref}")

            _, poles = I.digital_zeros_poles(design)
            require(np.max(np.abs(poles)) < 1.0,
                    f"unstable digital pole: {name}, fc={fc}")

            zeros, _ = I.digital_zeros_poles(design)
            if len(zeros):
                require(np.max(np.abs(zeros)) <= 1.0 + 1e-9,
                        f"digital zero outside unit circle: {name}, fc={fc}")
    print(f"PASS  {len(cases) * 3} 组低通原型：DC = 1、Nyquist = A(fs/sigma)、稳定、零点在单位圆内")


def check_four_types_response() -> None:
    """高阶 LP/HP/BP/BS：与 scipy 模拟原型对照。

    Butterworth N=4；Chebyshev II N=8、rs=80 dB。Chebyshev II 的 zpk 增益
    必须经 normalized_gain 换算，否则通带 / 阻带电平会整体偏移。
    """
    wc = 2.0 * np.pi * F_REF
    specs = (
        ("lowpass", wc),
        ("highpass", wc),
        ("bandpass", [2.0 * np.pi * 4000.0, 2.0 * np.pi * 6000.0]),
        ("bandstop", [2.0 * np.pi * 4000.0, 2.0 * np.pi * 6000.0]),
    )
    families = (
        ("Butterworth-4", lambda wn, bt: signal.butter(
            4, wn, btype=bt, analog=True, output="zpk")),
        ("ChebyshevII-8-80dB", lambda wn, bt: signal.cheby2(
            8, 80.0, wn, btype=bt, analog=True, output="zpk")),
        ("ChebyshevII-7-80dB", lambda wn, bt: signal.cheby2(
            7, 80.0, wn, btype=bt, analog=True, output="zpk")),
    )
    f = np.linspace(20.0, FS / 2.0, 20000)
    f_check = np.array([20.0, 1000.0, F_REF, FS / 3.0])
    for family, prototype in families:
        print(f"       {family}:")
        for btype, wn in specs:
            z0, p0, k0 = prototype(wn, btype)
            z0 = np.asarray(z0, dtype=complex)
            p0 = np.asarray(p0, dtype=complex)
            z = z0 / wc
            p = p0 / wc
            k_s = k0 * wc ** (len(z0) - len(p0))
            gain = I.normalized_gain(z, p, k_s)

            _, h_ref = signal.freqs_zpk(z0, p0, k0, 2.0 * np.pi * f_check)
            ha_check = I.analog_magnitude(z, p, f_check, F_REF, gain=gain)
            require(np.max(np.abs(ha_check - np.abs(h_ref))) <=
                    1e-9 * max(np.max(np.abs(h_ref)), 1e-300),
                    f"{family} {btype}: normalized gain mismatch")

            design = I.design(z, p, FS, F_REF, ZOOM_SIGMA, gain=gain)
            ha = I.analog_magnitude(z, p, f, F_REF, gain=gain)
            hd = I.digital_magnitude(design, f, FS)
            err = np.max(np.abs(hd - ha))
            _, poles = I.digital_zeros_poles(design)
            require(np.max(np.abs(poles)) < 1.0, f"unstable {family} {btype}")
            require(err < 0.1, f"{family} {btype} response deviation too large: {err:.3e}")
            print(f"         {btype:8}: gain={gain:.2e}, max|Hd-Ha| = {err:.3e}")
    print("PASS  高阶 LP / HP / BP / BS：Butterworth N=4、Chebyshev II N=8/7 与 scipy 一致")


def check_imaginary_axis_zeros() -> None:
    z, p, _ = signal.cheb2ap(8, 80.0)
    z = I.complete_conjugates(z)
    p = I.complete_conjugates(p)
    design = I.design(z, p, FS, 5000.0, ZOOM_SIGMA)
    zeros, _ = I.digital_zeros_poles(design)
    on_circle = zeros[np.abs(np.abs(zeros) - 1.0) <= 1e-12]
    real_axis = zeros[np.abs(zeros.imag) <= 1e-9]
    require(len(on_circle) == 6,
            f"expected 6 zeros on the unit circle, got {len(on_circle)}")
    require(len(real_axis) == 2,
            f"expected 2 real zeros, got {len(real_axis)}")
    require(np.max(np.abs(zeros)) <= 1.0 + 1e-12,
            "Chebyshev II zero outside the unit circle")

    analog_hz = np.sort(np.abs(z) * 5000.0)
    mapped_hz = np.sort(np.abs(np.angle(on_circle)) / (2.0 * np.pi) * FS)
    ref_hz = analog_hz[analog_hz < FS / ZOOM_SIGMA]
    err = float(np.max(np.abs(mapped_hz - ref_hz) / ref_hz))
    require(err < 1e-2,
            f"unit-circle null frequencies shifted by {err:.2%}")
    print("PASS  Chebyshev II 虚轴零点：6 个落在单位圆上，2 个（高于 fs/sigma）落到实轴；"
          f"零频最大相对误差 {err:.2%}")


def check_origin_pole_rejected() -> None:
    try:
        I.design(np.array([]), np.array([0.0]), FS, 1000.0, ZOOM_SIGMA)
    except ValueError:
        print("PASS  s = 0 处的极点（积分器）被明确拒绝")
        return
    raise AssertionError("s = 0 处的极点本应被拒绝")


def main() -> int:
    check_identity()
    check_first_order()
    check_second_order()
    check_four_pole_example()
    check_bilinear_limit()
    check_prewarp()
    check_properties()
    check_four_types_response()
    check_imaginary_axis_zeros()
    check_origin_pole_rejected()
    print("\n全部检查通过。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
