"""Run constrained rational Remez for the analog-matching filter families."""

import argparse
import csv
from dataclasses import dataclass
import importlib
from itertools import product
import json
from math import comb, fsum
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial import polynomial as poly
from scipy import optimize, signal

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "constraint_fit"), str(ROOT / "fit")]
from fit_curve import MagnitudeDb, _UnitDiskRoots
from remez import Exchange, Extrema, Level, SelectReference


@dataclass(frozen=True)
class Family:
    numerator_degree: int = 2
    denominator_degree: int = 2
    dc_zeros: int = 0
    notch: bool = False
    anchors: tuple = ()


FAMILIES = {
    "lowpass_1st": Family(),
    "highpass_1st": Family(1, dc_zeros=1),
    "lowpass_2nd": Family(anchors=("dc",)),
    "highpass_2nd": Family(0, dc_zeros=2),
    "bandpass_2nd": Family(1, dc_zeros=1),
    "peaking_2nd": Family(anchors=("dc", "center")),
    "highshelf_1st": Family(anchors=("dc", "nyquist")),
    "highshelf_2nd": Family(anchors=("dc",)),
    "notch_2nd": Family(0, notch=True, anchors=("dc",)),
    "notch_3rd": Family(1, 3, notch=True, anchors=("dc",)),
}


def PowerVariable(frequency):
    return np.sin(0.5 * np.pi * frequency) ** 2


def _ReducedTarget(name, module, frequency, fc, parameters):
    """Analog power divided by the mandatory digital zero factors."""
    wc = 2 * np.pi * fc / module.FS
    w = np.pi * frequency
    ratio = w / wc
    if name == "lowpass_1st":
        return 1 / (1 + ratio**2)
    if name in ("highpass_1st", "highpass_2nd", "bandpass_2nd"):
        r2_over_x = (2 / (wc * np.sinc(frequency / 2))) ** 2
        if name == "highpass_1st":
            return r2_over_x / (1 + ratio**2)
        q = parameters[0]
        denominator = (1 - ratio**2)**2 + (ratio / q)**2
        numerator = r2_over_x**2 if name == "highpass_2nd" else r2_over_x / q**2
        return numerator / denominator
    if FAMILIES[name].notch:
        q = parameters[0]
        xc = PowerVariable(2 * fc / module.FS)
        # Cancel the common zero using the sine difference identity, even at fc.
        ratio_over_notch = (2 * xc * (w + wc) /
                            (wc**2 * np.sin((w + wc) / 2) * np.sinc((w - wc) / (2 * np.pi))))
        return ratio_over_notch**2 / ((1 - ratio**2)**2 + (ratio / q)**2)
    return module.analog_power(w, wc, *parameters)


def _PowerCoefficients(model, family):
    if family.denominator_degree == 3:
        return np.array([1.0, model.c]), np.r_[1.0, model.d]

    def Convert(coefficients):
        start, end, middle = coefficients / model.alpha[0]
        return np.array([start, end - start + 4 * middle, -4 * middle])

    numerator = np.array([1.0]) if family.notch else Convert(model.beta)[family.dc_zeros:]
    return numerator, Convert(model.alpha)


def _ReversePowerVariable(coefficients):
    # Change P(x) to P(1-y), with accurate cancellation at the new origin.
    return np.array([(-1)**i * fsum(comb(j, i) * coefficients[j]
                                  for j in range(i, len(coefficients)))
                     for i in range(len(coefficients))])


def _Roots(coefficients, reverse):
    roots = poly.polyroots(coefficients)
    return _UnitDiskRoots(2 * roots - 1 if reverse else 1 - 2 * roots)


def _AsSos(model, family, center, reverse):
    numerator, denominator = model
    zeros, poles = _Roots(numerator, reverse), _Roots(denominator, reverse)
    if zeros is None or poles is None or numerator[0] <= 0 or denominator[0] <= 0:
        return None
    fixed = np.repeat(1.0, family.dc_zeros)
    if family.notch:
        notch_roots = np.exp(1j * np.pi * center * np.array([1.0, -1.0]))
        fixed = np.r_[fixed, notch_roots]
    sos = signal.zpk2sos(np.r_[zeros, fixed], poles, 1.0)
    # Normalize away from all exact zeros.
    ref = np.array([0.75 if abs(center - 0.25) < 0.1 else 0.25])
    x = PowerVariable(ref)
    y = PowerVariable(1 - ref) if reverse else x
    power = poly.polyval(y, numerator) / poly.polyval(y, denominator) * x**family.dc_zeros
    if family.notch:
        power *= (1 - x / PowerVariable(center))**2
    gain_db = 10 * np.log10(power[0]) - MagnitudeDb(sos, ref)[0]
    if not np.isfinite(gain_db):
        return None
    sos[0, :3] *= 10 ** (gain_db / 20)
    return sos


def _GridSeed(family, variable, target, grid, anchors, ceiling, valid):
    """Find a feasible starting model for another Remez exchange."""
    split = family.numerator_degree + 1
    numerator = poly.polyvander(variable(grid), family.numerator_degree)
    denominator = poly.polyvander(variable(grid), family.denominator_degree)[:, 1:]
    power = target(grid)
    fixed_power = target(anchors)
    equality = np.c_[poly.polyvander(variable(anchors), family.numerator_degree),
                     -fixed_power[:, None] *
                     poly.polyvander(variable(anchors), family.denominator_degree)[:, 1:]]
    active = np.unique(np.r_[np.arange(0, grid.size, 32), grid.size - 1])
    bands = np.array_split(np.arange(grid.size), 32)
    floor, best = 0.0, None
    for _ in range(28):
        if ceiling - floor < 2e-8:
            break
        limit = (floor + ceiling) / 2
        ratio = 10 ** (limit / 10)
        n, d, p = numerator[active], denominator[active], power[active]
        low, high = p / ratio, p * ratio
        rows = np.vstack((np.c_[n, -high[:, None] * d],
                          np.c_[-n, low[:, None] * d],
                          np.c_[np.zeros_like(n), -d]))
        rhs = np.r_[high, -low, np.full(active.size, 1 - 1e-9)]
        scale = np.maximum(1.0, np.abs(rhs))
        fit = optimize.linprog(np.zeros(rows.shape[1]), A_ub=rows / scale[:, None],
                               b_ub=rhs / scale, A_eq=equality, b_eq=fixed_power,
                               bounds=(None, None), method="highs",
                               options={"primal_feasibility_tolerance": 1e-9})
        if fit.status == 2:
            floor = limit
            continue
        if not fit.success:
            break
        model = fit.x[:split], np.r_[1.0, fit.x[split:]]
        n_eval, d_eval = numerator @ model[0], 1 + denominator @ fit.x[split:]
        with np.errstate(invalid="ignore", divide="ignore"):
            error = np.abs(10 * np.log10(n_eval / (d_eval * power)))
        error = np.where((n_eval > 0) & (d_eval > 0), error, np.inf)
        peak = np.max(error)
        if peak > limit + 1e-6:
            augmented = np.union1d(active, [band[np.argmax(error[band])] for band in bands])
            if augmented.size == active.size:
                break
            active = augmented
            continue
        if peak >= ceiling or not valid(model):
            break
        best, ceiling = model, peak
    return best


def Design(name, fc, parameters, *, grid_size=8193):
    """Return the BRASIL seed and Remez-refined constrained SOS with error metrics."""
    family = FAMILIES[name]
    module = importlib.import_module(f"{name}_interp")
    center = 2 * fc / module.FS
    reverse = center > 0.5
    locations = {"dc": 0.0, "center": center, "nyquist": 1.0}
    anchors = np.array([locations[key] for key in family.anchors])
    split = family.numerator_degree + 1
    count = split + family.denominator_degree + 1 - len(anchors)

    def Target(frequency):
        return _ReducedTarget(name, module, frequency, fc, parameters)

    def Variable(frequency):
        return PowerVariable(1 - frequency if reverse else frequency)

    grid = np.unique(np.r_[np.linspace(0, 1, grid_size),
                           np.geomspace(1e-8, 1.0, grid_size // 4),
                           1 - np.geomspace(1e-8, 1.0, grid_size // 4), center])
    reference = Target(grid)

    def Error(model, frequency):
        desired = reference if frequency is grid else Target(frequency)
        x = Variable(frequency)
        with np.errstate(divide="ignore", invalid="ignore"):
            return 10 * np.log10(poly.polyval(x, model[0]) / (poly.polyval(x, model[1]) * desired))

    def Orientation(frequency):
        # Every difference of two admissible responses has the anchored center
        # as a known zero. Remove its sign flip before selecting the extrema.
        return np.where(frequency < center, -1.0, 1.0) if "center" in family.anchors else 1.0

    def ExchangeError(model, frequency):
        return Orientation(frequency) * Error(model, frequency)

    def Candidates(frequency, signs, peak):
        all_frequency = np.r_[frequency, anchors]
        desired = Target(all_frequency)
        x = Variable(all_frequency)
        numerator_basis = poly.polyvander(x, family.numerator_degree)
        denominator_basis = poly.polyvander(x, family.denominator_degree)
        target_signs = np.r_[signs * Orientation(frequency), np.zeros(len(anchors))]

        def Matrix(ripple):
            target = desired * 10.0 ** (ripple * target_signs / 10.0)
            scale = np.maximum(1.0, target)
            return np.c_[numerator_basis / scale[:, None],
                         -denominator_basis * (target / scale)[:, None]]

        for ripple, vector in Level(Matrix, peak, samples=33):
            if vector[split] == 0.0:
                continue
            model = tuple(np.split(vector / vector[split], [split]))
            residual = Error(model, all_frequency)
            if (np.all(np.isfinite(residual))
                    and np.max(np.abs(residual - ripple * target_signs)) <= 1e-7
                    and _AsSos(model, family, center, reverse) is not None):
                yield model

    def Seed(match):
        model = _PowerCoefficients(match, family)
        return tuple(_ReversePowerVariable(c) for c in model) if reverse else model

    match, _, _ = module.design(fc, *parameters)
    if match is None:
        return None
    initial = Seed(match)
    best = Exchange(initial, ExchangeError, Candidates, count, grid, tol=1e-9)
    if hasattr(module, "DB_LIST"):
        _, peaks = Extrema(lambda f: ExchangeError(best, f), grid)
        peak = np.max(np.abs(peaks))
        refs = SelectReference(peaks, count)
        if peak > 1e-5 and (refs is None or np.min(np.abs(peaks[refs])) < 0.99 * peak):
            opposite, _, _ = module.design(fc, *parameters[:-1], -parameters[-1])
            if opposite is not None:
                reciprocal = Seed(opposite)[::-1]
                alternative = Exchange(reciprocal, ExchangeError, Candidates, count, grid, tol=1e-9)
                _, alternate_peaks = Extrema(lambda f: ExchangeError(alternative, f), grid)
                if np.max(np.abs(alternate_peaks)) < peak:
                    best = alternative
    _, peaks = Extrema(lambda f: ExchangeError(best, f), grid)
    peak = np.max(np.abs(peaks))
    refs = SelectReference(peaks, count)
    if peak > 1e-5 and (refs is None or np.min(np.abs(peaks[refs])) < 0.99 * peak):
        warm_start = _GridSeed(family, Variable, Target, grid, anchors, peak,
                               lambda model: _AsSos(model, family, center, reverse) is not None)
        if warm_start is not None:
            alternative = Exchange(warm_start, ExchangeError, Candidates, count, grid, tol=1e-9)
            _, alternate_peaks = Extrema(lambda f: ExchangeError(alternative, f), grid)
            if np.max(np.abs(alternate_peaks)) < peak:
                best = alternative

    check = np.unique(np.r_[np.linspace(0, 1, 65538), grid])
    _, initial_extrema = Extrema(lambda f: ExchangeError(initial, f), check)
    _, extrema = Extrema(lambda f: ExchangeError(best, f), check)
    peak = np.max(np.abs(extrema))
    initial_peak = np.max(np.abs(initial_extrema))
    if peak > initial_peak:
        best, extrema, peak = initial, initial_extrema, initial_peak
    sos = _AsSos(best, family, center, reverse)
    if sos is None:
        raise RuntimeError(f"Cannot factor {name} at {fc} Hz, {parameters}")
    selected = SelectReference(extrema, count)
    lower = np.min(np.abs(extrema[selected])) if selected is not None else 0.0
    # Check the cascade, including its declared anchors and preserved zeros.
    frequencies = np.r_[anchors, np.linspace(0.0001, 0.9999, 1024)]
    if family.notch:
        frequencies = frequencies[np.abs(frequencies - center) > 1e-9]
    x = PowerVariable(frequencies)
    y = Variable(frequencies)
    expected_db = 10 * np.log10(poly.polyval(y, best[0]) / poly.polyval(y, best[1]))
    if family.dc_zeros:
        expected_db += 10 * family.dc_zeros * np.log10(x)
    if family.notch:
        expected_db += 20 * np.log10(np.abs(1 - x / PowerVariable(center)))
    np.testing.assert_allclose(MagnitudeDb(sos, frequencies), expected_db, rtol=0, atol=2e-6)
    if len(anchors):
        np.testing.assert_allclose(Error(best, anchors), 0, rtol=0, atol=2e-7)
    zeros, poles, _ = signal.sos2zpk(sos)
    if family.notch:
        np.testing.assert_allclose(np.min(np.abs(zeros - np.exp(1j * np.pi * center))), 0,
                                   rtol=0, atol=1e-9)
    if family.dc_zeros and np.count_nonzero(np.abs(zeros - 1) < 1e-6) < family.dc_zeros:
        raise RuntimeError(f"Missing DC zero for {name}")
    radius = np.max(np.abs(poles))
    if not np.isfinite(peak) or radius >= 1:
        raise RuntimeError(f"Invalid SOS for {name} at {fc} Hz, {parameters}")
    return sos, initial_peak, peak, lower, radius, count


def Rows(families):
    yield ("family", "fc_hz", "q", "gain_db", "brasil_db", "remez_db",
           "alternation_lower_db", "reference_count", "max_pole", "sos")
    for name in families:
        module = importlib.import_module(f"{name}_interp")
        for fc, q, gain in product(module.FC_LIST, getattr(module, "Q_LIST", (None,)),
                                   getattr(module, "DB_LIST", (None,))):
            parameters = tuple(value for value in (q, gain) if value is not None)
            result = Design(name, fc, parameters)
            if result is None:
                raise RuntimeError(f"No usable seed for {name} at {fc} Hz, {parameters}")
            sos, initial, peak, lower, radius, count = result
            yield (name, fc, q, gain, initial, peak, lower, count,
                   radius, json.dumps(sos.tolist(), separators=(",", ":")))


def Main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", choices=FAMILIES, action="append")
    args = parser.parse_args()
    families = tuple(dict.fromkeys(args.family)) if args.family else tuple(FAMILIES)
    output = Path(__file__).resolve().parent / "output"
    output.mkdir(exist_ok=True)
    suffix = "_" + "_".join(families) if args.family else ""
    path = output / f"constrained{suffix}.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        csv.writer(stream).writerows(Rows(families))
    print(path)


if __name__ == "__main__":
    Main()
