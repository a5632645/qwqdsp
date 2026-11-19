"""Fit a finite dB magnitude curve with free biquad coefficients."""

from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial import chebyshev as cheb
from scipy import optimize, signal
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "constraint_fit"))
from brasil import equalize, segment_maxima


@dataclass
class FitResult:
    sos: np.ndarray
    initial_sos: np.ndarray | None
    nodes: np.ndarray


def MagnitudeDb(sos, frequency):
    """Evaluate SOS at frequencies normalized to Nyquist."""
    _, response = signal.sosfreqz(sos, worN=np.pi * np.asarray(frequency))
    return 20.0 * np.log10(np.abs(response))


def _Crossings(grid, error):
    indices = np.flatnonzero(np.diff(np.signbit(error)))
    left, right = error[indices], error[indices + 1]
    return grid[indices] - left * np.diff(grid)[indices] / (right - left)


def RipplePeaks(error):
    """Return the signed extremum of each consecutive sign lobe."""
    cuts = np.flatnonzero(np.diff(np.signbit(error))) + 1
    return np.array([span[np.argmax(np.abs(span))] for span in np.split(error, cuts)])


def _SpectralRoots(coefficients, warp):
    """Factor a strictly positive polynomial in cos(w) on the unit circle."""
    sign = np.sign(cheb.chebval(0.0, coefficients))
    if sign == 0.0:
        return None
    roots = cheb.chebroots(coefficients).astype(complex)
    if warp != 1.0:
        # Undo the fractional-linear change of variable before factoring in z.
        divisor = (1.0 + warp) + (1.0 - warp) * roots
        if np.any(divisor == 0.0):
            return None
        roots = ((1.0 - warp) + (1.0 + warp) * roots) / divisor
    mapped = _UnitDiskRoots(roots)
    return (sign, mapped) if mapped is not None else None


def _UnitDiskRoots(cosine_roots):
    # t = (z + 1/z)/2. Compute the outer root first to avoid cancellation.
    cosine_roots = np.asarray(cosine_roots, dtype=complex)
    half_difference = np.sqrt(cosine_roots - 1.0) * np.sqrt(cosine_roots + 1.0)
    plus, minus = cosine_roots + half_difference, cosine_roots - half_difference
    outer = np.where(np.abs(plus) >= np.abs(minus), plus, minus)
    if np.any(np.abs(outer) <= 1.0 + 1e-12):
        return None
    return 1.0 / outer


def _Interpolate(nodes, target_db, order):
    desired_db = np.asarray(target_db(nodes), dtype=float)
    if desired_db.shape != nodes.shape or not np.all(np.isfinite(desired_db)):
        return None
    small_power = 10.0 ** (-np.abs(desired_db) / 10.0)
    x = np.sin(0.5 * np.pi * nodes) ** 2
    for warp in (1.0, 0.03, 0.001):
        sos = _SolveInterpolant(x, desired_db, small_power, nodes, order, warp)
        if sos is not None:
            return sos
    return None


def _SolveInterpolant(x, desired_db, small_power, nodes, order, warp):
    basis = cheb.chebvander(1.0 - 2.0 * x / (x + warp * (1.0 - x)), order)
    # Solve N(t) = target_power * D(t) without fixing a denominator coefficient.
    # Divide each row by max(1, target_power) to avoid overflow.
    numerator_scale = np.where(desired_db > 0.0, small_power, 1.0)
    denominator_scale = np.where(desired_db < 0.0, small_power, 1.0)
    matrix = np.hstack((basis * numerator_scale[:, None],
                        -basis * denominator_scale[:, None]))
    _, singular, vectors = np.linalg.svd(matrix, full_matrices=True)
    if singular[-1] <= 1e-13 * singular[0]:
        return None
    numerator, denominator = np.split(vectors[-1], 2)
    if warp != 1.0 and (numerator[-1] == 0.0 or denominator[-1] == 0.0):
        return None
    zeros = _SpectralRoots(numerator, warp)
    poles = _SpectralRoots(denominator, warp)
    if zeros is None or poles is None or zeros[0] != poles[0]:
        return None
    sos = signal.zpk2sos(zeros[1], poles[1], 1.0)
    actual_db = MagnitudeDb(sos, nodes)
    if not np.all(np.isfinite(actual_db)):
        return None
    gain_db = np.mean(desired_db - actual_db)
    if np.max(np.abs(actual_db + gain_db - desired_db)) > 1e-5:
        return None
    sos[0, :3] *= 10.0 ** (gain_db / 20.0)
    return sos


def _RefineNodes(nodes, measure, *, max_iter=30):
    """Use the sensitivity of the interval peaks to finish equalizing them."""
    peaks = measure(nodes)
    best_nodes, best_peak = nodes, peaks.max()
    for _ in range(max_iter):
        residual = np.diff(peaks)
        if np.ptp(peaks) < 1e-7:
            break
        gaps = np.diff(np.r_[0.0, nodes, 1.0])
        steps = 1e-3 * np.minimum(gaps[:-1], gaps[1:])
        jacobian = np.empty((nodes.size, nodes.size))
        for i, step in enumerate(steps):
            offset = np.zeros_like(nodes)
            offset[i] = step
            left, right = measure(nodes - offset), measure(nodes + offset)
            if left is None or right is None:
                return best_nodes, best_peak
            jacobian[:, i] = np.diff(right - left) / (2.0 * step)
        delta, _, rank, _ = np.linalg.lstsq(jacobian, -residual, rcond=1e-9)
        if rank < nodes.size:
            break
        for fraction in 0.5 ** np.arange(20):
            candidate = nodes + fraction * delta
            if np.min(np.diff(np.r_[0.0, candidate, 1.0])) <= 1e-12:
                continue
            candidate_peaks = measure(candidate)
            if (candidate_peaks is not None
                    and np.linalg.norm(np.diff(candidate_peaks)) < np.linalg.norm(residual)):
                nodes, peaks = candidate, candidate_peaks
                if peaks.max() < best_peak:
                    best_nodes, best_peak = nodes, peaks.max()
                break
        else:
            break
    return best_nodes, best_peak


def _MinimaxSeed(grid, reference, order, upper_error, warp):
    """Bound the sample-wise dB error by bisection over linear feasibility."""
    x = np.sin(0.5 * np.pi * grid) ** 2
    basis = cheb.chebvander(1.0 - 2.0 * x / (x + warp * (1.0 - x)), order)
    if np.max(np.abs(reference)) + upper_error > 3000.0:
        return None
    power = 10.0 ** (reference / 10.0)
    active = np.unique(np.r_[np.arange(0, grid.size, 32), grid.size - 1])
    lower_error = 0.0
    best_sos = None
    for _ in range(24):
        limit = (lower_error + upper_error) / 2.0
        ratio = 10.0 ** (limit / 10.0)
        rows = basis[active]
        low = power[active] / ratio
        high = power[active] * ratio
        zeros = np.zeros_like(rows)
        lhs = np.vstack((np.c_[rows, -high[:, None] * rows[:, 1:]],
                         np.c_[-rows, low[:, None] * rows[:, 1:]],
                         np.c_[zeros, -rows[:, 1:]]))
        rhs = np.r_[high, -low, np.full(active.size, 1.0 - 1e-9)]
        result = optimize.linprog(np.zeros(2 * order + 1), A_ub=lhs, b_ub=rhs,
                                  bounds=(None, None), method="highs",
                                  options={"primal_feasibility_tolerance": 1e-9})
        if result.status == 2:
            lower_error = limit
            continue
        if not result.success:
            break
        numerator = result.x[:order + 1]
        denominator = np.r_[1.0, result.x[order + 1:]]
        num, den = basis @ numerator, basis @ denominator
        with np.errstate(divide="ignore", invalid="ignore"):
            signed_error = 10.0 * (np.log10(num) - np.log10(den)) - reference
        error = np.abs(signed_error)
        error = np.where(np.isfinite(error), error, np.inf)
        if error.max() > limit + 5e-5:
            worst = []
            for band in np.array_split(np.arange(grid.size), min(grid.size, 32)):
                invalid = band[~np.isfinite(signed_error[band])]
                if invalid.size:
                    worst.append(invalid[0])
                else:
                    worst.extend((band[np.argmin(signed_error[band])],
                                  band[np.argmax(signed_error[band])]))
            augmented = np.union1d(active, worst)
            if augmented.size == active.size:
                break
            active = augmented
            continue
        zeros, poles = _SpectralRoots(numerator, warp), _SpectralRoots(denominator, warp)
        if zeros is None or poles is None or zeros[0] != 1 or poles[0] != 1:
            break
        if warp != 1.0 and (zeros[1].size != order or poles[1].size != order):
            break
        sos = signal.zpk2sos(zeros[1], poles[1], 1.0)
        mid = grid.size // 2
        gain_db = 10.0 * np.log10(num[mid] / den[mid]) - MagnitudeDb(sos, grid[mid])[0]
        sos[0, :3] *= 10.0 ** (gain_db / 20.0)
        actual_error = np.max(np.abs(MagnitudeDb(sos, grid) - reference))
        if not np.isfinite(actual_error) or actual_error >= upper_error:
            break
        upper_error, best_sos = actual_error, sos
        if upper_error - lower_error < 1e-5:
            break
    return best_sos


def FitCurve(target_db, sections, *, grid_size=8193, starts=4, max_iter=400):
    """Fit a callable returning finite dB values at normalized frequencies [0, 1].

    Return a minimum-phase SOS, its optional fixed-node baseline, and the
    estimated error crossings. Return None if no realizable SOS was found.
    """
    if sections < 1 or grid_size < 4 * sections + 2 or starts < 1:
        raise ValueError("need positive sections and starts, and enough grid points")
    order = 2 * sections
    count = 2 * order + 1
    grid = np.unique(np.r_[np.linspace(0.0, 1.0, grid_size),
                           np.geomspace(1e-6, 1.0, grid_size // 4)])
    reference = np.asarray(target_db(grid), dtype=float)
    if reference.shape != grid.shape or not np.all(np.isfinite(reference)):
        raise ValueError("target_db must return a finite dB value per frequency")

    variation = np.concatenate(([0.0], np.cumsum(np.abs(np.diff(reference)))))
    if variation[-1] == 0.0:
        sos = np.tile([1.0, 0.0, 0.0, 1.0, 0.0, 0.0], (sections, 1))
        sos[0, 0] = 10.0 ** (reference[0] / 20.0)
        return FitResult(sos, sos.copy(), np.empty(0))

    # Reserve more interval length for features that a linear grid compresses.
    coordinate = 0.25 * grid + 0.75 * variation / variation[-1]

    def Measure(nodes):
        sos = _Interpolate(np.interp(nodes, coordinate, grid), target_db, order)
        if sos is None:
            return None
        error = MagnitudeDb(sos, grid) - reference
        if not np.all(np.isfinite(error)):
            return None
        return segment_maxima(coordinate, np.abs(error), nodes)

    regular = (np.arange(count) + 0.5) / count
    frequency_seeds = [regular, (1.0 - np.cos(np.pi * regular)) / 2.0]
    for exponent in (1.5, 2.0, 3.0, 4.0):
        frequency_seeds.extend((regular ** exponent, 1.0 - regular[::-1] ** exponent))
    seeds = [regular, *(np.interp(seed, grid, coordinate) for seed in frequency_seeds)]
    rng = np.random.default_rng(17)
    for _ in range(48):
        random_nodes = np.sort(rng.uniform(0.0, 1.0, count))
        seeds.extend((random_nodes, np.interp(random_nodes, grid, coordinate)))

    feasible = []
    for nodes in seeds:
        peaks = Measure(nodes)
        if peaks is not None:
            feasible.append((float(peaks.max()), nodes))
    feasible.sort(key=lambda entry: entry[0])

    best_error = feasible[0][0] if feasible else np.ptp(reference)
    initial_sos = (_Interpolate(np.interp(feasible[0][1], coordinate, grid), target_db, order)
                   if feasible else None)
    best_sos = initial_sos
    if best_error > 1e-8:
        for _, nodes in feasible[:starts]:
            fitted_nodes, _, error = equalize(nodes, Measure, max_iter=max_iter)
            refined_nodes, refined_error = _RefineNodes(fitted_nodes, Measure)
            if refined_error < error:
                fitted_nodes, error = refined_nodes, refined_error
            if error < best_error:
                best_error = error
                best_sos = _Interpolate(np.interp(fitted_nodes, coordinate, grid), target_db, order)
        for warp in (1.0, 0.03, 0.01):
            if best_sos is not None:
                peaks = RipplePeaks(MagnitudeDb(best_sos, grid) - reference)
                near_peak = peaks[np.abs(peaks) >= best_error - max(1e-5, 2e-4 * best_error)]
                if near_peak.size and 1 + np.count_nonzero(np.diff(np.signbit(near_peak))) >= count + 1:
                    break
            seed = _MinimaxSeed(grid, reference, order, best_error, warp)
            if seed is None:
                continue
            signed_error = MagnitudeDb(seed, grid) - reference
            crossings = _Crossings(grid, signed_error)
            error = np.max(np.abs(signed_error))
            if crossings.size == count and np.min(np.diff(np.r_[0.0, crossings, 1.0])) > 1e-12:
                nodes = np.interp(crossings, grid, coordinate)
                if Measure(nodes) is not None:
                    nodes, _, node_error = equalize(nodes, Measure, max_iter=max_iter)
                    nodes, node_error = _RefineNodes(nodes, Measure)
                    if node_error < error:
                        error = node_error
                        seed = _Interpolate(np.interp(nodes, coordinate, grid), target_db, order)
            if error < best_error:
                best_error, best_sos = error, seed
    if best_sos is None:
        return None
    crossings = _Crossings(grid, MagnitudeDb(best_sos, grid) - reference)
    return FitResult(best_sos, initial_sos, crossings)


def SampledCurve(frequency_hz, gain_db, sample_rate):
    """Interpolate tabulated dB values spanning DC through Nyquist."""
    if not np.isfinite(sample_rate) or sample_rate <= 0:
        raise ValueError("sample_rate must be positive and finite")
    frequency = np.asarray(frequency_hz, dtype=float) / (sample_rate / 2.0)
    if frequency.ndim != 1 or frequency.size < 2 or frequency[0] != 0 or frequency[-1] != 1:
        raise ValueError("frequency samples must span DC through Nyquist")
    return PchipInterpolator(frequency, gain_db, extrapolate=False)
