"""Rational Remez fit for free-coefficient, minimum-phase digital SOS."""

from pathlib import Path
import sys

import numpy as np
from numpy.polynomial import chebyshev as cheb
from scipy import signal

from fit_curve import FitCurve, MagnitudeDb, _SpectralRoots

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from remez import Exchange, Extrema, Level, SelectReference


def _LevelSos(frequency, target_db, signs, order, warp, ceiling):
    x = np.sin(0.5 * np.pi * frequency) ** 2
    basis = cheb.chebvander(1.0 - 2.0 * x / (x + warp * (1.0 - x)), order)

    def Matrix(ripple):
        db = target_db + ripple * signs
        scale = 10.0 ** (-np.abs(db) / 10.0)
        return np.c_[basis * np.where(db > 0.0, scale, 1.0)[:, None],
                     -basis * np.where(db < 0.0, scale, 1.0)[:, None]]

    for ripple, vector in Level(Matrix, ceiling):
        numerator, denominator = np.split(vector, 2)
        zeros, poles = _SpectralRoots(numerator, warp), _SpectralRoots(denominator, warp)
        if zeros is None or poles is None or zeros[0] != poles[0]:
            continue
        if warp != 1.0 and (zeros[1].size != order or poles[1].size != order):
            continue
        sos = signal.zpk2sos(zeros[1], poles[1], 1.0)
        gain_db = target_db + ripple * signs - MagnitudeDb(sos, frequency)
        if np.all(np.isfinite(gain_db)) and np.ptp(gain_db) <= 1e-6:
            sos[0, :3] *= 10.0 ** (gain_db.mean() / 20.0)
            yield sos


def RefineRemez(target_db, initial_sos, grid):
    """Refine a feasible SOS seed on a normalized DC-to-Nyquist grid."""
    order = 2 * len(initial_sos)
    reference = target_db(grid)

    def Error(sos, frequency):
        target = reference if frequency is grid else target_db(frequency)
        return MagnitudeDb(sos, frequency) - target

    def Candidates(frequency, signs, peak):
        target = target_db(frequency)
        for warp in (1.0, 0.03, 0.01, 0.001):
            yield from _LevelSos(frequency, target, signs, order, warp, peak)

    return Exchange(initial_sos, Error, Candidates, 2 * order + 2, grid)


def FitRemez(target_db, sections, *, grid_size=8193):
    """Fit a finite dB callable on [0, 1] and return a minimum-phase SOS."""
    from compare_lpnorm import FitLeastP

    baseline = FitCurve(target_db, sections, grid_size=grid_size)
    grid = np.unique(np.r_[np.linspace(0, 1, grid_size),
                           np.geomspace(1e-6, 1.0, grid_size // 4)])
    candidates = []
    if baseline is not None:
        candidates.append(RefineRemez(target_db, baseline.sos, grid))
        peak, lower = RemezBounds(target_db, candidates[0], grid)
        if peak - lower <= 1e-8:
            return candidates[0]
    _, least_p = FitLeastP(target_db, sections)
    candidates.append(RefineRemez(target_db, least_p, grid))
    reference = target_db(grid)
    return min(candidates, key=lambda sos: np.max(np.abs(MagnitudeDb(sos, grid) - reference)))


def RemezBounds(target_db, sos, grid):
    """Return the sampled/refined peak and the alternating-error lower bound."""
    _, values = Extrema(lambda frequency: MagnitudeDb(sos, frequency) - target_db(frequency),
                        grid)
    chosen = SelectReference(values, 4 * len(sos) + 2)
    lower = np.min(np.abs(values[chosen])) if chosen is not None else 0.0
    return np.max(np.abs(values)), lower
