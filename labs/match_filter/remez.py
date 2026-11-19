"""Scalar leveling and extremum exchange for rational dB minimax fits."""

import numpy as np
from scipy import optimize


def Level(matrix, ceiling, *, samples=65):
    """Yield the null vectors of the leveling system as the dB ripple varies."""
    def Scaled(ripple):
        coefficients = matrix(ripple)
        scales = np.maximum(np.linalg.norm(coefficients, axis=0), 1e-300)
        return coefficients / scales, scales

    def Determinant(ripple):
        sign, log_det = np.linalg.slogdet(Scaled(ripple)[0])
        return sign * np.exp(log_det)

    levels = np.linspace(0.0, ceiling * 1.001, samples)
    determinants = [Determinant(ripple) for ripple in levels]
    for i, (left, right) in enumerate(zip(determinants[:-1], determinants[1:])):
        if left == right == 0.0:
            continue
        if left != 0.0 and right != 0.0 and np.signbit(left) == np.signbit(right):
            continue
        ripple = optimize.brentq(Determinant, levels[i], levels[i + 1], xtol=1e-14)
        scaled, scales = Scaled(ripple)
        _, _, vectors = np.linalg.svd(scaled, full_matrices=False)
        yield ripple, vectors[-1] / scales


def Extrema(error, grid, signed=None):
    """Find the largest signed extremum of each sign lobe."""
    if signed is None:
        signed = error(grid)
    starts = np.r_[0, np.flatnonzero(np.diff(np.signbit(signed))) + 1, grid.size]
    indices = [a + np.argmax(np.abs(signed[a:b])) for a, b in zip(starts[:-1], starts[1:])]
    points, values = grid[indices].copy(), signed[indices].copy()
    for k, i in enumerate(indices):
        if i == 0 or i == grid.size - 1:
            continue
        sign = np.sign(values[k])

        def Objective(frequency):
            return -sign * error(np.array([frequency]))[0]

        refined = optimize.minimize_scalar(Objective, bounds=(grid[i - 1], grid[i + 1]),
                                           method="bounded", options={"xatol": 1e-14})
        if -refined.fun > abs(values[k]):
            points[k], values[k] = refined.x, -sign * refined.fun
    return points, values


def SelectReference(values, count):
    """Select the strongest alternating set including the global extremum."""
    heights = np.abs(values)
    global_index = np.argmax(heights)
    signs = np.signbit(values)
    for threshold in np.unique(heights)[::-1]:
        kept = []
        for i in np.flatnonzero(heights >= threshold):
            if not kept or signs[i] != signs[kept[-1]]:
                kept.append(i)
            elif heights[i] > heights[kept[-1]]:
                kept[-1] = i
        if len(kept) < count:
            continue
        possible = [kept[i:i + count] for i in range(len(kept) - count + 1)
                    if global_index in kept[i:i + count]]
        if possible:
            return max(possible, key=lambda seq: heights[seq].sum())
    return None


def Exchange(initial, error, level, count, grid, *, tol=1e-8, max_iter=25):
    """Exchange extrema and keep the feasible iterate with the lowest peak."""
    model, best = initial, initial
    best_peak = np.inf
    for iteration in range(max_iter + 1):
        signed = error(model, grid)
        grid_peak = np.max(np.abs(signed))
        if not np.isfinite(grid_peak):
            break
        if grid_peak <= tol:
            return model
        points, values = Extrema(lambda f: error(model, f), grid, signed)
        peak = np.max(np.abs(values))
        if peak < best_peak:
            best, best_peak = model, peak
        chosen = SelectReference(values, count)
        if (iteration == max_iter or chosen is None
                or peak - np.min(np.abs(values[chosen])) <= tol):
            break
        candidate = None
        for fit in level(points[chosen], np.sign(values[chosen]), peak):
            score = np.max(np.abs(error(fit, grid)))
            if np.isfinite(score) and (candidate is None or score < candidate[0]):
                candidate = score, fit
        if candidate is None or candidate[0] > 2.0 * peak:
            break
        model = candidate[1]
    return best
