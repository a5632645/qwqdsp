"""Compare the BRASIL experiment with independent nonlinear SOS least-p fitting."""

import argparse
import csv
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy import optimize, signal

from experiment import SAMPLE_RATE, TARGETS, RippleStats
from fit_curve import FitCurve, MagnitudeDb, SampledCurve


KDb = 20.0 / np.log(10.0)
KReflectionLimit = 1.0 - 1e-8


def _Coefficients(reflection):
    return np.stack((np.ones_like(reflection[..., 0]),
                     reflection[..., 0] * (1.0 + reflection[..., 1]),
                     reflection[..., 1]), axis=-1)


def _ToSos(parameters):
    reflection = parameters[1:].reshape(-1, 2, 2)
    sos = _Coefficients(reflection).reshape(-1, 6)
    sos[0, :3] *= 10.0 ** (parameters[0] / 20.0)
    return sos


def _Seed(frequency, width, sign, scale):
    # Start with nearby pole/zero pairs and allow all coefficients to move.
    radius = np.exp(-np.pi * width * scale * np.array([0.6, 1.4]))
    if sign > 0:
        radius = radius[::-1]
    second = radius ** 2
    first = -2.0 * radius * np.cos(np.pi * frequency) / (1.0 + second)
    return np.column_stack((first, second))


def _InitialParameters(grid, target_db, sections):
    features = []
    for sign in (-1, 1):
        indices, properties = signal.find_peaks(sign * target_db,
                                                prominence=0.05, width=1)
        samples = np.arange(grid.size)
        left = np.interp(properties["left_ips"], samples, grid)
        right = np.interp(properties["right_ips"], samples, grid)
        features.extend(zip(properties["prominences"], grid[indices], right - left,
                            np.full(indices.size, sign)))
    features.sort(reverse=True)
    regular = [(0.0, f, 0.15, (-1) ** j)
               for j, f in enumerate((np.arange(sections) + 0.5) / sections)]
    middle = np.median(target_db)
    for shift in range(min(len(features), 4) + 1):
        chosen = (features[shift:] + features[:shift] + regular)[:sections]
        for scale in (0.6, 1.4):
            seed = [_Seed(frequency, width, sign, scale)
                    for _, frequency, width, sign in chosen]
            yield np.r_[middle, np.clip(np.ravel(seed), -KReflectionLimit, KReflectionLimit)]
    yield np.r_[middle, np.zeros(4 * sections)]


def FitLeastP(target, sections):
    """Fit stable SOS in dB with independent multistart and p=2 through p=128."""
    grid = np.unique(np.r_[np.linspace(0.0, 1.0, 769), np.geomspace(1e-6, 1.0, 257)])
    target_db = target(grid)
    widths = np.diff(grid)
    weights = 0.5 * (np.r_[widths, 0.0] + np.r_[0.0, widths])
    root_weights = np.sqrt(weights)
    z = np.exp(-1j * np.pi * grid)[:, None, None]
    z2 = z * z
    signs = np.array([1.0, -1.0])
    lower = np.r_[-np.inf, np.full(4 * sections, -KReflectionLimit)]
    upper = np.r_[np.inf, np.full(4 * sections, KReflectionLimit)]

    def Evaluate(parameters):
        reflection = parameters[1:].reshape(sections, 2, 2)
        first, second = reflection[..., 0], reflection[..., 1]
        response = 1.0 + first * (1.0 + second) * z + second * z2
        error = parameters[0] + KDb * (np.log(np.abs(response)) * signs).sum((1, 2)) - target_db
        d_first = KDb * signs * np.real((1.0 + second) * z / response)
        d_second = KDb * signs * np.real((first * z + z2) / response)
        jacobian = np.column_stack((np.ones(grid.size),
                                    np.stack((d_first, d_second), axis=-1).reshape(grid.size, -1)))
        return error, jacobian

    def Solve(parameters, power):
        scale = max(1e-8, np.abs(Evaluate(parameters)[0]).max())

        def Residual(candidate):
            error, _ = Evaluate(candidate)
            normalized = error / scale
            return root_weights * normalized * np.abs(normalized) ** (power / 2 - 1)

        def Jacobian(candidate):
            error, jacobian = Evaluate(candidate)
            factor = (power / 2) * (np.abs(error) / scale) ** (power / 2 - 1) / scale
            return jacobian * (root_weights * factor)[:, None]

        result = optimize.least_squares(Residual, parameters, jac=Jacobian,
                                        bounds=(lower, upper), max_nfev=120)
        error, _ = Evaluate(result.x)
        return np.max(np.abs(error)), np.dot(weights, error ** 2), result.x

    starts = [Solve(seed, 2) for seed in _InitialParameters(grid, target_db, sections)]
    best_rms = min(starts, key=lambda row: row[1])
    starts.sort(key=lambda row: row[0])
    best_peak = starts[0]
    for _, _, parameters in starts[:3]:
        for power in (4, 8, 16, 32, 64, 128):
            row = Solve(parameters, power)
            parameters = row[2]
            if row[0] < best_peak[0]:
                best_peak = row
    return _ToSos(best_rms[2]), _ToSos(best_peak[2])


def Rows(target_index, sections_count):
    grid = np.unique(np.r_[np.linspace(0.0, 1.0, 65538), np.geomspace(1e-7, 1.0, 8193)])
    yield "Target", "SOS", "Method", "Peak dB", "RMS dB", "Alternation", "Max pole", "Seconds"
    for index, (name, frequency, gain) in enumerate(TARGETS, start=1):
        if target_index is not None and index != target_index:
            continue
        target = SampledCurve(frequency, gain, SAMPLE_RATE)
        reference = target(grid)
        for sections in ((sections_count,) if sections_count is not None else (1, 2, 3)):
            start = perf_counter()
            least_squares, least_p = FitLeastP(target, sections)
            least_p_time = perf_counter() - start
            start = perf_counter()
            brasil = FitCurve(target, sections)
            rows = [("least-squares", least_squares, ""),
                    ("least-p", least_p, f"{least_p_time:.2f}")]
            if brasil is not None:
                rows.append(("BRASIL + grid", brasil.sos, f"{perf_counter() - start:.2f}"))
            for method, sos, elapsed in rows:
                error = MagnitudeDb(sos, grid) - reference
                radius = np.abs(signal.sos2zpk(sos)[1]).max()
                if radius >= 1.0 or not np.all(np.isfinite(error)):
                    raise RuntimeError(f"Invalid fit: {name}, {sections} SOS, {method}")
                _, alternating = RippleStats(error, 4 * sections + 2)
                rms = np.sqrt(np.trapezoid(error ** 2, grid))
                yield (name, sections, method, f"{np.abs(error).max():.5f}", f"{rms:.5f}",
                       f"{alternating}/{4 * sections + 2}", f"{radius:.8f}", elapsed)


def Main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=int, choices=range(1, len(TARGETS) + 1))
    parser.add_argument("--sos", type=int, choices=(1, 2, 3))
    args = parser.parse_args()
    output = Path(__file__).resolve().parent / "output"
    output.mkdir(exist_ok=True)
    target_suffix = f"_target_{args.target}" if args.target is not None else ""
    sos_suffix = f"_sos_{args.sos}" if args.sos is not None else ""
    path = output / f"industry_comparison{target_suffix}{sos_suffix}.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        csv.writer(stream).writerows(Rows(args.target, args.sos))
    print(path)


if __name__ == "__main__":
    Main()
