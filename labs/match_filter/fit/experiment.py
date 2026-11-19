"""Compare fixed interpolation nodes and BRASIL on hand-specified dB curves."""

import argparse
from pathlib import Path

import numpy as np
from scipy import signal

from fit_curve import FitCurve, MagnitudeDb, RipplePeaks, SampledCurve


SAMPLE_RATE = 48000
TARGETS = (
    ("Broad contour",
     [0, 400, 1500, 3400, 6800, 12300, 19400, 24000],
     [2, 1, 3, 7, 3, -2, -4, -3]),
    ("Alternating bands",
     [0, 500, 1800, 4200, 7000, 10500, 15300, 19500, 24000],
     [0, 1, 5, -3, 4, 1, -5, 2, 0]),
    ("Low-frequency detail",
     [0, 40, 110, 320, 900, 2200, 7100, 13000, 24000],
     [1, 0, 5, 2, -3, 2, 4, -2, 0]),
    ("Narrow dip",
     [0, 2000, 6200, 9800, 10100, 10400, 14500, 20000, 24000],
     [0, 1, 2, 1, -22, 1, 3, 0, -1]),
)


def PlotResults(path, grid, rows):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(rows), 2, figsize=(11, 9), sharex=True)
    frequency = grid * (SAMPLE_RATE / 2)
    for (response_axis, error_axis), (name, anchor_hz, anchor_db, reference, fits) in zip(
            axes, rows):
        response_axis.plot(frequency, reference, color="black", label="Target")
        response_axis.scatter(anchor_hz, anchor_db, s=14, color="black")
        for sections, response in fits:
            response_axis.plot(frequency, response, label=f"{sections} SOS")
            line, = error_axis.plot(frequency, response - reference)
            peak = np.max(np.abs(response - reference))
            for sign in (-1, 1):
                error_axis.axhline(sign * peak, color=line.get_color(), alpha=0.28, linewidth=0.7)
        response_axis.set_title(name, loc="left")
        response_axis.set_ylabel("Magnitude (dB)")
        error_axis.axhline(0, color="black", linewidth=0.7)
        error_axis.set_ylabel("Signed error (dB)")
        for axis in (response_axis, error_axis):
            axis.set_xscale("symlog", linthresh=100)
            axis.set_xlim(0, SAMPLE_RATE / 2)
            axis.grid(alpha=0.25)
    axes[0, 0].legend(ncol=2)
    for axis in axes[-1]:
        axis.set_xlabel("Frequency (Hz)")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def RippleStats(error, expected):
    extrema = RipplePeaks(error)
    heights = np.abs(extrema)
    peak = heights.max()
    alternating = 0
    # Look for the closest-to-equal-height set with the required alternation.
    for height in np.unique(heights)[::-1]:
        signs = np.signbit(extrema[heights >= height])
        count = 1 + np.count_nonzero(np.diff(signs))
        if height >= 0.99 * peak:
            alternating = count
        if count >= expected:
            return peak - height, alternating
    return np.nan, alternating


def Main():
    parser = argparse.ArgumentParser(description=__doc__)
    default_plot = Path(__file__).resolve().parent / "output" / "curves.png"
    parser.add_argument("--plot", type=Path, nargs="?", const=default_plot,
                        help="write the response and error figure (default: output/curves.png)")
    args = parser.parse_args()

    # Use a separate, denser linear grid, plus a logarithmic grid to verify low frequencies.
    grid = np.unique(np.concatenate((np.linspace(0.0, 1.0, 65538),
                                     np.geomspace(1e-7, 1.0, 8193))))
    rows = []
    print(f"{'Target':22} {'SOS':>3} {'Fixed (dB)':>12} {'Fit (dB)':>12} "
          f"{'Spread (dB)':>12} {'Alternation':>12} {'Max |pole|':>12}")
    for name, frequency, gain in TARGETS:
        target = SampledCurve(frequency, gain, SAMPLE_RATE)
        reference = target(grid)
        fits = []
        for sections in (1, 2, 3):
            result = FitCurve(target, sections)
            if result is None:
                print(f"{name:22} {sections:3d}   no realizable fit")
                continue
            response = MagnitudeDb(result.sos, grid)
            initial_error = (np.max(np.abs(MagnitudeDb(result.initial_sos, grid) - reference))
                             if result.initial_sos is not None else np.nan)
            error = np.max(np.abs(response - reference))
            expected = 4 * sections + 2
            spread, alternating = RippleStats(response - reference, expected)
            spread_text = f"{spread:.4f}" if np.isfinite(spread) else "-"
            _, poles, _ = signal.sos2zpk(result.sos)
            radius = np.max(np.abs(poles))
            if radius >= 1.0 or not np.isfinite(error):
                raise RuntimeError(f"invalid SOS for {name}, {sections} sections")
            print(f"{name:22} {sections:3d} {initial_error:12.4f} {error:12.4f} "
                  f"{spread_text:>12} {f'{alternating}/{expected}':>12} {radius:12.8f}",
                  flush=True)
            fits.append((sections, response))
        rows.append((name, frequency, gain, reference, fits))
    if args.plot:
        PlotResults(args.plot, grid, rows)
        print(args.plot)


if __name__ == "__main__":
    Main()
