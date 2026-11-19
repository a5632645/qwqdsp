"""Run the free-biquad rational Remez benchmark."""

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy import signal

from experiment import SAMPLE_RATE, TARGETS
from fit_curve import SampledCurve
from remez_fit import FitRemez, RemezBounds


def Rows(target_index):
    grid = np.unique(np.r_[np.linspace(0, 1, 65538), np.geomspace(1e-7, 1, 8193)])
    yield ("target", "sections", "peak_db", "alternation_lower_db", "gap_db", "max_pole")
    for index, (name, hz, db) in enumerate(TARGETS, 1):
        if target_index is not None and target_index != index:
            continue
        target = SampledCurve(hz, db, SAMPLE_RATE)
        for sections in (1, 2, 3):
            sos = FitRemez(target, sections)
            peak, lower = RemezBounds(target, sos, grid)
            radius = np.max(np.abs(signal.sos2zpk(sos)[1]))
            if len(sos) != sections or radius >= 1 or not np.isfinite(peak):
                raise RuntimeError(f"invalid SOS for {name}, {sections} sections")
            yield name, sections, f"{peak:.9f}", f"{lower:.9f}", f"{peak-lower:.9f}", f"{radius:.9f}"


def Main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=int, choices=range(1, len(TARGETS) + 1))
    args = parser.parse_args()
    output = Path(__file__).resolve().parent / "output"
    output.mkdir(exist_ok=True)
    suffix = f"_target_{args.target}" if args.target is not None else ""
    path = output / f"remez{suffix}.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        csv.writer(stream).writerows(Rows(args.target))
    print(path)


if __name__ == "__main__":
    Main()
