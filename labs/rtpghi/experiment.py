"""Compare causal and bidirectional phase integration on the same STFT."""

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import wavfile

from rtpghi import (
    Stft, analyze, anchor_frame_range, causal, detect_onsets, phase_anchors, phase_edges, propagate,
    phase_loop_residual, propagate_global, wola,
)


METHODS = ("direct", "centered", "reset", "forward", "bidirectional")
WORKSPACE = Path(__file__).resolve().parents[3]


@dataclass(frozen=True)
class Case:
    name: str
    signal: np.ndarray
    onsets: tuple[int, ...] | None
    duration: int = 1  # Intended transient support in samples at the output.


@dataclass(frozen=True)
class RenderResult:
    outputs: dict
    weak_seeds: dict
    consistency: dict
    anchor_bins: np.ndarray
    local_coverage: np.ndarray
    global_coverage: np.ndarray | None
    heap_metrics: dict
    propagation_scheme: str
    loop_metrics: dict


def impulse_case(rate, offset, name):
    signal = np.zeros(round(0.46 * rate))
    onset = round(0.22 * rate) + offset
    if not 0 <= onset < len(signal):
        raise ValueError("impulse offset puts the event outside the test signal")
    signal[onset] = 0.8
    return Case(name, signal, (onset,))


def builtin_cases(rate):
    impulse = impulse_case(rate, 37, "impulse")
    length = len(impulse.signal)
    onset, = impulse.onsets

    duration = round(0.007 * rate)
    rng = np.random.default_rng(42)
    burst = np.zeros(length)
    burst[onset:onset + duration] = 0.32 * rng.uniform(-1, 1, duration) * np.exp(-5 * np.arange(duration) / duration)
    burst[onset] = 0.8

    t = np.arange(length) / rate
    tone = 0.09 * np.sin(2 * np.pi * 441 * t) + 0.04 * np.sin(2 * np.pi * 882 * t + 0.3)
    tone *= np.minimum(1, np.minimum(t, len(t) / rate - t) / 0.03)
    second = onset + round(0.020 * rate)
    mixed = tone + 0.85 * burst
    mixed[second:second + duration] += 0.6 * burst[onset:onset + duration]
    return {
        "impulse": impulse,
        "burst": Case("burst", burst, (onset,), duration),
        "mixed": Case("mixed", mixed, (onset, second), duration),
        "tone": Case("tone", tone, ()),
    }


def input_wav(path):
    rate, signal = wavfile.read(path)
    if np.issubdtype(signal.dtype, np.unsignedinteger):
        midpoint = (np.iinfo(signal.dtype).max + 1) // 2
        signal = (signal.astype(float) - midpoint) / midpoint
    elif np.issubdtype(signal.dtype, np.signedinteger):
        signal = signal.astype(float) / -int(np.iinfo(signal.dtype).min)
    else:
        signal = signal.astype(float)
    if signal.ndim == 2:
        signal = signal.mean(axis=1)
    if signal.size == 0 or not np.all(np.isfinite(signal)):
        raise ValueError("expected a nonempty WAV with finite audio samples")
    return rate, signal


def measure_loops(magnitude, edges, tolerance):
    strength = np.minimum(np.minimum(magnitude[:-1, :-1], magnitude[:-1, 1:]),
                          np.minimum(magnitude[1:, :-1], magnitude[1:, 1:]))
    peak = magnitude.max()
    valid = strength > tolerance * peak
    if not np.any(valid):
        return dict(loop_cells=0, loop_rms_rad=None, loop_nonclosing_fraction=None)
    residual = phase_loop_residual(edges)[valid]
    weight = (strength[valid] / peak) ** 2
    return dict(loop_cells=len(residual), loop_rms_rad=float(np.sqrt(np.average(residual ** 2, weights=weight))),
                loop_nonclosing_fraction=float(np.mean(np.abs(residual) > 1e-8)))


def render(case, rate, config, onsets, radius, tolerance, *, global_heap=False, propagation_scheme="centered"):
    length = max(1, round(len(case.signal) * config.ratio))
    rows = max(2, 1 + (length + config.synthesis_hop - 1) // config.synthesis_hop)
    centers = np.arange(rows) * config.analysis_hop
    analysis = analyze(case.signal, config.window_size, config.fft_size, centers)
    magnitude, phase = np.abs(analysis), np.angle(analysis)
    edge_sets = {"direct": phase_edges(phase, config, centered=False),
                 "centered": phase_edges(phase, config, centered=True)}
    propagation = edge_sets[propagation_scheme]
    loop_metrics = {scheme: measure_loops(magnitude, edges, tolerance) for scheme, edges in edge_sets.items()}
    empty = np.full(phase.shape, np.nan)
    anchors = phase_anchors(case.signal, analysis, config, onsets, rate, tolerance)
    forward = propagate(magnitude, propagation, anchors, radius, tolerance, backward=False)
    bidirectional = propagate(magnitude, propagation, anchors, radius, tolerance, backward=True)
    fixed = (empty, empty, anchors, forward, bidirectional)
    outputs = {}
    weak_seeds = {}
    consistency = {}

    def synthesize(method, estimate):
        target = magnitude * np.exp(1j * estimate)
        target[:, [0, -1]] = target[:, [0, -1]].real
        outputs[method] = wola(target, config, length)
        observed = analyze(outputs[method], config.window_size, config.fft_size,
                           np.arange(rows) * config.synthesis_hop)
        consistency[method] = float(np.linalg.norm(observed - target) / max(1e-30, np.linalg.norm(target)))

    for method, locked in zip(METHODS, fixed):
        estimate, weak_seeds[method] = causal(magnitude, phase, edge_sets.get(method, propagation), locked, tolerance)
        synthesize(method, estimate)

    anchor_bins = np.count_nonzero(np.isfinite(anchors), axis=1)
    local_coverage = np.count_nonzero(np.isfinite(bidirectional), axis=1)
    global_coverage = None
    heap_metrics = {}
    if global_heap:
        anchor_range = anchor_frame_range(anchor_bins > 0, radius)
        methods = (
            ("global-magnitude", empty, False), ("global", empty, True),
            ("global-magnitude-anchored", anchors, False), ("global-anchored", anchors, True),
        )
        for method, locked, bottleneck in methods:
            integrated = propagate_global(magnitude, phase, propagation, locked, tolerance, bottleneck=bottleneck)
            known = np.isfinite(integrated.phase)
            estimate = np.random.default_rng(0).uniform(-np.pi, np.pi, magnitude.shape)
            estimate[known] = integrated.phase[known]
            synthesize(method, estimate)
            coverage = integrated.anchored_by_frame
            heap_metrics[method] = dict(heap_seeds=integrated.seed_count, heap_integrated_bins=int(known.sum()),
                                        heap_anchor_bins=int(coverage.sum()),
                                        heap_anchor_bins_outside_radius=int(coverage[~anchor_range].sum()))
            if method == "global-anchored":
                global_coverage = coverage
    return RenderResult(outputs, weak_seeds, consistency, anchor_bins, local_coverage, global_coverage,
                        heap_metrics, propagation_scheme, loop_metrics)


def energy_db(numerator, denominator):
    return 10 * np.log10(max(1e-15, numerator / max(1e-30, denominator)))


def measure(result, case, config, targets):
    records = []
    for method, output in result.outputs.items():
        weak_seeds = result.weak_seeds.get(method)
        scheme = method if method in result.loop_metrics else result.propagation_scheme
        for onset in targets or (None,):
            record = dict(method=method, edge_scheme=scheme, onset_input=onset, onset_output=None, pre_db=None, post_db=None,
                          peak_offset_samples=None, main_peak=None, weak_time_seeds=None,
                          stft_residual=result.consistency[method], heap_seeds=None, heap_integrated_bins=None,
                          heap_anchor_bins=None, heap_anchor_bins_outside_radius=None, **result.loop_metrics[scheme])
            record.update(result.heap_metrics.get(method, {}))
            if onset is not None:
                target = min(len(output) - 1, round(onset * config.ratio))
                index = targets.index(onset)
                # Clip at neighboring onset midpoints so another attack cannot
                # be counted as this event's pre/post energy.
                left = target - config.window_size // 2
                right = target + case.duration + config.window_size // 2
                if index:
                    left = max(left, (target + round(targets[index - 1] * config.ratio)) // 2)
                if index + 1 < len(targets):
                    right = min(right, (target + round(targets[index + 1] * config.ratio)) // 2)
                left, right = max(0, min(left, target)), min(len(output), max(right, target + 1))
                end = min(right, target + case.duration)
                local = output[left:right]
                total = np.dot(local, local)
                pre = output[left:target]
                post = output[end:right]
                row = round(onset / config.analysis_hop)
                near = 1 + config.window_size // (2 * config.analysis_hop)
                record.update(onset_output=target, pre_db=energy_db(np.dot(pre, pre), total),
                              post_db=energy_db(np.dot(post, post), total),
                              peak_offset_samples=int(np.argmax(np.abs(local))) + left - target,
                              main_peak=float(np.max(np.abs(output[target:end]))))
                if weak_seeds is not None:
                    record["weak_time_seeds"] = int(weak_seeds[max(0, row - near):row + near + 1].sum())
            records.append(record)
    return records


def plot_results(path, case, rate, config, result, targets):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (waveform, graph) = plt.subplots(2, 1, figsize=(11, 6), constrained_layout=True)
    length = len(next(iter(result.outputs.values())))
    marked = [min(length - 1, round(onset * config.ratio)) for onset in targets[:2]]
    origin = marked[0] if marked else length // 2
    reference = "first onset" if targets else "segment midpoint"
    left = max(0, origin - config.window_size // 2)
    last = marked[-1] if marked else origin
    right = min(length, last + config.window_size // 2 + case.duration)
    ms = (np.arange(left, right) - origin) * (1000 / rate)
    for method, signal in result.outputs.items():
        if result.global_coverage is not None and method in ("direct", "reset", "forward"):
            continue
        waveform.plot(ms, signal[left:right], label=method, linewidth=1, alpha=0.8)
    for onset in marked:
        waveform.axvspan(1000 * (onset - origin) / rate, 1000 * (onset - origin + case.duration) / rate,
                         color="0.8", alpha=0.6)
    waveform.set(xlabel=f"Output time relative to {reference} (ms)", ylabel="Amplitude",
                 title=f"{case.name} · time scale {config.ratio:.4g} · {result.propagation_scheme} propagation edges")
    waveform.legend(ncol=4, loc="upper right")

    reached = result.local_coverage
    locked = result.anchor_bins
    time = (np.arange(len(reached)) * config.synthesis_hop - origin) * (1000 / rate)
    graph.step(time, reached - locked, where="mid", label="Local anchor propagation")
    if result.global_coverage is not None:
        graph.step(time, result.global_coverage - locked, where="mid", label="Global anchor propagation", color="C2")
        reached = result.global_coverage
    if locked.any():
        graph.stem(time[locked > 0], locked[locked > 0], linefmt="C1-", markerfmt="C1o", basefmt=" ",
                   label="Fixed phase bins")
        active = time[reached > 0]
        margin = 1000 * config.synthesis_hop / rate
        graph.set_xlim(active[0] - margin, active[-1] + margin)
    elif len(ms) > 1:
        graph.set_xlim(ms[0], ms[-1])
    graph.set(xlabel=f"Synthesis frame center relative to {reference} (ms)", ylabel="Bins per frame")
    graph.legend(loc="upper right")
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="WAV input, mixed down to mono")
    parser.add_argument("--case", choices=("impulse", "burst", "mixed", "tone"), nargs="+")
    parser.add_argument("--impulse-offset", type=int, nargs="+", help="run only impulses at 220 ms plus these sample offsets")
    parser.add_argument("--ratio", type=float, nargs="+", default=[1.5, 2.0])
    parser.add_argument("--window", type=int, default=2048)
    parser.add_argument("--fft", type=int, default=4096)
    parser.add_argument("--hop", type=int, default=256, help="synthesis hop in samples")
    parser.add_argument("--radius", type=int, help="anchor propagation radius in synthesis frames")
    parser.add_argument("--global-heap", action="store_true", help="also compare whole-STFT integration with and without anchors")
    parser.add_argument("--propagation-edges", choices=("centered", "direct"), default="centered",
                        help="edges for local resets and local/global heap propagation")
    parser.add_argument("--tolerance", type=float, default=1e-6)
    parser.add_argument("--anchors", choices=("known", "detected", "both"), default="both")
    parser.add_argument("--onset-ms", type=float, nargs="+", help="reference onsets for a WAV input")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    if args.out is None:
        folder = "runs" if args.propagation_edges == "centered" else "runs-direct"
        args.out = WORKSPACE / "temp" / "rtpghi" / folder
    if sum(value is not None for value in (args.input, args.case, args.impulse_offset)) > 1:
        parser.error("choose one of --input, --case, or --impulse-offset")
    if args.onset_ms and not args.input:
        parser.error("--onset-ms requires --input")
    if (args.window < 32 or args.window % 2 or args.fft < args.window or args.fft % 2
            or not 0 < args.hop <= args.window // 2 or not 0 < args.tolerance < 1):
        parser.error("use an even FFT >= even window >= 32, hop in [1, window/2], and 0 < tolerance < 1")
    radius = args.radius if args.radius is not None else (args.window - 1) // args.hop
    if radius < 0 or any(not np.isfinite(ratio) or ratio <= 0 for ratio in args.ratio):
        parser.error("radius must be nonnegative and ratios must be finite and positive")
    hops = [round(args.hop / ratio) for ratio in args.ratio]
    if any(not 1 <= hop <= args.window // 2 for hop in hops):
        parser.error("requested ratio puts the analysis hop outside [1, window/2]")
    if args.input:
        rate, signal = input_wav(args.input)
        if args.onset_ms and any(not np.isfinite(ms) or not 0 <= round(ms * rate / 1000) < len(signal)
                                for ms in args.onset_ms):
            parser.error("reference onsets must fall inside the input WAV")
        onsets = tuple(sorted({round(ms * rate / 1000) for ms in args.onset_ms})) if args.onset_ms else None
        cases = [Case(args.input.stem, signal, onsets, max(1, round(0.007 * rate)))]
    else:
        rate = 48_000
        if args.impulse_offset:
            try:
                cases = [impulse_case(rate, offset, f"impulse-{offset}") for offset in args.impulse_offset]
            except ValueError as error:
                parser.error(str(error))
        else:
            builtins = builtin_cases(rate)
            cases = [builtins[name] for name in args.case or builtins]
    if args.anchors == "known" and any(case.onsets is None for case in cases):
        parser.error("--anchors known requires --onset-ms for a WAV input")

    args.out.mkdir(parents=True, exist_ok=True)
    summary = []
    for case in cases:
        detected = tuple(detect_onsets(case.signal, rate))
        wavfile.write(args.out / f"{case.name}-input.wav", rate, case.signal.astype(np.float32))
        modes = {"known": case.onsets, "detected": detected}
        for hop in dict.fromkeys(hops):
            config = Stft(args.window, args.fft, hop, args.hop)
            for mode, onsets in modes.items():
                if onsets is None or args.anchors not in ("both", mode):
                    continue
                targets = case.onsets if case.onsets is not None else detected
                name = f"{case.name}-ha{hop}-hs{args.hop}-{mode}"
                result = render(case, rate, config, onsets, radius, args.tolerance, global_heap=args.global_heap,
                                propagation_scheme=args.propagation_edges)
                for method, output in result.outputs.items():
                    wavfile.write(args.out / f"{name}-{method}.wav", rate, output.astype(np.float32))
                if not args.no_plot:
                    plot_results(args.out / f"{name}.png", case, rate, config, result, targets)
                counts = result.anchor_bins
                print(f"{name} ({config.ratio:.6g}) onsets={onsets}, fixed bins/frame={counts[counts > 0].tolist()}",
                      flush=True)
                for metrics in measure(result, case, config, targets):
                    summary.append(dict(case=case.name, ratio=config.ratio, anchors=mode, radius=radius,
                                        supplied_onsets=" ".join(map(str, onsets)), fixed_bins=int(counts.sum()),
                                        propagated_bins=int(result.local_coverage.sum() - counts.sum()), **metrics))
    with (args.out / "metrics.csv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=summary[0].keys())
        writer.writeheader()
        writer.writerows(summary)


if __name__ == "__main__":
    main()
