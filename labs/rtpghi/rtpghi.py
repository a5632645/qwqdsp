"""Offline RTPGHI with centered STFT frames and optional fixed phase bins."""

from dataclasses import dataclass
from heapq import heapify, heappop, heappush

import numpy as np
from scipy.signal import find_peaks


def wrap(phase):
    return (phase + np.pi) % (2 * np.pi) - np.pi


@dataclass(frozen=True)
class Stft:
    window_size: int
    fft_size: int
    analysis_hop: int
    synthesis_hop: int

    @property
    def ratio(self):
        return self.synthesis_hop / self.analysis_hop


@dataclass(frozen=True)
class Edges:
    # Radians per synthesis step, oriented toward increasing time/frequency.
    time: np.ndarray
    frequency: np.ndarray


@dataclass(frozen=True)
class GridIntegration:
    phase: np.ndarray
    anchored_by_frame: np.ndarray
    seed_count: int


def window(size):
    return np.hanning(size + 1)[:-1]


def analyze(signal, window_size, fft_size, centers, *, weights=None):
    half = window_size // 2
    starts = centers - half
    left = max(0, -int(starts[0]))
    right = max(0, int(starts[-1]) + window_size - len(signal))
    padded = np.pad(signal, (left, right))
    if weights is None:
        weights = window(window_size)
    frames = padded[starts[:, None] + left + np.arange(window_size)] * weights
    frames = np.pad(frames, ((0, 0), (0, fft_size - window_size)))
    return np.fft.rfft(np.roll(frames, -half, axis=1), axis=1)


def wola(spectrum, config, length, *, weights=None):
    half = config.window_size // 2
    if weights is None:
        weights = window(config.window_size)
    frames = np.fft.irfft(spectrum, n=config.fft_size, axis=1)
    frames = np.roll(frames, half, axis=1)[:, :config.window_size] * weights
    squared = weights ** 2
    size = (len(frames) - 1) * config.synthesis_hop + config.window_size
    signal = np.zeros(size)
    norm = np.zeros(size)
    for row, frame in enumerate(frames):
        start = row * config.synthesis_hop
        segment = slice(start, start + config.window_size)
        signal[segment] += frame
        norm[segment] += squared
    signal = signal[half:half + length]
    norm = norm[half:half + length]
    return np.divide(signal, norm, out=np.zeros_like(signal), where=norm > 0)


def phase_edges(phase, config, *, centered):
    carrier = 2 * np.pi * config.analysis_hop / config.fft_size * np.arange(phase.shape[1])
    time = wrap(np.diff(phase, axis=0) - carrier) + carrier
    frequency = wrap(np.diff(phase, axis=1))
    if centered:
        # Average the unwrapped steps on either side of each node. Use the
        # available one-sided difference at the outer boundaries.
        dt = np.concatenate((time[:1], 0.5 * (time[:-1] + time[1:]), time[-1:]))
        df = np.concatenate((frequency[:, :1], 0.5 * (frequency[:, :-1] + frequency[:, 1:]),
                             frequency[:, -1:]), axis=1)
        time = 0.5 * (dt[:-1] + dt[1:])
        frequency = 0.5 * (df[:, :-1] + df[:, 1:])
    return Edges(config.ratio * time, config.ratio * frequency)


def phase_loop_residual(edges):
    """Wrapped circulation around each elementary time-frequency cell."""
    return wrap(edges.time[:, :-1] + edges.frequency[1:] - edges.time[:, 1:] - edges.frequency[:-1])


def causal(magnitude, analysis_phase, edges, fixed, tolerance):
    """Algorithm 1, with fixed bins entered alongside previous-frame sources."""
    phase = np.random.default_rng(0).uniform(-np.pi, np.pi, magnitude.shape)
    phase[0] = np.where(np.isfinite(fixed[0]), fixed[0], analysis_phase[0])
    weak_seeds = np.zeros(len(phase), dtype=int)
    peaks = magnitude.max(axis=1)
    bins = phase.shape[1]
    for row in range(1, len(phase)):
        previous = magnitude[row - 1]
        current = magnitude[row]
        known = np.isfinite(fixed[row])
        phase[row, known] = fixed[row, known]
        pending = (current > tolerance * max(peaks[row - 1], peaks[row])) & ~known
        remaining = np.count_nonzero(pending)
        heap = [(-previous[k], 0, int(k)) for k in np.flatnonzero(pending)]
        heap.extend((-current[k], 1, int(k)) for k in np.flatnonzero(known))
        heapify(heap)
        while remaining:
            _, is_current, k = heappop(heap)
            if is_current:
                neighbors = []
                if k + 1 < bins:
                    neighbors.append((k + 1, phase[row, k] + edges.frequency[row, k]))
                if k > 0:
                    neighbors.append((k - 1, phase[row, k] - edges.frequency[row, k - 1]))
            else:
                neighbors = [(k, phase[row - 1, k] + edges.time[row - 1, k])]
            for target, estimate in neighbors:
                if pending[target]:
                    pending[target] = False
                    remaining -= 1
                    phase[row, target] = wrap(estimate)
                    heappush(heap, (-current[target], 1, target))
                    if not is_current and previous[target] < 0.1 * current[target]:
                        weak_seeds[row] += 1
    return phase, weak_seeds


def propagate(magnitude, edges, anchors, radius, tolerance, *, backward):
    """Integrate from immutable anchors within the local frame ranges."""
    locked = np.isfinite(anchors)
    active_rows = anchor_frame_range(locked.any(axis=1), radius)
    pending = (magnitude > tolerance * magnitude.max()) & active_rows[:, None] & ~locked
    return _integrate_grid(magnitude, edges, anchors, pending, locked, backward=backward).phase


def anchor_frame_range(anchor_rows, radius):
    active_rows = np.zeros_like(anchor_rows)
    for row in np.flatnonzero(anchor_rows):
        active_rows[max(0, row - radius):row + radius + 1] = True
    return active_rows


def propagate_global(magnitude, analysis_phase, edges, anchors, tolerance, *, bottleneck=True):
    """Integrate the full grid, using a peak seed for each unanchored component."""
    locked = np.isfinite(anchors)
    pending = (magnitude > tolerance * magnitude.max()) & ~locked
    return _integrate_grid(magnitude, edges, anchors, pending, locked, backward=True,
                           seed_phase=analysis_phase, bottleneck=bottleneck)


def _integrate_grid(magnitude, edges, anchors, pending, locked, *, backward, seed_phase=None, bottleneck=True):
    rows, bins = magnitude.shape
    fixed = anchors.copy()
    strength = np.zeros(magnitude.shape)
    distance = np.zeros(magnitude.shape, dtype=int)
    heap = [(-magnitude[row, k], 0, int(row), int(k)) for row, k in zip(*np.nonzero(locked))]
    heapify(heap)

    def drain_heap():
        while heap:
            priority, steps, row, k = heappop(heap)
            if bottleneck and steps and (-priority != strength[row, k] or steps != distance[row, k]):
                continue
            pending[row, k] = False
            neighbors = []
            if k > 0:
                neighbors.append((row, k - 1, -edges.frequency[row, k - 1]))
            if k + 1 < bins:
                neighbors.append((row, k + 1, edges.frequency[row, k]))
            if row + 1 < rows:
                neighbors.append((row + 1, k, edges.time[row, k]))
            if backward and row > 0:
                neighbors.append((row - 1, k, -edges.time[row - 1, k]))
            for next_row, next_k, delta in neighbors:
                target = (next_row, next_k)
                if not pending[target]:
                    continue
                capacity = magnitude[target]
                if bottleneck:
                    capacity = min(-priority, capacity)
                    if capacity < strength[target] or (capacity == strength[target] and steps + 1 >= distance[target]):
                        continue
                    strength[target] = capacity
                    distance[target] = steps + 1
                else:
                    # PGHI assigns a phase on first discovery and queues by bin magnitude.
                    pending[target] = False
                fixed[target] = wrap(fixed[row, k] + delta)
                heappush(heap, (-capacity, steps + 1, *target))

    drain_heap()
    anchored_by_frame = np.count_nonzero(np.isfinite(fixed), axis=1)
    seed_count = 0
    if seed_phase is not None:
        # Pick the strongest still unvisited bin each time a component ends.
        candidates = np.flatnonzero(pending)
        order = np.argsort(-magnitude.flat[candidates], kind="stable")
        for candidate in candidates[order]:
            row, k = divmod(int(candidate), bins)
            if pending[row, k]:
                fixed[row, k] = seed_phase[row, k]
                heap.append((-magnitude[row, k], 0, row, k))
                seed_count += 1
                drain_heap()
    return GridIntegration(fixed, anchored_by_frame, seed_count)


def detect_onsets(signal, sample_rate):
    """Broadband positive spectral flux, refined on the input waveform."""
    size = max(32, 2 * round(0.0027 * sample_rate))
    hop = size // 4
    centers = np.arange(-hop, len(signal) + size // 2, hop)
    magnitude = np.abs(analyze(signal, size, size, centers))[:, 1:]
    previous = np.pad(magnitude[:-1], ((0, 0), (1, 1)))
    reference = np.maximum(np.maximum(previous[:, :-2], previous[:, 1:-1]), previous[:, 2:])
    rise = np.maximum(magnitude[1:] - reference, 0)
    significant = rise > 0.001 * magnitude.max()
    flux = np.sum(rise * significant, axis=1)
    flux[np.mean(significant, axis=1) < 0.12] = 0
    if not np.any(flux):
        return []
    peaks, _ = find_peaks(np.pad(flux, 1), height=0.08 * flux.max(),
                          distance=max(1, round(0.012 * sample_rate / hop)))
    attack = np.abs(np.diff(signal, prepend=0))
    onsets = []
    for peak in peaks:
        center = int(centers[peak])
        start = max(0, center - size // 2)
        stop = min(len(signal), center + size // 2 + 1)
        local = attack[start:stop]
        if local.size and np.any(local):
            onsets.append(start + int(np.flatnonzero(local >= 0.7 * local.max())[0]))
    return sorted(set(onsets))


def phase_anchors(signal, spectrum, config, onsets, sample_rate, tolerance):
    """Choose transient-dominated bins and translate their phase to output time."""
    anchors = np.full(spectrum.shape, np.nan)
    confidence = np.zeros(spectrum.shape)
    size = min(config.window_size, max(32, 2 * round(0.0027 * sample_rate)))
    omega = 2 * np.pi / config.fft_size * np.arange(spectrum.shape[1])
    for onset in onsets:
        row = round(onset / config.analysis_hop)
        output_onset = round(config.ratio * onset)
        displacement = output_onset - row * config.synthesis_hop - (onset - row * config.analysis_hop)
        short = np.abs(analyze(signal, size, config.fft_size, np.array([onset - size, onset])))
        before, event = short
        current = np.abs(spectrum[row])
        dominant = (event > 2 * before) & (event >= 0.5 * current) & (current > tolerance * current.max())
        strength = np.where(dominant, event, 0)
        mask = (strength > 0.15 * strength.max()) & (strength > confidence[row])
        anchors[row, mask] = wrap(np.angle(spectrum[row, mask]) - omega[mask] * displacement)
        confidence[row, mask] = strength[mask]
    return anchors
