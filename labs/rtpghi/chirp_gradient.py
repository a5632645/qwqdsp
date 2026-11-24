"""Compare time-stretch phase gradients against sampled Gaussian-window chirps."""

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import wavfile

from rtpghi import Edges, Stft, analyze, phase_edges, phase_loop_residual, propagate_global, wola, wrap


SAMPLE_RATE = 48_000
TOLERANCE = 1e-6
WORKSPACE = Path(__file__).resolve().parents[3]


@dataclass(frozen=True)
class Chirp:
    frequency_hz: float
    rate_hz_per_second: float
    amplitude: float = 1.0
    phase: float = 0.0


@dataclass(frozen=True)
class MomentFit:
    time: np.ndarray
    frequency: np.ndarray
    beta: np.ndarray
    detuning_sigma: np.ndarray
    confidence: np.ndarray


CASES = {
    "single": (Chirp(4000, 20_000),),
    "crossing": (Chirp(4000, 20_000, 0.55), Chirp(4000, -16_000, 0.40, 1.1)),
}


def chirp_signal(components, samples, ratio):
    seconds = samples / SAMPLE_RATE
    return sum(chirp.amplitude * np.cos(2 * np.pi * chirp.frequency_hz * seconds
                                       + np.pi * chirp.rate_hz_per_second / ratio * seconds ** 2 + chirp.phase)
               for chirp in components)


def chirp_gradients(omega, instantaneous_frequency, beta, sigma):
    """Physical phase derivatives in radians/sample and samples."""
    detuning = omega - instantaneous_frequency
    denominator = 1 + beta ** 2
    time = instantaneous_frequency + beta ** 2 / denominator * detuning
    frequency = -sigma ** 2 * beta / denominator * detuning
    return time, frequency


def gradient_edges(time, frequency, config):
    bin_width = 2 * np.pi / config.fft_size
    return Edges(0.5 * config.synthesis_hop * (time[:-1] + time[1:]),
                 0.5 * bin_width * (frequency[:, :-1] + frequency[:, 1:]))


def estimate_moments(c0, c1, c2, omega, sigma):
    """Fit an isolated Gaussian-windowed linear chirp with three STFTs."""
    magnitude = np.abs(c0)
    active = magnitude > TOLERANCE * magnitude.max()
    mu = np.divide(c1, c0, out=np.zeros_like(c0), where=active)
    nu = np.divide(c2, c0, out=np.zeros_like(c0), where=active)
    # c1 and c2 use (u/sigma) and (u/sigma)^2 times the Gaussian.
    # h is the normalized second frequency derivative of log(C).
    h = mu ** 2 - nu
    inverse = np.divide(1, h, out=np.zeros_like(h), where=np.abs(h) > TOLERANCE)
    detuning_sigma = -1j * mu * inverse
    time = omega + mu.imag / sigma
    frequency = -sigma * mu.real
    fit_error = np.maximum(np.abs(inverse.real + 1), np.abs(detuning_sigma.imag))
    confidence = active * np.clip((0.15 - fit_error) / 0.10, 0, 1)
    return MomentFit(time, frequency, inverse.imag, detuning_sigma.real, confidence)


def moment_gradients(fit, omega, sigma, ratio):
    instantaneous_frequency = omega - fit.detuning_sigma / sigma
    target_time, target_frequency = chirp_gradients(omega, instantaneous_frequency, fit.beta / ratio, sigma)
    trusted = fit.confidence > 0
    return (np.where(trusted, target_time, fit.time), np.where(trusted, target_frequency, ratio * fit.frequency))


def phase_potential_edges(edges, correction):
    return Edges(edges.time + np.diff(correction, axis=0), edges.frequency + np.diff(correction, axis=1))


def chirp_spectrum_correction(fit, ratio):
    """Return the phase offset and magnitude gain for a Gaussian-windowed chirp."""
    beta_a = fit.beta
    beta_s = beta_a / ratio
    denominator_a = 1 + beta_a ** 2
    denominator_s = 1 + beta_s ** 2
    detuning_power = fit.detuning_sigma ** 2
    phase_a = 0.5 * (np.arctan(beta_a) - beta_a / denominator_a * detuning_power)
    phase_s = 0.5 * (np.arctan(beta_s) - beta_s / denominator_s * detuning_power)
    correction = fit.confidence * (phase_s - ratio * phase_a)
    log_gain = 0.25 * (np.log(denominator_a) - np.log(denominator_s))
    log_gain -= 0.5 * detuning_power * (1 / denominator_s - 1 / denominator_a)
    # Limit expansion of poorly supported spectral tails during compression.
    gain = np.exp(np.minimum(fit.confidence * log_gain, np.log(4)))
    return correction, gain


def sample_grid(components, centers, ratio, config, weights):
    half = config.window_size // 2
    samples = np.arange(centers[0] - half, centers[-1] + half + 1)
    signal = chirp_signal(components, samples, ratio)
    spectrum = analyze(signal, config.window_size, config.fft_size, centers - samples[0], weights=weights)
    return signal, centers - samples[0], spectrum


def measure_edges(edges, reference, magnitude, common):
    time = common[:-1] & common[1:]
    frequency = common[:, :-1] & common[:, 1:]
    cells = time[:, :-1] & time[:, 1:]
    time_error = wrap((edges.time - reference.time)[time])
    freq_error = wrap((edges.frequency - reference.frequency)[frequency])
    time_magnitude = np.minimum(magnitude[:-1], magnitude[1:])
    freq_magnitude = np.minimum(magnitude[:, :-1], magnitude[:, 1:])
    cell_magnitude = np.minimum(time_magnitude[:, :-1], time_magnitude[:, 1:])
    loop_error = phase_loop_residual(edges)[cells]
    return dict(time_edge_rms_rad=float(np.sqrt(np.average(time_error ** 2, weights=time_magnitude[time] ** 2))),
                frequency_edge_rms_rad=float(np.sqrt(np.average(freq_error ** 2, weights=freq_magnitude[frequency] ** 2))),
                loop_rms_rad=float(np.sqrt(np.average(loop_error ** 2, weights=cell_magnitude[cells] ** 2))))


def snr(signal, reference):
    energy = np.sum(reference ** 2)
    error = np.sum((signal - reference) ** 2)
    return float(-10 * np.log10(max(error / energy, 1e-30)))


def run_case(name, sigma_ms, config, weights, samples, out):
    components = CASES[name]
    sigma = sigma_ms * SAMPLE_RATE / 1000
    half_frames = round(0.11 * SAMPLE_RATE / config.analysis_hop)
    frames = np.arange(-half_frames, half_frames + 1)
    ta, ts = frames * config.analysis_hop, frames * config.synthesis_hop
    source, centers, analysis = sample_grid(components, ta, 1, config, weights)
    target_signal, _, target = sample_grid(components, ts, config.ratio, config, weights)
    moment1 = analyze(source, config.window_size, config.fft_size, centers, weights=weights * samples)
    moment2 = analyze(source, config.window_size, config.fft_size, centers, weights=weights * samples ** 2)
    omega = 2 * np.pi * np.fft.rfftfreq(config.fft_size)
    fit = estimate_moments(analysis, moment1, moment2, omega, sigma)
    retuned = moment_gradients(fit, omega, sigma, config.ratio)
    correction, gain = chirp_spectrum_correction(fit, config.ratio)
    reference_phase = np.angle(target)
    reference = Edges(wrap(np.diff(reference_phase, axis=0)), wrap(np.diff(reference_phase, axis=1)))

    methods = {
        "direct": phase_edges(np.angle(analysis), config, centered=False),
        "centered": phase_edges(np.angle(analysis), config, centered=True),
        "copied-physical": gradient_edges(fit.time, fit.frequency, config),
        "paper-moments": gradient_edges(fit.time, config.ratio * fit.frequency, config),
        "retuned-moments": gradient_edges(*retuned, config),
    }
    methods["retuned-potential"] = phase_potential_edges(methods["paper-moments"], correction)
    methods["retuned-direct-potential"] = phase_potential_edges(methods["direct"], correction)
    if name == "single":
        chirp, = components
        rate = 2 * np.pi * chirp.rate_hz_per_second / SAMPLE_RATE ** 2
        instantaneous_frequency = (2 * np.pi * chirp.frequency_hz / SAMPLE_RATE + rate * ta)[:, None]
        beta = rate * sigma ** 2
        gradients = chirp_gradients(omega, instantaneous_frequency, beta, sigma)
        methods["paper-analytic"] = gradient_edges(gradients[0], config.ratio * gradients[1], config)
        methods["target-analytic"] = gradient_edges(
            *chirp_gradients(omega, instantaneous_frequency, beta / config.ratio, sigma), config,
        )
    methods["target-edges"] = reference

    magnitude = np.abs(analysis)
    target_magnitude = np.abs(target)
    common = ((magnitude > TOLERANCE * magnitude.max())
              & (target_magnitude > TOLERANCE * target_magnitude.max()))
    anchors = np.full(analysis.shape, np.nan)
    peak_bin = np.argmax(magnitude[half_frames])
    anchors[half_frames, peak_bin] = reference_phase[half_frames, peak_bin]
    length = int(ts[-1] - ts[0]) + 1
    half_window = config.window_size // 2
    target_signal = target_signal[half_window:half_window + length]
    trim = min(round(3 * sigma), length // 4)
    interior = slice(trim, length - trim)
    wave_base = f"{name}-sigma{sigma_ms:g}-ha{config.analysis_hop}-hs{config.synthesis_hop}"
    wavfile.write(out / f"{wave_base}-reference.wav", SAMPLE_RATE, target_signal.astype(np.float32))
    exact = wola(target, config, length, weights=weights)
    oracle_snr = snr(exact[interior], target_signal[interior])
    fitted_energy = float(np.sum(magnitude[fit.confidence > 0] ** 2) / np.sum(magnitude ** 2))
    random_phase = np.random.default_rng(0).uniform(-np.pi, np.pi, analysis.shape)

    records = []
    for method, edges in methods.items():
        # Give every method the same true seed per disconnected component.
        integration = propagate_global(magnitude, reference_phase, edges, anchors, TOLERANCE, bottleneck=False)
        phase = np.where(np.isfinite(integration.phase), integration.phase, random_phase)
        signed_error = wrap(phase[common] - reference_phase[common])
        phase_rms = np.sqrt(np.average(signed_error ** 2, weights=target_magnitude[common] ** 2))
        spectrum = magnitude * np.exp(1j * phase)
        measurements = dict(case=name, ratio=config.ratio, sigma_ms=sigma_ms, phase_rms_rad=float(phase_rms),
                            target_magnitude_and_phase_snr_db=oracle_snr, fitted_analysis_energy=fitted_energy,
                            **measure_edges(edges, reference, target_magnitude, common))
        variants = {method: spectrum}
        if method in ("retuned-potential", "retuned-direct-potential"):
            variants[f"{method}-magnitude"] = spectrum * gain
        for variant, coefficients in variants.items():
            output = wola(coefficients, config, length, weights=weights)
            records.append(dict(method=variant, waveform_snr_db=snr(output[interior], target_signal[interior]),
                                **measurements))
            if method in ("direct", "centered", "retuned-moments", "retuned-potential", "retuned-direct-potential",
                          "target-edges"):
                wavfile.write(out / f"{wave_base}-{variant}.wav", SAMPLE_RATE, output.astype(np.float32))
    return records


def plot_rows(rows, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = list(dict.fromkeys((row["case"], row["ratio"]) for row in rows))
    figure, axes = plt.subplots(len(panels), 2, figsize=(11, 3 * len(panels)), constrained_layout=True, squeeze=False)
    series = (
        ("direct", "Direct edges"),
        ("centered", "Centered edges"),
        ("retuned-moments", "Retuned gradients"),
        ("retuned-direct-potential", "Phase correction"),
        ("retuned-direct-potential-magnitude", "Phase and magnitude correction"),
        ("target-edges", "Reference phase, input magnitude"),
    )
    for (name, ratio), (phase_axis, wave_axis) in zip(panels, axes):
        for index, (method, label) in enumerate(series):
            selected = [row for row in rows if (row["case"], row["ratio"], row["method"]) == (name, ratio, method)]
            sigma = [row["sigma_ms"] for row in selected]
            color = f"C{index}"
            if not method.endswith("-magnitude"):
                phase_axis.semilogy(sigma, [max(row["phase_rms_rad"], 1e-12) for row in selected], "o-", color=color)
            wave_axis.plot(sigma, [row["waveform_snr_db"] for row in selected], "o-", label=label, color=color)
        phase_axis.set(title=f"{name} · {ratio:.4g}×", ylabel="Phase RMS error (rad)", xlabel="Gaussian σ (ms)")
        wave_axis.set(ylabel="Waveform SNR (dB)", xlabel="Gaussian σ (ms)")
    axes[0, 1].legend(fontsize=8)
    figure.savefig(out / "chirp-gradient.png", dpi=140)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=tuple(CASES), nargs="+", default=list(CASES))
    parser.add_argument("--sigma-ms", type=float, nargs="+", default=[1, 4, 8])
    parser.add_argument("--ratio", type=float, nargs="+", default=[2])
    parser.add_argument("--out", type=Path, default=WORKSPACE / "temp" / "rtpghi" / "chirp-gradient")
    args = parser.parse_args()
    if any(not 0.5 <= width <= 12 for width in args.sigma_ms):
        parser.error("Gaussian sigma must be between 0.5 and 12 ms")
    if any(not 0.5 <= ratio <= 4 for ratio in args.ratio):
        parser.error("ratio must be between 0.5 and 4")
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    synthesis_hop = 128
    for ratio in args.ratio:
        analysis_hop = round(synthesis_hop / ratio)
        for sigma_ms in sorted(set(args.sigma_ms)):
            sigma = sigma_ms * SAMPLE_RATE / 1000
            half = int(np.ceil(7 * sigma))
            samples = np.arange(-half, half) / sigma
            weights = np.exp(-0.5 * samples ** 2)
            fft_size = 1 << int(2 * half - 1).bit_length()
            config = Stft(2 * half, fft_size, analysis_hop, synthesis_hop)
            for name in args.case:
                rows.extend(run_case(name, sigma_ms, config, weights, samples, args.out))
                print(f"{name}: {config.ratio:.6g}x, sigma={sigma_ms:g} ms", flush=True)
    with (args.out / "metrics.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    plot_rows(rows, args.out)


if __name__ == "__main__":
    main()
