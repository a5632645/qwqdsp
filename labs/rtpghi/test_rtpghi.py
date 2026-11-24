import unittest

import numpy as np

from chirp_gradient import (
    Chirp, chirp_gradients, chirp_spectrum_correction, estimate_moments, gradient_edges,
    moment_gradients, phase_potential_edges, sample_grid,
)
from rtpghi import (
    Edges, Stft, analyze, causal, detect_onsets, phase_anchors, phase_edges, phase_loop_residual, propagate,
    propagate_global, wola, wrap,
)


class PhaseIntegrationTest(unittest.TestCase):
    def test_physical_frequency_rescaling_restores_chirp_gradient_closure(self):
        config = Stft(32, 64, 5, 10)
        sigma, beta = 5, 0.8
        rate = beta / sigma ** 2
        omega = (2 * np.pi / config.fft_size) * np.arange(10)
        instantaneous = 0.3 + rate * config.analysis_hop * np.arange(4)[:, None]
        time, frequency = chirp_gradients(omega, instantaneous, beta, sigma)
        copied = gradient_edges(time, frequency, config)
        paper = gradient_edges(time, config.ratio * frequency, config)
        expected = -(1 - 1 / config.ratio) * beta ** 2 / (1 + beta ** 2) * config.synthesis_hop * (omega[1] - omega[0])
        np.testing.assert_allclose(phase_loop_residual(copied), expected, atol=1e-12)
        np.testing.assert_allclose(phase_loop_residual(paper), 0, atol=1e-12)
        wide_time, wide_frequency = chirp_gradients(omega, instantaneous, beta, sigma * np.sqrt(config.ratio))
        np.testing.assert_allclose(wide_time, time)
        np.testing.assert_allclose(wide_frequency, config.ratio * frequency)

    def test_gaussian_moment_fit_matches_sampled_target_chirp(self):
        sigma = 32
        config = Stft(14 * sigma, 512, 48, 96)
        samples = np.arange(-7 * sigma, 7 * sigma) / sigma
        window = np.exp(-0.5 * samples ** 2)
        components = (Chirp(4000, 100_000),)
        rows = np.arange(-4, 5)
        signal, centers, source = sample_grid(components, rows * config.analysis_hop, 1, config, window)
        _, _, target = sample_grid(components, rows * config.synthesis_hop, config.ratio, config, window)
        moment1 = analyze(signal, config.window_size, config.fft_size, centers, weights=window * samples)
        moment2 = analyze(signal, config.window_size, config.fft_size, centers, weights=window * samples ** 2)
        omega = 2 * np.pi * np.fft.rfftfreq(config.fft_size)
        fit = estimate_moments(source, moment1, moment2, omega, sigma)
        retuned = moment_gradients(fit, omega, sigma, config.ratio)
        strong = np.abs(target) > 0.1 * np.abs(target).max()
        self.assertTrue(np.all(fit.confidence[strong] > 0))
        edges = gradient_edges(*retuned, config)
        phase = np.angle(target)
        time_error = wrap(edges.time - np.diff(phase, axis=0))
        frequency_error = wrap(edges.frequency - np.diff(phase, axis=1))
        np.testing.assert_allclose(time_error[strong[:-1] & strong[1:]], 0, atol=1e-8)
        np.testing.assert_allclose(frequency_error[strong[:, :-1] & strong[:, 1:]], 0, atol=1e-8)

        correction, gain = chirp_spectrum_correction(fit, config.ratio)
        base = phase_edges(np.angle(source), config, centered=False)
        edges = phase_potential_edges(base, correction)
        time_error = wrap(edges.time - np.diff(phase, axis=0))
        frequency_error = wrap(edges.frequency - np.diff(phase, axis=1))
        np.testing.assert_allclose(time_error[strong[:-1] & strong[1:]], 0, atol=1e-8)
        np.testing.assert_allclose(frequency_error[strong[:, :-1] & strong[:, 1:]], 0, atol=1e-8)
        np.testing.assert_allclose((np.abs(source) * gain)[strong], np.abs(target)[strong], rtol=1e-8)

    def test_phase_potential_preserves_cell_circulation(self):
        rng = np.random.default_rng(16)
        base = Edges(rng.uniform(-1, 1, (3, 6)), rng.uniform(-1, 1, (4, 5)))
        potential = rng.normal(size=(4, 6)) * (np.arange(6) % 2)
        updated = phase_potential_edges(base, potential)
        np.testing.assert_allclose(wrap(phase_loop_residual(updated) - phase_loop_residual(base)), 0, atol=1e-12)

    def test_centered_difference_unwraps_each_step_before_averaging(self):
        config = Stft(32, 32, 8, 12)
        carrier = 2 * np.pi * config.analysis_hop / config.fft_size
        time_residual = np.array([0, 2.8, 2.8, 0.5])
        frequency_phase = np.array([0, 2.8, 5.6, 6.1])
        phase = wrap(time_residual[:, None] + np.arange(4)[:, None] * carrier * np.arange(4) + frequency_phase)
        edges = phase_edges(phase, config, centered=True)
        np.testing.assert_allclose(edges.time[:, 0], config.ratio * np.array([2.1, 0.125, -1.725]), atol=1e-12)
        np.testing.assert_allclose(edges.frequency[0], config.ratio * np.array([2.8, 2.225, 1.075]))

    def test_direct_edge_wrap_turn_closes_only_at_integer_scale(self):
        phase = np.array([[0, np.pi / 2], [-np.pi / 2, -np.pi]])
        for synthesis_hop, expected in ((12, np.pi), (16, 0)):
            with self.subTest(synthesis_hop=synthesis_hop):
                edges = phase_edges(phase, Stft(32, 32, 8, synthesis_hop), centered=False)
                np.testing.assert_allclose(np.abs(phase_loop_residual(edges)), [[expected]], atol=1e-12)

    def test_integer_direct_edges_integrate_from_one_seed_independently_of_heap_priority(self):
        config = Stft(32, 32, 8, 16)
        rng = np.random.default_rng(3)
        phase = rng.uniform(-np.pi, np.pi, (7, config.fft_size // 2 + 1))
        magnitude = rng.uniform(0.1, 0.9, phase.shape)
        magnitude[3, 2] = 1.0
        edges = phase_edges(phase, config, centered=False)
        anchors = np.full(phase.shape, np.nan)
        expected = config.ratio * (phase - phase[3, 2]) + phase[3, 2]
        for bottleneck in (False, True):
            with self.subTest(bottleneck=bottleneck):
                integrated = propagate_global(magnitude, phase, edges, anchors, 1e-6, bottleneck=bottleneck)
                self.assertEqual(integrated.seed_count, 1)
                np.testing.assert_allclose(wrap(integrated.phase - expected), 0, atol=1e-12)

    def test_reverse_propagation_subtracts_time_edges(self):
        expected = np.arange(5)[:, None] * 0.17 + np.arange(4) * 0.24
        magnitude = np.ones_like(expected)
        anchors = np.full_like(expected, np.nan)
        anchors[2, 1] = expected[2, 1]
        edges = Edges(np.diff(expected, axis=0), np.diff(expected, axis=1))
        forward = propagate(magnitude, edges, anchors, 1, 1e-6, backward=False)
        backward = propagate(magnitude, edges, anchors, 1, 1e-6, backward=True)
        np.testing.assert_allclose(forward[2:4], expected[2:4])
        np.testing.assert_allclose(backward[1:4], expected[1:4])
        self.assertTrue(np.isnan(forward[:2]).all())
        self.assertTrue(np.isnan(backward[[0, 4]]).all())
        self.assertEqual(backward[2, 1], anchors[2, 1])

    def test_strong_path_wins_without_overwriting_competing_anchors(self):
        magnitude = np.array([[8, 0.25, 4, 4, 5]])
        edges = Edges(np.zeros((0, 5)), np.zeros((1, 4)))
        anchors = np.array([[0, np.nan, np.nan, np.nan, 1]])
        fixed = propagate(magnitude, edges, anchors, 0, 1e-6, backward=True)
        np.testing.assert_allclose(fixed, [[0, 0, 1, 1, 1]])

    def test_global_heap_seeds_each_disconnected_component_at_its_peak(self):
        magnitude = np.array([[1, 1, 0, 2], [3, 8, 0, 3], [5, 1, 0, 7]], dtype=float)
        expected = 0.17 * np.arange(3)[:, None] + 0.24 * np.arange(4)
        edges = Edges(np.diff(expected, axis=0), np.diff(expected, axis=1))
        analysis = np.full_like(expected, -1)
        analysis[1, 1], analysis[2, 3] = expected[1, 1], expected[2, 3]
        anchors = np.full_like(expected, np.nan)
        for bottleneck in (False, True):
            with self.subTest(bottleneck=bottleneck):
                integrated = propagate_global(magnitude, analysis, edges, anchors, 1e-6, bottleneck=bottleneck)
                np.testing.assert_allclose(integrated.phase[magnitude > 0], expected[magnitude > 0])
                self.assertTrue(np.isnan(integrated.phase[:, 2]).all())
                self.assertEqual(integrated.seed_count, 2)
                np.testing.assert_array_equal(integrated.anchored_by_frame, [0, 0, 0])

        anchors[2, 0] = expected[2, 0] + 0.5
        for bottleneck in (False, True):
            with self.subTest(anchored=True, bottleneck=bottleneck):
                integrated = propagate_global(magnitude, analysis, edges, anchors, 1e-6, bottleneck=bottleneck)
                np.testing.assert_allclose(integrated.phase[:, :2], expected[:, :2] + 0.5)
                np.testing.assert_allclose(integrated.phase[:, 3], expected[:, 3])
                self.assertEqual(integrated.phase[2, 0], anchors[2, 0])
                self.assertEqual(integrated.seed_count, 1)
                np.testing.assert_array_equal(integrated.anchored_by_frame, [2, 2, 2])

    def test_wola_identity(self):
        config = Stft(256, 512, 48, 48)
        signal = np.random.default_rng(1).normal(size=1537)
        centers = np.arange(34) * config.analysis_hop
        spectrum = analyze(signal, config.window_size, config.fft_size, centers)
        np.testing.assert_allclose(wola(spectrum, config, len(signal)), signal, atol=1e-12)

    def test_off_grid_impulse_anchor_and_preringing(self):
        config = Stft(256, 512, 31, 48)
        signal = np.zeros(1024)
        onset = 521
        signal[onset] = 1
        self.assertEqual(detect_onsets(signal, 48_000), [onset])
        centers = np.arange(35) * config.analysis_hop
        spectrum = analyze(signal, config.window_size, config.fft_size, centers)
        phase, magnitude = np.angle(spectrum), np.abs(spectrum)
        edges = phase_edges(phase, config, centered=True)
        anchors = phase_anchors(signal, spectrum, config, [onset], 48_000, 1e-6)
        row = round(onset / config.analysis_hop)
        target = round(onset * config.ratio)
        offset = target - row * config.synthesis_hop
        omega = 2 * np.pi / config.fft_size * np.arange(spectrum.shape[1])
        np.testing.assert_allclose(wrap(anchors[row] + omega * offset), 0, atol=1e-12)

        def render(backward):
            fixed = propagate(magnitude, edges, anchors, 5, 1e-6, backward=backward)
            estimate, _ = causal(magnitude, phase, edges, fixed, 1e-6)
            return wola(magnitude * np.exp(1j * estimate), config, round(len(signal) * config.ratio))

        forward, bidirectional = render(False), render(True)
        self.assertEqual(np.argmax(np.abs(bidirectional)), target)
        self.assertLess(np.linalg.norm(bidirectional[:target]), 0.1 * np.linalg.norm(forward[:target]))


if __name__ == "__main__":
    unittest.main()
