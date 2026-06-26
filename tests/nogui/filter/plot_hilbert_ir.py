"""Plot magnitude response of the Hilbert filter impulse response from WAV."""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.wavfile as wav
from pathlib import Path

wav_path = Path(__file__).parent / "../../work_dir/output/cheby_hilbert_impulse.wav"

sr, data = wav.read(wav_path)
print(f"Sample rate: {sr} Hz,  shape: {data.shape}")

if data.dtype == np.int16:
    data = data.astype(np.float32) / 32768.0
elif data.dtype == np.int32:
    data = data.astype(np.float32) / 2147483648.0

# complex impulse response: ch0=real, ch1=imag
h = data[:, 0] + 1j * data[:, 1]
N = len(h)
print(f"Impulse response length: {N}")

# FFT
H = np.fft.fft(h)
H = np.fft.fftshift(H)
freq = np.fft.fftshift(np.fft.fftfreq(N, d=1.0 / sr))
mag_db = 20 * np.log10(np.abs(H) + 1e-15)

# plot
fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(freq, mag_db, lw=1)
ax.set_title("Magnitude Response (from WAV impulse response)")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Magnitude (dB)")
ax.grid(alpha=0.3)
ax.set_ylim(-180, 5)
ax.set_xlim(-sr / 2, sr / 2)
plt.tight_layout()
plt.show()
