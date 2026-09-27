# -*- coding: utf-8 -*-
"""使用 FCPE、RMVPE 和 swiftf0 对 wormhole.wav 做 F0 对比。"""
from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
from pathlib import Path
import os
import shutil
import subprocess
import sys
import types
import urllib.request
import librosa
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import soundfile as sf


SCRIPT_DIR = Path(__file__).resolve().parent
LIB_DIR = SCRIPT_DIR.parents[1]
REPO_DIR = LIB_DIR.parents[1]
DEFAULT_INPUT = LIB_DIR / "work_dir" / "input" / "wormhole.wav"
DEFAULT_OUTPUT = SCRIPT_DIR / "output" / "wormhole_f0.png"
RMVPE_MODEL = SCRIPT_DIR / "models" / "rmvpe.pt"
RMVPE_URL = "https://huggingface.co/lj1995/VoiceConversionWebUI/resolve/main/rmvpe.pt"


# ------------------------------------------------------------
# 音频
# ------------------------------------------------------------


def read_audio(path: Path) -> tuple[np.ndarray, int]:
    """读取单声道浮点音频。"""
    audio, sample_rate = sf.read(path, always_2d=True, dtype="float32")
    mono = np.mean(audio, axis=1, dtype=np.float32)
    return mono, int(sample_rate)


def resample_to_16k(audio: np.ndarray, sample_rate: int) -> np.ndarray:
    """将模型输入转换为 16 kHz。"""
    if sample_rate == 16000:
        return np.asarray(audio, dtype=np.float32)
    return np.asarray(
        librosa.resample(audio, orig_sr=sample_rate, target_sr=16000), dtype=np.float32
    )


# ------------------------------------------------------------
# FCPE
# ------------------------------------------------------------


def _install_torchaudio_fallback() -> None:
    """为仅需要 Resample 名称的 FCPE 导入提供最小兼容模块。"""
    import torch
    import torch.nn.functional as functional

    class Resample(torch.nn.Module):
        def __init__(self, orig_freq: int, new_freq: int, **_: object) -> None:
            super().__init__()
            self.orig_freq = orig_freq
            self.new_freq = new_freq

        def forward(self, waveform: "torch.Tensor") -> "torch.Tensor":
            if self.orig_freq == self.new_freq:
                return waveform
            size = round(waveform.shape[-1] * self.new_freq / self.orig_freq)
            return functional.interpolate(
                waveform.unsqueeze(1), size=size, mode="linear", align_corners=False
            ).squeeze(1)

    transforms = types.ModuleType("torchaudio.transforms")
    transforms.Resample = Resample
    torchaudio = types.ModuleType("torchaudio")
    torchaudio.transforms = transforms
    sys.modules["torchaudio"] = torchaudio
    sys.modules["torchaudio.transforms"] = transforms


def load_fcpe_model():
    """加载 FCPE bundled 模型，兼容当前环境中不可加载的 torchaudio。"""
    try:
        from torchfcpe import spawn_bundled_infer_model
    except (ImportError, OSError, RuntimeError):
        for name in list(sys.modules):
            if name == "torchaudio" or name.startswith("torchaudio."):
                del sys.modules[name]
        _install_torchaudio_fallback()
        from torchfcpe import spawn_bundled_infer_model
    return spawn_bundled_infer_model(device="cpu")


def run_fcpe(audio_16k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """运行 FCPE，返回帧时间和 Hz。"""
    import torch

    model = load_fcpe_model()
    waveform = torch.from_numpy(audio_16k).float().unsqueeze(0).unsqueeze(-1)
    with torch.inference_mode():
        f0 = model.infer(
            waveform,
            sr=16000,
            decoder_mode="local_argmax",
            threshold=0.006,
            f0_min=32.70,
            f0_max=1975.5,
            interp_uv=False,
        )
    values = f0[0, :, 0].detach().cpu().numpy().astype(np.float32)
    times = np.arange(values.size, dtype=np.float32) * (160.0 / 16000.0)
    return times, values


# ------------------------------------------------------------
# RMVPE
# ------------------------------------------------------------


def ensure_rmvpe_model() -> Path:
    """下载 RMVPE 权重到本目录的未跟踪 models/ 文件夹。"""
    RMVPE_MODEL.parent.mkdir(parents=True, exist_ok=True)
    if not RMVPE_MODEL.exists():
        print(f"下载 RMVPE 权重: {RMVPE_URL}")
        urllib.request.urlretrieve(RMVPE_URL, RMVPE_MODEL)
    return RMVPE_MODEL


def load_rmvpe_class():
    """只加载 rvc-python 中的 RMVPE 文件，避免导入完整 RVC 推理栈。"""
    try:
        package_root = Path(importlib.metadata.distribution("rvc-python").locate_file("rvc_python"))
    except importlib.metadata.PackageNotFoundError as exc:
        raise RuntimeError(
            "未安装 rvc-python；请按 README 安装 RMVPE 依赖。"
        ) from exc

    # rmvpe.py 仅在 use_jit=True 时使用 jit 的具体函数；这里的离线推理走默认模型。
    package = types.ModuleType("rvc_python")
    package.__path__ = [str(package_root)]
    lib_package = types.ModuleType("rvc_python.lib")
    lib_package.__path__ = [str(package_root / "lib")]
    jit_module = types.ModuleType("rvc_python.lib.jit")
    sys.modules["rvc_python"] = package
    sys.modules["rvc_python.lib"] = lib_package
    sys.modules["rvc_python.lib.jit"] = jit_module

    module_name = "rvc_python.lib.rmvpe_for_qwqdsp"
    source = package_root / "lib" / "rmvpe.py"
    spec = importlib.util.spec_from_file_location(module_name, source)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载 RMVPE 实现: {source}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module.RMVPE


def run_rmvpe(audio_16k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """运行 RMVPE，返回帧时间和 Hz。"""
    rmvpe_class = load_rmvpe_class()
    model = rmvpe_class(str(ensure_rmvpe_model()), is_half=False, device="cpu")
    values = np.asarray(model.infer_from_audio(audio_16k, thred=0.03), dtype=np.float32)
    times = np.arange(values.size, dtype=np.float32) * (160.0 / 16000.0)
    return times, values


# ------------------------------------------------------------
# swiftf0
# ------------------------------------------------------------


def find_clang() -> str:
    """查找仓库使用的 clang++。"""
    candidates = [
        os.environ.get("CXX", ""),
        r"C:\Program Files\LLVM\bin\clang++.exe",
        shutil.which("clang++") or "",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return candidate
    raise RuntimeError("找不到 clang++，swiftf0 helper 无法编译。")


def build_swiftf0_helper(rebuild: bool) -> Path:
    """编译使用仓库内 swiftf0 模型权重的离线 helper。"""
    source = SCRIPT_DIR / "swiftf0_runner.cpp"
    executable = SCRIPT_DIR / "swiftf0_runner.exe"
    if not rebuild and executable.exists() and executable.stat().st_mtime >= source.stat().st_mtime:
        return executable

    command = [
        find_clang(),
        "-std=c++20",
        "-O2",
        f"-I{REPO_DIR / 'qwqdsp' / 'include'}",
        f"-I{REPO_DIR / 'eigen'}",
        f"-I{REPO_DIR / 'eigen' / 'unsupported'}",
        str(source),
        f"-o{executable}",
    ]
    print("编译 swiftf0 helper...")
    subprocess.run(command, cwd=REPO_DIR, check=True)
    return executable


def run_swiftf0(wav_path: Path, rebuild: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """运行仓库内 swiftf0 helper，返回帧时间、Hz 和置信度。"""
    executable = build_swiftf0_helper(rebuild)
    completed = subprocess.run(
        [str(executable), str(wav_path)],
        cwd=REPO_DIR,
        check=True,
        capture_output=True,
        text=True,
    )
    rows = []
    for line in completed.stdout.splitlines():
        if line and not line.startswith("#"):
            rows.append([float(value) for value in line.split(",")])
    if not rows:
        raise RuntimeError("swiftf0 helper 没有输出任何分析帧。")
    data = np.asarray(rows, dtype=np.float32)
    return data[:, 0], data[:, 1], data[:, 2]


# ------------------------------------------------------------
# 绘图
# ------------------------------------------------------------


def plot_results(
    audio: np.ndarray,
    sample_rate: int,
    fcpe: tuple[np.ndarray, np.ndarray],
    rmvpe: tuple[np.ndarray, np.ndarray],
    swiftf0: tuple[np.ndarray, np.ndarray, np.ndarray],
    output_path: Path,
) -> None:
    """绘制波形与三个模型的 F0 四子图。"""
    duration = audio.size / sample_rate
    time = np.arange(audio.size, dtype=np.float32) / sample_rate
    swift_time, swift_values, swift_confidence = swiftf0
    swift_values = np.where(swift_confidence >= 0.5, swift_values, np.nan)

    figure, axes = plt.subplots(
        4,
        1,
        sharex=True,
        figsize=(16, 10),
        gridspec_kw={"height_ratios": [1.2, 1.0, 1.0, 1.0]},
    )
    axes[0].plot(time, audio, color="black", linewidth=0.35)
    axes[0].set_ylabel("waveform")
    axes[0].set_title("wormhole.wav")
    axes[0].grid(True, alpha=0.25)

    plots = [
        (axes[1], fcpe[0], fcpe[1], "FCPE", "tab:blue"),
        (axes[2], rmvpe[0], rmvpe[1], "RMVPE", "tab:orange"),
        (axes[3], swift_time, swift_values, "swiftf0", "tab:green"),
    ]
    for axis, times, values, name, color in plots:
        values = np.where(values > 0.0, values, np.nan)
        axis.plot(times, values, color=color, linewidth=0.9)
        axis.set_ylabel("Hz")
        axis.set_ylim(20.0, 2200.0)
        axis.grid(True, alpha=0.25)
        axis.set_title(name, loc="left", fontweight="bold")

    axes[-1].set_xlabel("Time (s)")
    axes[-1].set_xlim(0.0, duration)
    figure.suptitle("F0 analysis: FCPE / RMVPE / swiftf0", fontsize=14)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.965))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


# ------------------------------------------------------------
# main
# ------------------------------------------------------------


def main() -> None:
    """执行三种 F0 分析并保存四子图。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wav", nargs="?", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--rebuild-swiftf0", action="store_true")
    args = parser.parse_args()

    wav_path = args.wav.resolve()
    if not wav_path.exists():
        raise FileNotFoundError(wav_path)

    audio, sample_rate = read_audio(wav_path)
    audio_16k = resample_to_16k(audio, sample_rate)
    print(f"输入: {wav_path}")
    print(f"采样率: {sample_rate} Hz, 时长: {audio.size / sample_rate:.3f} s")

    print("运行 FCPE...")
    fcpe = run_fcpe(audio_16k)
    print("运行 RMVPE...")
    rmvpe = run_rmvpe(audio_16k)
    print("运行 swiftf0...")
    swiftf0 = run_swiftf0(wav_path, args.rebuild_swiftf0)
    plot_results(audio, sample_rate, fcpe, rmvpe, swiftf0, args.out.resolve())
    print(f"输出: {args.out.resolve()}")


if __name__ == "__main__":
    main()
