# wormhole.wav F0 三模型对比

使用 FCPE、RMVPE 与仓库自带的 swiftf0 模型，对 `tests/work_dir/input/wormhole.wav` 执行离线基频分析，输出从上到下依次为波形、FCPE、RMVPE、swiftf0。

## 运行

在仓库根目录运行：

```powershell
python qwqdsp/tests/hybrid/f0_compare/run.py
```

默认输出 `qwqdsp/tests/hybrid/f0_compare/output/wormhole_f0.png`。也可以指定 WAV 和输出路径：

```powershell
python qwqdsp/tests/hybrid/f0_compare/run.py path/to/input.wav --out path/to/result.png
```

首次运行会自动编译 swiftf0 helper，并将 RMVPE 官方 `rmvpe.pt` 权重下载至本目录 `models/`。helper 直接复用仓库 `swift_f0_rt` 的 16 kHz 前端、抗混叠抽取器及 SwiftF0 推理模型；分析结果时间戳已补偿抽取器群延迟。swiftf0 曲线隐藏置信度低于 0.5 的帧。

## 依赖

Python 包：`numpy`、`soundfile`、`matplotlib`、`librosa`、`torch`、`torchfcpe`、`rvc-python`。C++ helper 使用仓库配置的 `clang++` 和已有 Eigen 目录。FCPE 与 RMVPE 默认在 CPU 上运行，不要求 GPU。
