# PGHI 离线拉伸示例

从工作区根目录构建。运行时所在的目录就是输出目录。

```powershell
cmake --build build --config Release --target qwqdsp-example-pghi_stretch
.\build\qwqdsp\example\offline\Release\qwqdsp-example-pghi_stretch.exe
```

先输入音频路径，支持 WAV、FLAC 和 MP3，可粘贴带引号的路径。再设置分析参数，直接回车使用默认值。

| 选项 | 默认值 | 说明 |
| --- | --- | --- |
| 分析窗长 | `2048` | `256–8192` 的 2 的幂 |
| 瞬态锚点 | `y` | 用短窗宽带起音定位，固定瞬态频点的合成相位 |
| Gaussian chirp 修正 | `n` | 单声道选项，在单条主导脊线处联合修正相位和幅度 |

分析完成后，可在同一会话中反复渲染。每次输入时长倍率 `0.25–4` 和移调半音数 `−12–+12`。首次默认 `2`、`0`，后续回车沿用上一次的值。时长输入 `r` 可直接重采样，此时时长随音高变化。输入 `q` 退出。

例如输入 `1` 和 `+12`，输出升高一个八度、时长与输入相同。输入 `r` 和 `+12`，输出同样升高一个八度，时长缩短为一半。保时长移调时，PGHI 先按 `时长倍率 × 音高倍率` 拉伸，再用 Kaiser 窗 sinc 重采样器输出。

结果为 Float32 WAV，保留采样率和声道数。文件名形如 `drumloop_pghi_2x.wav`、`drumloop_pghi_1x_p+12st.wav` 或 `drumloop_resampled_p+12st.wav`，同名时自动添加编号。

默认使用直接相位差分、全 STFT 幅度优先堆和双向瞬态锚点。程序在内存中复用传播路径、分析相位步长、瞬态定位和 chirp 拟合。复数频谱缓存上限为 128 MiB，超过时在合成中按帧从原音频取频谱。多声道共用相位修正，并以 WOLA 按窗平方和归一化。全图积分的内存和处理时间随文件时长增长。Gaussian 模式适合孤立扫频和单音。
