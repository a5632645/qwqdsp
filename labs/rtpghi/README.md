# 双向 RTPGHI 离线实验

线性 chirp 的目标梯度、窗宽效应和输入 STFT 矩估计见 [时间拉伸后的目标相位梯度](chirp_gradient.md)。

## 相位差分

博客的 `propagatePghi()` 用当前帧减前一帧，减去 bin 中心频率对应的相位步长后解缠，再换算到合成 hop。频率方向直接用相邻 bin 的相位差，按伸缩率换算。自然顺序窗会先扣除窗中心引入的线性相位，传播后再补回。

《Phase Vocoder Done Right》第 III 节先分别解缠前向与后向差分，再平均成节点的中心梯度。时间边和频率边都用两端梯度的平均值做梯形积分。中心时间梯度需要后一帧。算法 1 的旧帧入堆条件按**当前** bin 的幅度判定。旧 bin 即使很弱，同一频带缺少更可靠的入口时仍可能先传入，再沿新帧的强 bin 扩散。实验中的 `direct` 与 `centered` 共用堆传播和 WOLA，单独比较两套差分。实测中心差分的收益随信号和落点而变。[论文，第 III 节与算法 1](https://arxiv.org/html/2202.07382v1)

## 瞬态相位与传播

1. 预先计算整个信号的中心化 STFT、幅度、相位和梯度。瞬态位置可人工给定，也可用短窗宽带正谱通量检测，再在波形上细化。
2. 选最接近起音的分析帧。在短窗起音谱中找相对起音前上升、且足以解释长窗幅度的 bin。固定这些 bin 的合成相位。在长窗分析相位上加延时斜坡 `−2πkΔ/M`，其中 `Δ` 是起音相对合成帧中心和分析帧中心的偏移之差。这样，偏离 hop 网格的冲激也能对齐到 `round(α·onset)`。
3. 从所有锚点入堆，在有合成窗重叠的帧范围内积分。路径按沿途最弱 bin 的幅度排序，相同优先级先处理较短的路径。沿正方向加边梯度，沿负方向减边梯度。锚点相位始终固定。最后用逐帧 RTPGHI 补齐其余 bin，并把已求得的相位续接到后续帧。

离线 PGHI 已有从已知相位掩码出发的积分。本实验以相位声码器的分析相位差分构造梯度，并以映射到合成时间的瞬态相位作锚点。[PGHI 论文，第 IV-B 节](https://ltfat.org/notes/ltfatnote040.pdf)

| 方法 | 差分与相位来源 |
| --- | --- |
| `direct` | 博客式直接边差分，逐帧单向积分 |
| `centered` | 论文式中心差分与梯形积分，逐帧单向积分 |
| `reset` | 同 `centered`，在指定 bin 固定相位，再逐帧向前积分 |
| `forward` | 同一批锚点在时频网格上向前传播 |
| `bidirectional` | 同 `forward`，加入向过去的时间边 |

## 实验结果

48 kHz，2048 点周期 Hann，4096 点 FFT，合成 hop 256，默认传播半径 7 帧（约 37.3 ms）。表中为 2× 时的前振铃能量占比，分母是事件附近的总输出能量。噪声突发由单点起音和 7 ms 衰减噪声组成。

| 输入 | `direct` | `centered` | `reset` | `forward` | `bidirectional` |
| --- | ---: | ---: | ---: | ---: | ---: |
| 单冲激 | −3.3 dB | −5.7 dB | −5.8 dB | −5.8 dB | −150.0 dB |
| 噪声突发 | −9.8 dB | −12.1 dB | −14.8 dB | −14.7 dB | −29.8 dB |

- **落点扫描。** 扫描 8 个冲激位置，`bidirectional` 在两个拉伸率下的峰值都落在目标采样点。在实际 1.497× 下，`forward` 的前振铃为 −24.9 至 −6.9 dB，`bidirectional` 为 −111.2 至 −55.5 dB。2× 下分别为 −17.1 至 −6.2 dB 和 −150.0 至 −88.4 dB。−150 dB 是统计下限。
- **半径扫描。** 在默认冲激落点、2× 下，双向传播半径为 0、1、2、3、7 帧时，前振铃为 −5.8、−17.4、−32.2、−54.1、−150.0 dB。噪声突发在 2–3 帧时已接近 7 帧的结果，2× 下尾部能量为 −31.1 dB，论文式单向为 −32.7 dB。
- **弱 bin。** 冲激落点扫描中，`centered` 在起音附近有 1–2 次直接从幅度低于目标十分之一的旧 bin 取相位。单次薄弱入口就能沿强频带扩散。双向传播对冲激的改善同时包含相位重新定标和更换积分路径。
- **持续音和混合信号。** 合成测试的自动检测命中了全部 4 个起音，在持续音对照中未产生锚点。混合信号在 1.497× 下，双向传播后第二个起音附近的最大峰仍提前 477 个采样点（约 9.9 ms）。混合输出的前后能量还包含持续音本身。两个相邻锚点可以通过重叠的传播区域竞争，按幅度传播也可能改动持续音的相位。
- **鼓循环。** 仓库鼓循环的前 2 秒在 2× 下自动检测到 10 个候选起音，共固定 3040 个 bin。STFT 相对残差为 `centered` 0.391、`forward` 0.434、`bidirectional` 0.394。音频无人工起音标注，检测准确率与听感尚未验证。
- **WOLA 一致性。** 各帧冲激完全对齐时，WOLA 的峰值增益为 `Σ(w_analysis·w_synthesis) / Σ(w_synthesis²)`，窗值分别按冲激在分析帧和合成帧内的位置取样。2× 时，双向传播的主峰在目标样点，幅度为 0.986，输入为 0.8。沿用分析 STFT 的幅度使重分析与目标 STFT 的相对残差仍约为 0.49。

## 全 STFT 堆

把上一节的积分域从瞬态附近扩大到整张 STFT。每个 bin 有前后两条时间边和上下两条频率边。堆保存已发现的传播前沿。堆空后，从剩余 bin 中选幅度最大的一个，以它的分析相位初始化，再积分下一个连通区域。普通 bin 的幅度门限取整张 STFT 峰值的 `tolerance` 倍。带锚点的版本先从所有固定相位传播，再补其余连通区域。原始离线 PGHI 已采用全时频堆来从幅度重建相位。这里共用上一节的中心相位差分和伸缩后的合成步长。 [PGHI 论文，算法 1 与第 IV-B 节](https://ltfat.org/notes/ltfatnote040.pdf)

| 方法 | 堆优先级 | 瞬态锚点 |
| --- | --- | --- |
| `global-magnitude` | 当前 bin 幅度，首次到达即定相位 | 无 |
| `global` | 路径最小幅度优先，同优先级先走步数短的路径 | 无 |
| `global-magnitude-anchored` | 当前 bin 幅度 | 有 |
| `global-anchored` | 路径最小幅度 | 有 |

中心差分得到的时间边与频率边绕一格相加，可能得到非零的模 `2π` 残差。堆选出的路径因此会影响最终相位，尤其在长持续音和多个锚点共享的区域。CSV 记录自动种子数、全局积分 bin 数、锚点可达 bin 数，以及其中落在原局部半径之外的数量。

同样设置下的 2× 结果。STFT 残差衡量 WOLA 输出重分析后与目标复数 STFT 的差异。

| 方法 | 冲激前振铃 | 突发前振铃 | 混合信号残差 | 持续音残差 |
| --- | ---: | ---: | ---: | ---: |
| `centered` | −5.7 dB | −12.1 dB | 0.133 | 0.119 |
| `bidirectional` | −150.0 dB | −29.8 dB | 0.203 | 0.119 |
| `global-magnitude` | −3.6 dB | −15.0 dB | 0.133 | 0.119 |
| `global` | −3.6 dB | −14.9 dB | 0.264 | 0.372 |
| `global-magnitude-anchored` | −150.0 dB | −27.5 dB | 0.133 | 0.119 |
| `global-anchored` | −150.0 dB | −29.8 dB | 0.602 | 0.372 |

- 无锚点的两个全局堆都把冲激的主峰放晚了 1 个样点。自动选中的种子继承分析相位，在该样例里与冲激按目标时刻平移后的相位相差 −1.20 rad。扩大积分域后仍需要给瞬态定相。
- 两种全局锚定方式在 2× 混合样本中都使两个起音附近的主峰对齐，而且都传播到了原局部半径外的 46,201 个 bin。两者的 STFT 残差从 0.133 到 0.602，说明同一覆盖范围内，路径选择也会改变合成结果。1.497× 下，第二次起音附近最大峰提前约 9.9 ms 的问题仍然存在。
- 将门限从 `1e-6` 提到 `1e-3`，路径最小幅度版在混合样本中传播到半径外的 bin 减至 6,585 个，STFT 残差仍约 0.602。
- 鼓循环的前 2 秒在 2× 下，`centered`、`bidirectional`、`global-magnitude-anchored` 的 STFT 残差分别为 0.391、0.394、0.396。此样本使用自动检测，听感与起音准确率仍需单独评估。

## 直接边差分

`--propagation-edges direct` 把重置、局部双向和全局堆都切到直接差分。每条时间边用相邻帧的解缠相位差，每条频率边用相邻 bin 的主值相位差，分别乘以实际伸缩比 `α = Hs/Ha`。逐帧 `direct` 和 `centered` 仍提供两套差分的基线。

直接边等于分析端点相位之差加整数圈 `2πk`。绕闭环时这些端点差相消，伸缩后只剩 `2π·α·Σk`。因此，**整数 `α` 下直接边绕任意闭环都在模 `2π` 意义下闭合**，等价于对相位场 `α·φ_analysis` 取边差。同一个连通分量只有一个种子时，堆走哪条路径都得到相同相位。多个硬锚点仍可能给同一分量设定不同的相位常数。非整数 `α` 下，解缠的整圈经过缩放后可能变成非整圈。中心差分还会平均相邻边，因此不具备整数倍率的闭合保证。闭合度描述目标相位场，WOLA 重分析残差还取决于相位和幅度能否共同合成一致的 STFT。

CSV 的 `loop_rms_rad` 是时频小格的模 `2π` 局部闭环误差，取四角中最弱幅度的平方加权。四角都超过全局门限才参与统计，`loop_nonclosing_fraction` 是其中闭环误差大于 `1e-8` rad 的比例。

- **2× 混合样本。** 闭环 RMS 从中心差分的 0.129 rad 降到直接差分的约 `6e-15` rad。`global-anchored` 的 STFT 残差从 0.602 降至 0.109，`global-magnitude-anchored` 从 0.133 降至 0.109，两次起音的主峰仍对齐。四个合成样本的无锚点全局输出在两种堆排序下相差不超过 `1.2e-13` 样本幅值。
- **1.497×。** 混合样本的闭环 RMS 从 0.132 降至 0.055 rad，`global-anchored` 的 STFT 残差从 0.510 降至 0.103。持续音的路径最弱幅度版残差从 0.196 升至 0.380，直接差分的当前 bin 幅度版为 0.059。非整数倍率下，堆排序仍会影响结果。
- **瞬态。** 直接差分保留了 2× 锚定冲激的样点对齐和 −150 dB 前振铃下限。噪声突发的全局锚定结果变化很小，路径最弱幅度版在 2× 下两套差分均约 −29.8 dB。1.497× 混合样本的第二次起音附近最大峰仍提前 477 个样点。
- **真实鼓循环。** 2× 时，`global-magnitude-anchored` 的残差由 0.3955 降至 0.3885，`global-anchored` 由 0.3922 降至 0.3883。局部 `forward` 的残差由 0.4336 升至 0.4448。听感尚未验证。

## 复现

从工作区根目录运行，需要 Python 3.10+、NumPy、SciPy 和 Matplotlib。WAV 使用浮点单声道，保留输出电平供比较。所有实验在内存中保存整个 STFT，适合短片段。

```powershell
python qwqdsp/labs/rtpghi/experiment.py
python qwqdsp/labs/rtpghi/experiment.py --impulse-offset 0 21 42 63 84 105 126 147 --anchors known --no-plot --out temp/rtpghi/impulse-sweep
python qwqdsp/labs/rtpghi/experiment.py --case impulse burst --ratio 2 --radius 3 --anchors known --out temp/rtpghi/radius-3
python qwqdsp/labs/rtpghi/experiment.py --input qwqdsp/work_dir/input/sine_with_drum.wav --ratio 2 --out temp/rtpghi/wav
python qwqdsp/labs/rtpghi/experiment.py --global-heap --anchors known --out temp/rtpghi/global
python qwqdsp/labs/rtpghi/experiment.py --global-heap --propagation-edges direct --anchors known --out temp/rtpghi/global-direct
python -m unittest discover -s qwqdsp/labs/rtpghi -v
```

默认输出 `temp/rtpghi/runs` 下的 WAV、波形和锚点覆盖图、`metrics.csv`，直接边传播默认写到 `temp/rtpghi/runs-direct`。`--onset-ms` 可给 WAV 提供人工标注，`--anchors known` 或 `detected` 可单独运行一种锚点来源。CSV 记录起音位置、峰值偏移、前后能量、弱旧帧入口数和 STFT 重分析残差。对自动检测的 WAV，前后能量以检测位置和 7 ms 名义持续时间计，适合作为定位线索。测试覆盖相位解缠、整数与非整数倍闭环、双向边符号、锚点竞争、WOLA 恒等重构和离网格冲激。
