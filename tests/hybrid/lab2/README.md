# lab2 — 无窗函数 DFT 实时音乐分析（NC 方法）

本文件夹是用 Python 复现论文 **《Window Function-less DFT with Reduced Noise
and Latency for Real-Time Music Analysis》**（arXiv:2410.07982v3，Biesinger 等，
京都大学）的算法实现，输入一个 `.wav`，输出一张上/下双面板的频谱对比图，
并附一个针对「时间重分配」有效性的对照实验。

## 论文核心行为

传统 FFT/DFT 输出线性间隔的频率区间，不适合音乐分析；且必须加窗函数抑制旁瓣，
加窗又会加宽主瓣、劣化频率分辨率。本文提出「邻域频谱分量合成」(NC) 方法：

- 把每个输出 bin 的中心频率按十二平均律指数公式铺设，直接对应音高：
  `f(n) = 440 Hz × 2^((n-69)/12)`。
- 每个 bin 取左右两个相邻 DFT 分量，用其实/虚部做点积得到 NC 幅度：
  `max(0, -(Re_L·Re_R + Im_L·Im_R))`。
- 这样无需窗函数即压低旁瓣（图例中 STFT 可见的频谱泄漏条纹几乎消失），
  主瓣更窄、频分辨率翻倍，且窗长可减半，提升时间分辨率与延迟表现。

## 文件说明

| 文件 | 说明 |
|------|------|
| `nc_dft.py` | 算法核心：bin 构造、NC 幅度计算(FFT 卷积等价)、时间抽取与 IIR 平滑 |
| `run.py` | 驱动脚本：读 wav → 上(Blackman-Harris STFT)/下(NC) 双面板对比图 |
| `run_sliding_dft.py` | 改进版：用**递归滑动 DFT** 逐样本流式计算，输出**无需相位校正**(式8) |
| `verify_time_reassignment.py` | 对照实验：验证「时间重分配」对低频瞬态是否真的聚焦(见下文) |
| `output/` | 生成的图片(已 gitignore，本地生成) |
| `Window Function-less DFT*` | 论文 HTML 与 Fig.4 参考图(已 gitignore，本地对照) |

`nc_dft.py` 中每个公式都标注了论文对应的编号（(1)(5)(6)(7)(8)），方便对照。

## 两种实现的区别

论文原实现(IV-B/IV-C 节)用「全局累积相位」的参考正弦波，每样本把相位
θ_L/θ_R 累加，因此输出前必须用旋转矩阵 R(θ) 校正(式8)。本文件夹提供两种
等价实现：

- **`run.py`**(基于 `nc_dft.py`)：用 FFT 卷积等价实现逐样本累加器，结果与
  论文算法一致(离线处理)。
- **`run_sliding_dft.py`**：改用**递归滑动 DFT** `X(n) = W·X(n-1) + x[n] − x[n−N]·W^N`，
  其中 `W = exp(-j·2π·f/F_S)` 恒定。因为 X(n) 恒等于「以当前窗口起点为相位
  零点的窗口 DFT」，相位不随时间漂移，**任意时刻读出的 X 已正确，无需相位
  校正**，且逐样本流式、每样本仅 O(m) 复杂度(论文 IV-E 节)。

两种实现在 drumloop.wav 上的 NC 输出已校验为**逐点完全一致(最大差异 0)**。

## 运行

测试输入在 `../work_dir/input/`（如 `drumloop.wav`、`sine.wav`、`sweep.wav`、
`carry.wav`）。默认处理 `drumloop.wav`：

```bash
python run.py                       # 默认 drumloop.wav → output/drumloop_combined.png
python run.py sine --out output/sine_combined.png
python run.py sweep --out output/sweep_combined.png
python run.py path/to/any.wav --out output/any.png

# 递归滑动 DFT 版(无相位校正)
python run_sliding_dft.py           # 默认 drumloop.wav → output/drumloop_sliding.png
python run_sliding_dft.py sine --out output/sine_sliding.png

# 时间重分配有效性对照实验
python verify_time_reassignment.py  # → output/tr_plain_vs_reassign.png
```

参数 `<wav>` 若在 `../work_dir/input/` 下，可直接写文件名（可省略 `.wav`），
否则写文件路径；`--out` 指定输出图片。

## 依赖

Python 3.11+，`numpy`、`scipy`、`soundfile`、`matplotlib`。

## 结果说明

上方为 Blackman-Harris 加窗 STFT 参照，下方为本文 NC 方法；纵轴为对数频率
27.5–3520 Hz，幅度为归一化 dB（0 最强，-60 为显示下限）。不同音频的对比行为：
- `sine.wav`（单音）：NC 只留下一根干净的细线，其余被压到 -60 dB 地板；
  STFT 同频带但边缘有扩散和伪影。
- `sweep.wav`（扫频）：NC 输出一根跟随扫频的细线；STFT 带更粗且有窗函数条纹。
- `drumloop.wav`（鼓点）：STFT 有丰富宽带瞬态，NC 主要对低频鼓点响应。

## 时间重分配验证结论

`verify_time_reassignment.py` 构造一颗 60Hz 低频瞬态，对比「plain(无时间重分配)」
与「reassign(用左右两 DFT 分量交叉相位 `arg(X_R·conj(X_L))` 算群延迟，再折算成
整数列偏移写入持久缓冲)」两种频谱：

- **plain**：一个轮廓清晰的连续脉冲块，能量集中在其真实 onset 附近。
- **reassign**：连续块被撕成多条竖直细条，低频出现空隙，高频冒出伪影。

**结论**：群延迟偏移的时间重分配对低频瞬态**不聚焦，反而引入伪影**。原因是
windowless NC 是逐箱滑动 DFT，左右分量交叉相位算出的群延迟不稳定、多成分下
被错误地分数列，而非聚焦到真实时刻。因此该方案**未采用**——C++ 侧的
`WindowlessNcTimeFrame` 已被移除，仅保留本实验作为证据。
