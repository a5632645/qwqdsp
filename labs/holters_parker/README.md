# holters_parker — Holters–Parker 非整数倍重采样

本目录研究论文 *A Combined Model for a Bucket Brigade Device and its Input and Output Filters*
(DAFx-18) 里的**重采样内核**，并把结论落回仓库（`resample_iir.hpp` / `elliptic_blep.hpp` /
`iir_design.hpp`）。

- 算法本身的逐步描述、推导与约束 → **[algorithm.md](algorithm.md)**
- 本仓库实现与论文的逐式对应、实测 → **[repo_impl_notes.md](repo_impl_notes.md)**
- 本文件 = 索引 + 结论摘要 + 跑法。

## 结论速览

1. **算法 = 论文 §3.1 的输入滤波器 + §3 的 modified impulse-invariant transform**（式 (2)–(11)）：
   单极点状态迭代 `X ← X·exp(λ) + c·x[k]` + 在任意分数时刻取值 `Σ Re(X·exp(λ·d))`，
   **不需要插值器**。§3.2 的输出滤波器（矩形保持形式）是另一半，本仓库没用。
2. **原型必须严格真**：`d = |H(j∞)| ≠ 0` 会让镜像以 `d` 的电平通过，且 `h(t)` 含 `δ`。
   偶数阶椭圆的修正已实现为 `IIRDesign::Elliptic(..., even_order_modify = true)`
   （`IIRDesign::EllipticEvenOrderModify`，见 [algorithm.md](algorithm.md) §9.1）。
3. **带直接项的正确实现**（§9.4）：状态迭代里每个输入样本注入**原始样本**（不乘任何上采样系数，
   §9.5），然后在"输出第 `q·n` 点恰好落在输入第 `p·n` 点"处叠加 `d·x[p·n]`
   （`f_in : f_out = p : q` 既约）。
4. **改变截止频率**（§9.3）：极点与留数**同乘** `scale = w/fpass`，直接项 `d` **不变**。
5. **精度主导项**：重采样后的杂散通常由原型在**镜像频率**处的衰减决定（本例第一个镜像落在
   过渡带 → −46 dB），而不是直接项；要压它得提高阶数/`rs`。
6. **NonIntegerSRC** 可当系数表与 SIMD 写法的参考，但它的滤波部分接线是失效的
   （活动分支的 `x` 被写入后从未读取，实测输出 ≈ 零阶保持），见
   [NonIntegerSRC_notes.md](NonIntegerSRC_notes.md)。

## 论文

| | |
|---|---|
| 标题 | *A Combined Model for a Bucket Brigade Device and its Input and Output Filters* |
| 作者 | Martin Holters（Helmut Schmidt University）、Julian D. Parker（Native Instruments GmbH） |
| 会议 | DAFx-18，第 21 届，Aveiro, Portugal，2018-09-04～08，pp. DAFx-11–16 |
| 官方页 | <https://www.dafx.de/paper-archive/2018/papers/DAFx2018_paper_12.pdf>（全文） |
| 作者页 | <https://www.hsu-hh.de/ant/en/team/martin-holters/dafx2018-bbd>（含试听样例） |

> 论文全文只引官方链接（DAFx archive 开放获取），**不在仓库内放副本**。

## 核心思想

BBD 就是一个**定长、可变采样率**的延迟线。输入侧要先从音频率 `f_s` 采到 BBD 时钟率 `f_BBD`，
输出侧再从 `f_BBD` 回到 `f_s`。通常这两步要额外插值，论文的要点是：

- 用**电路的输入低通 `H_in(s)` / 输出低通 `H_out(s)`** 兼作抗混叠与重构滤波器，
  经「modified impulse-invariant transform」（对部分分式展开的冲激响应直接采样）离散化；
- 因为变换保留了连续时间极点 `p_m` 的**分数采样间隔**信息，取某个时刻的样本不再需要
  插值，只要把各一阶子系统 `p̄_m = e^{p_m T_s}` 递归到对应整数拍、再乘一个依赖分数
  `d_n` 的权重即可（论文式 (7)–(11) 输入侧、(19)–(25) 输出侧）；
- 好处：不引入额外插值滤波带来的幅频失配，也没有「插值滤波器随延迟变化」导致的调制失真。

- 共轭极点对可合并为实系数二阶子系统（式 (26)–(33)）；论文也提到可用多项式近似或查表算 `b` 系数。
- 论文 Table 1 给了 Juno-60 chorus 的 `H_in` / `H_out`（各为 1 阶高通 + 5 阶低通，模型只用 5 阶低通）
  的极点 / 留数：`H_in` 1 实极点 + 2 对共轭，`H_out` 1 实极点 + 2 对共轭。
- 稳态近似（式 (34)–(38)）：常数 `f_BBD` 时 BBD 等效为
  `H_BBD(iω) = e^{-iω N/(2 f_BBD)} · sinc(ω/(2π f_BBD))`，`N` 为级数。
- 论文验证：正弦对照解析式；`f_BBD` 阶跃时输出连续（普通数字延迟线会跳变）；与真机 Juno-60
  录音对比，44.1 kHz 下多出的混叠可用**过采样**消除（计算量大部分在 BBD 时钟率，与 `f_s` 无关）。

## 与仓库现有代码的关系

| 位置 | 关系 |
|---|---|
| `include/qwqdsp/fx/resample_iir.hpp` | 已实现的 HP 重采样器（离线、`source_fs`→`target_fs`），插值滤波器用 `EllipticBlep` |
| `include/qwqdsp/fx/resample_iir_dynamic.hpp` | 同上的**流式/可变比率**版本（`SetRatio` / `SetPitchShift`，`Push`/`Read`）。⚠ 调用了不存在的 `EllipticBlep::SetCutoff`，一旦实例化就编译失败（详见 [repo_impl_notes.md](repo_impl_notes.md) 末尾） |
| `include/qwqdsp/fx/elliptic_blep.hpp` | 实际用的滤波器：`EllipticBlep`，把论文的「部分分式 + 分数幂 `p̄^{d_n}`」做成 partial LUT（`SetScale` 改截止、`Get(fraction)` 取分数拍样本） |
| `include/qwqdsp/filter/iir_design.hpp` | 原型设计器；本目录促成的改动：椭圆原型的**偶数阶修正** `even_order_modify`（§9.1） |
| `notebooks/holters_parker_coeff.ipynb` | 生成 `TCoeffs`（`fpass`/`fstop`/极点/留数）的设计脚本，四族 ellip/cheby1/cheby2/butter 可选 |

> 仓库现有实现只保留了**低通**部分（注释写明「移除了高通滤波器系数」）；论文原始模型还含
> 1 阶高通（偏置，与低通**级联**，所以整体仍严格真）与 BBD 增益 2.3 dB。

## 文件

| 文件 | 说明 |
|---|---|
| `algorithm.md` | **算法描述**：式 (1)–(8)、伪代码、截止约束、代价/误差；§9 专讲直接项（§9.1 偶数阶修正、§9.2 两个模型、§9.3 缩放规则、§9.4 正确实现、§9.5 输入是否乘上采样倍数） |
| `repo_impl_notes.md` | **本仓库** `ResampleIIR`/`EllipticBlep` 对应论文哪些式子、与论文的差异、实测数据 |
| `direct_term_study.py` | **直接项四步研究**：step1 偶数阶椭圆分解/重建 → step2 类连续频谱（线性轴 + 全带宽对数轴）→ step3 抽取到低速率 → step4 直接构建 HP 重采样器 + 精确命中点直接项 |
| `plot_direct_term.py` | 三个比率的直接项图（对应 §9.2 的"两个模型"对比） |
| `plot_resample_spectrum.py` | 三个比率的输入/输出频谱（用仓库真实 `BestCoeffs` 的系数 + 模型） |
| `NonIntegerSRC_notes.md` | NonIntegerSRC 源码阅读笔记：与论文公式的映射、系数核对、失效点、实测 |
| `noninteger_src_repro.py` | 上面笔记的复现脚本（NonIntegerSRC 数值路径的逐行复刻 + 行为测量） |
| `audio_developer_guide.md` | 面向普通音频开发者的 Holters–Parker 任意倍重采样算法总结与实战避坑指南 |
| `output/` | 生成的图片（本地，不在版本管理内） |
| `README.md` | 本文件 |

## 怎么跑

```bash
# 四步研究（每一步都会把图写进 output/）
python qwqdsp/labs/holters_parker/direct_term_study.py 1     # 分解/重建
python qwqdsp/labs/holters_parker/direct_term_study.py 2     # 类连续频谱（10 周期 + 全带宽对数轴）
python qwqdsp/labs/holters_parker/direct_term_study.py 3     # 抽取到低速率
python qwqdsp/labs/holters_parker/direct_term_study.py 4     # HP 重采样器 + 直接项
python qwqdsp/labs/holters_parker/direct_term_study.py       # = all

python qwqdsp/labs/holters_parker/plot_direct_term.py        # 三比率直接项图
python qwqdsp/labs/holters_parker/plot_resample_spectrum.py  # 三比率输入/输出频谱
python qwqdsp/labs/holters_parker/noninteger_src_repro.py    # NonIntegerSRC 复刻
```

依赖：`numpy` / `scipy` / `matplotlib`（图中中文用 `Microsoft YaHei`）。参数都在各脚本顶部或
函数签名里（`fs_in / tone / n / upsample / fs_out`）。
