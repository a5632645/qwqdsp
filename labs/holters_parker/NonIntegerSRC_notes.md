# NonIntegerSRC 阅读笔记

对象：<https://github.com/jatinchowdhury18/NonIntegerSRC>（GPLv3，38★，最后一次看是 1 次提交量级的实验仓库）
目的：判断它作为 Holters–Parker 重采样参考实现的价值。复现脚本见 [noninteger_src_repro.py](noninteger_src_repro.py)。

## 结论速览

- 它把论文**输出滤波器那半边**（式 (18)–(25)）拿出来做成了通用 SRC：把输入信号当成 BBD 的
  零阶保持输出，用连续差分 `Δ(n)=y[n]−y[n−1]` 当阶跃幅度，用滤波器冲激响应的分段解析式
  直接在新采样率的时间点上求值。
- **但它不是一个可信的实现**：活动分支里 `outFilter->x` 被 `process()` 写入、从未被读取，
  输出用的是驱动项累加器 `xOutAccum`；`calcH0()` 定义了却从未调用；常量项硬编码成 1。
  实测该分支 ≈ 零阶保持 + 微小修正，96k→48k 时 30 kHz 输入完全不衰减（见下）。
- 仓库自带的 benchmark **只打印耗时**（`src_test.cpp` 里误差统计那两行被注释掉了），
  所以这类问题不会被它自己的测试发现。它宣传的是「比 libsamplerate 快 10–40×」，
  没说质量。
- 有价值的部分：**论文 Table 1 的系数被逐位抄进去了**，以及一套把 4 个复一阶滤波器
  塞进 2 个 `__m128` 的 SIMD 写法 + 快速 sin/cos；这两块可以借用，滤波器的接线要自己写。

## 文件 → 论文公式的映射

| 源码 | 论文 |
|---|---|
| `FilterSpec::iFiltRoot/iFiltPole` | Table 1 的 `H_in` 残差/极点 |
| `FilterSpec::oFiltRoot/oFiltPole` | Table 1 的 `H_out` 残差/极点 |
| `InputFilterBank`（只被 `#if 0` 分支用） | §3.1 输入滤波器，式 (2)–(11) |
| `OutputFilterBank`（活动分支） | §3.2 输出滤波器，式 (15)–(25) |
| `HPResampler::process` 的 `#if 0` 分支 | 式 (10)：`u_BBD = Σ_m g_in,m(d_n)·x_in,m(l_n)`，**读状态，形式正确** |
| `HPResampler::process` 的 `#else` 分支 | 式 (23)(25)：`x_out,m(k)=p̄ x(k−1)+Σ g Δ`，`y=H₀·y_old+Σ x`，**输出接错** |

### FilterSpec 与论文 Table 1 的逐项核对

| | 核对结果 |
|---|---|
| `iPole` / `oPole` | 与论文 `p2..p5` **完全一致**（`H_in`：-55482±25082i、-26292±59437i；`H_out`：-51468±21437i、-26276±59699i） |
| `oRoot` | 与论文 `H_out` 的 `r2..r5` **一致**，**但** -11256 那一对的实部符号相反（论文写 `+11256`，代码写 `-11256`） |
| `iRoot` | = 论文 `H_in` 的 `r2..r5` **整体乘 1/12.6271**（极点未缩放，所以是纯增益缩放） |
| 实极点对 | 论文的 `p1` / `r1`（`H_in`: −46580/251589，`H_out`: −176261/5092）**被丢弃**——`N_filt=4` 正好填满一个 SSE 寄存器 |
| DC 增益 | 按论文 Table 1 取 `−Σr/p`：`H_in≈0.855`、`H_out≈−1.833`（都不是 1）；用代码那 4 项是 `H_in≈0.360`、`H_out≈−2.607`。Table 1 显然不是 DC 归一化的 |

### 截止频率怎么定

`set_freq(sample_rate * 0.5)`，`originalCutoff` 是 9900（输入）/ 9500（输出）：

```
freqFactor = (fs/2) / originalCutoff
root_corr  = roots * freqFactor                       // 模拟域频率缩放
pole_corr  = exp(poles * freqFactor * Ts_filter)
```

即把论文那套固定 ~10 kHz 的 Juno-60 滤波器**缩放到 fs/2**当成通用 SRC 的滤波器。
注意 `sample_rate` 是**输入**采样率，所以截止被放在**输入** Nyquist；降采样时不会抗混叠
（本仓库 `resample_iir.hpp` 用的是 `min(source_fs, target_fs)/2`，是对的版本）。

## 时间与比率的约定（容易被变量名骗）

`ratio = 输出采样率 / 输入采样率`（`src_test.cpp`：48k→96k 用 2.0，96k→48k 用 0.5）。

```
Ts     = 1/fs               = 输入采样周期 T_in
Ts_in  = 1/(fs*ratio)       = 输出采样周期 T_out        （名字最反直觉的一个）
Ts_out = 1/(fs/ratio)       = ratio * T_in = ratio^2 * T_out
```

- 活动分支的内层 `while (tn < Ts) { ... tn += Ts_out; }`：每个输出样本消费
  `Ts/Ts_out = 1/ratio` 个输入样本 —— 这个**计数是对的**。
- 但 `outFilter->set_delta(Ts_out)` 把 `Ts_out` 当作 p̄ 的指数增量。滤波器的内部 `Ts = Ts_in = T_out`，
  所以真正的增量应当是 `T_in/T_out = ratio`，代码给的是 `Ts_out/T_out = ratio²` —— **差一个 ratio**。
  对照：输入分支用 `set_delta(Ts_in)`，而它的滤波器内部 `Ts = T_in`，增量 = `T_out/T_in = 1/ratio` ✓。
- `exp(i*angle*delta)` 又只保留相位（`|Aplus| ≡ 1`），把 `p̄^δ` 的模衰减丢掉了，
  所以 `Gcalc` 既不会随步数衰减、也不会携带 `p̄^{k−l_n}` 那种"越老越小"的记忆。

## 活动分支的数据流

```
每个输入样本 y：  Δ = y - y_old；y_old = y
                  Gcalc *= Aplus            // calcG()：把 d_n 前进一个步长
                  acc += Gcalc * Δ          // 驱动项 Σ g_out,m(d_n)·Δ(n)
每个输出样本：     x = p̄*x + acc             // process()：论文式 (23) 的递归  ← 结果没被用
                  out = y_old + Re(acc)     // ← 实际用的是 acc，不是 x
```

论文式 (25) 要求 `y(k) = H₀·y_BBD,old + Σ_m Re(x_out,m(k))`。源码里：

- `outFilter->x` 由 `process()` 更新后**再没被读过**；
- `calc_h0()`（= `−Σ Re(r/p)`，正是 H₀）定义了但从未被调用；
- `y_old` 的系数硬编码成 `1`。

即：递归状态 `x`、DC 项 `H₀` 两个都白算了。

## 复刻实测（`python noninteger_src_repro.py`）

(A) 降采样 96k→48k，输入 30 kHz（输出 Nyquist 24 kHz，必须滤掉；不滤会折到 18 kHz）：

| 变体 | 18 kHz 处幅度（满幅输入 = 1.0） |
|---|---|
| 活动分支（as-coded） | **1.0057** |
| 改成读状态 `x` | **1.0055** |
| `#if 0` 输入滤波器分支 | 0.9761 |
| 朴素 ZOH 抽取（对照） | 1.0000 |

→ 活动分支基本等于不抗混叠的零阶保持。`#if 0` 那个分支也只衰减 2%（因为截止放在输入 Nyquist，见上），
但**它的接线形式是对的**。

(B) 升采样 48k→96k，单音最小二乘拟合（幅度 / 残差 rms）：

| f | as-coded | 读状态 | ZOH |
|---|---|---|---|
| 1 kHz | amp 1.000 / resid 0.0248 | 1.000 / 0.0233 | 0.999 / 0.0231 |
| 5 kHz | 0.989 / 0.1233 | 1.009 / 0.1153 | 0.987 / 0.1152 |
| 10 kHz | 0.954 / 0.2421 | 1.024 / 0.2229 | 0.947 / 0.2273 |

→ as-coded 全程贴着 ZOH；残差里绝大部分是 ZOH 的镜像，滤波器没起到重构作用。

## 可以直接借用的部分

1. **论文 Table 1 的系数**（`FilterSpec`）——已经过逐项核对；本仓库要做 Holters–Parker 建模时可直接用，
   注意补回被丢掉的实极点对，并把它当成"未归一化"的绝对增益。
2. **4 复一阶 / 2 个 `__m128` 的 SIMD 布局**：`SSEComplex` + `vSum` + `fastsinSSE`/`fastcosSSE`
   （多项式 + `clampToPiRangeSSE`），以及"用递推乘法代替逐样本 `pow`"的算 `p̄^{d_n}` 的思路。
3. **`#if 0` 输入滤波器分支**（式 (10)）的接线方式：滤波递归跑在输入率，输出点用
   `Σ_m g_in,m(d_n)·x_in,m(l_n)` 求值——这是"不插值就能在任意时刻取样"的关键，形式是对的。

## 不能照抄的部分

- 活动分支的滤波（见上）：状态与 H₀ 都没接进输出。
- `set_delta` 的标度（差一个 ratio）。
- 截止放在 `fs/2`（应为 `min(fs_in,fs_out)/2`）。
- 只用了输入/输出滤波器中的**一个**；论文模型是两个都用（输入滤波抗混叠、输出滤波抗镜像）。
