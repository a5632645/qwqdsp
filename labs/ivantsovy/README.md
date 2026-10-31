# ivantsovy — 理想双线性变换 §4 一般形式的数值验证

验证 Yuriy Ivantsov, *On the Ideal Bilinear and Biquadratic Digital Filter*
（2026，<https://ivantsovy.com/research/paper1.pdf>）第 4 节的一般形式：
把连续传递函数的分子、分母分解成常数项为 1 的实系数一次 / 二次因子，
逐因子映射到 z 域；分子阶数亏欠（`s → ∞` 的零点）按 `a = 0` 的退化因子处理，
`s = 0` 的零点按 (4.2) 第二行映射。因此 LP / HP / BP / BS 四种经典类型都覆盖。
论文第 4 节的 Python 数值验证，不改动库代码。

本 lab 对照两种配置：

- **默认**：`σ = π√(2/3) ≈ 2.5651`，不做预畸变。论文幅度最优区间的上端点，
  整体幅度误差最小。
- **对照**：`σ = 2`，做预畸变。参考频率 `f_ref` 处精确匹配，Nyquist 的精确
  匹配点则移到 `f_ref·ω_d/σ`（本例 23.67 kHz），不再是 24 kHz。

## 配置

| 项目 | 值 |
|---|---|
| 采样率 | 48 kHz |
| 连续域参考频率 `f_ref` | 5 kHz |
| 阶数 | 主实验 8；四种类型图用 4 阶原型（BP / BS 为 8 阶） |
| 默认 `σ` | `π√(2/3) ≈ 2.5651`，无预畸变 |
| 对照 | `σ = 2`，预畸变 |
| Butterworth | 截止 5 kHz |
| Chebyshev II | 阻带边沿 5 kHz，阻带衰减 80 dB |
| High shelf | N=8，f0 = 1 kHz，±12 dB；Chebyshev II ripple=1 dB |

连续域按 `S = s/(2π f_ref)` 归一化。离散响应按因子串联求值，不展开高阶 z 多项式
（高阶展开系数病态）；论文公式与实现的对应关系见 `ivantsovy.py` 的模块说明。

## 默认配置的结果

| 指标 | Butterworth N=8 | Chebyshev II N=8, 80 dB |
|---|---|---|
| 全频段 `max(abs(Hd − Ha))`（线性） | `1.59e-3` | `1.47e-4` |
| 模拟幅度 ≥ −20 dB 内 `max(abs(Δ))` | 0.085 dB | 0.007 dB |
| ≥ −40 dB 内 | 0.279 dB | 0.020 dB |
| ≥ −60 dB 内 | 0.915 dB | 0.057 dB |
| ≥ −80 dB 内 | 3.05 dB | 0.214 dB |
| 最大数字极点模 | 0.880 | 0.943 |
| 数字零点 | 8 个阶数亏欠零点，`z = −0.1010` | 6 个在单位圆上，2 个在实轴 `z = −0.2838` |
| Nyquist 值 / 模拟 Nyquist 值 | `2.60e-5 / 3.55e-6` | `5.59e-5 / 1.08e-5` |

图（同一坐标系，实线模拟、虚线离散）：

- `output/butterworth_N8.png`
- `output/chebyshev2_N8.png`

## 预畸变：截止频率 vs Nyquist

截止频率和 `fs/σ` 的精确匹配互斥。预畸变把设计用的 `ω` 从 `fs/f_ref` 换成

```
ω_d = sqrt(σ² + π² cot²(π / (fs/f_ref)))
```

（即论文 (2.10) 的等价形式），于是 `F(ω_d) = f_ref`。数值上验证了：

- Butterworth（任意阶）：`Hd(f_ref) = 1/√2`，机器精度；
- Chebyshev II：`Hd(f_ref) = Ha(f_ref)`（阻带边沿电平），机器精度；
- Nyquist 的精确匹配点变为 `f_ref·ω_d/σ`：`σ=2` 时是 23.67 kHz，`σ=π√(2/3)`
  时是 18.72 kHz。因此 `σ=2 + 预畸变` 不再精确匹配 `A(24 kHz)`，而是
  `A(23.67 kHz)`。

| 配置 | Butter maxlin | Butter ≥−60 dB | Cheb2 maxlin | Cheb2 ≥−60 dB | `f_ref` 处 Δ（Butter / Cheb2） | Nyquist 匹配点 |
|---|---|---|---|---|---|---|
| σ=2.565，无预畸变 | `1.59e-3` | 0.915 dB | `1.47e-4` | 0.057 dB | +0.014 / +0.214 dB | 18.71 kHz |
| σ=2，无预畸变 | `4.66e-2` | 4.867 dB | `1.62e-2` | 2.432 dB | −0.506 / −12.30 dB | 24.00 kHz（精确） |
| σ=2，预畸变 | `1.00e-2` | 3.908 dB | `3.51e-2` | 0.708 dB | 0 / 0 dB | 23.67 kHz |
| σ=2.565，预畸变 | `7.37e-4` | 0.888 dB | `1.30e-3` | 0.029 dB | 0 / 0 dB | 18.72 kHz |

`σ=2 + 预畸变` 的图：

- `output/butterworth_N8_prewarp_sigma2.png`
- `output/chebyshev2_N8_prewarp_sigma2.png`

## 四种经典类型（Butterworth / Chebyshev II）

§4 的因子规则按零点来源分三类：

| 零点 | 映射 |
|---|---|
| 有限零点（含 ±j，如 Chebyshev II / band-stop） | 二次因子 (4.4) |
| 原点零点 `S^m`（HP / BP） | (4.2) 第二行：每个 `S → (1−z^-1)·ω/(2π)` |
| `s → ∞` 的阶数亏欠 | `a = 0` 的退化一次因子 (4.7) |

因此 LP / HP / BP / BS 都能覆盖；模拟增益原样保留。二阶特例与论文
Table 2.1 / 3.1 逐系数一致（`check_impl.py`）。两族原型（BP / BS 边沿 4–6 kHz，
σ=π√(2/3) + 预畸变）的幅度对照：

| 类型 | Butterworth N=4 maxlin / ≥−20 dB | Chebyshev II N=8, 80 dB maxlin / ≥−20 dB |
|---|---|---|
| LP | `3.2e-3` / 0.124 dB | `1.3e-3` / 0.029 dB |
| HP | `3.0e-4` / 0.012 dB | `1.5e-2` / 0.225 dB |
| BP | `3.5e-3` / 0.113 dB | `5.8e-3` / 0.159 dB |
| BS | `3.1e-3` / 0.058 dB | `1.2e-2` / 0.260 dB |

图：`output/four_types.png`（Butterworth）、`output/four_types_cheby2.png`
（Chebyshev II，带通 / 带阻为 8→16 阶）。`check_impl.py` 另有无预畸变的
N=4 / N=8 / N=7 高阶对照，max linear 分别不超过 `6.4e-3` / `1.9e-2` /
`2.7e-2`。奇数阶 Chebyshev II 的高通 / 带通会出现原点零点，同样由 (4.2)
第二行处理。

一个必须注意的增益细节：scipy 的 zpk 增益不能直接当成本实现的 `gain`。
使用 `normalized_gain` 换算成“原点零点 + 常数项 1 因子”形式后，模拟响应才与
`signal.freqs_zpk` 一致。Chebyshev II 低通的 zpk 增益是 `1e-4`；直接使用会
整体低 80 dB，高通 / 带通也会差一个常数。Chebyshev II 高通的直流阻带是
−80 dB（不是 z=1 处的零点），带阻的零陷和纹波也都保留。

s = 0 处的**极点**（积分器）没有对应规则，会被显式拒绝；经典四种类型的极点
都不在原点。

## High shelf（最严格）

scipy 没有 high-shelf 原型；本 lab 复刻 `iir_design_shelf.hpp` 的构造：

- Butterworth：`z = p/s`、`p = p·s`，`s = (10^(boost/20))^(1/N)`，于是
  `H(0)=1`、`H(∞)=10^(boost/20)`。
- Chebyshev II：取两个自然归一化 `cheb2ap` 原型（`eps1`、`eps2`），一个的极点
  作为 shelf 的极点，另一个的极点作为 shelf 的零点；两组原型的传输零点相同，
  对消后消失；再乘频率尺度 `m` 把中点放到 1 rad/s。scipy 没有偶数阶修正原型，
  所以固定 `even_modify=False`：高频端落在 `boost − sign(boost)·ripple` dB
  （棚带低谷），而棚带内的极值仍在有限频率处达到 `boost` dB。这是模拟原型
  本身的性质，不是离散化误差；要精确到棚高需要头文件里的修正原型。

fs=48 kHz、f0=1 kHz、σ=π√(2/3) + 预畸变：

| 案例 | `max(abs(ΔdB))`（20 Hz–Nyquist） |
|---|---|
| Butterworth N=8，+12 / −12 dB | `< 1e-4` |
| Chebyshev II N=8，ripple=1 dB，+12 / −12 dB | 0.062 |
| Chebyshev II N=8，ripple=1 dB，+12 / −12 dB，f0=5 kHz | 0.380 |
| Chebyshev II N=8，ripple=3 dB，+12 dB，f0=5 kHz | 1.229 |
| Chebyshev II N=16，+24 dB，ripple=3 dB，f0=5 kHz | 2.390（压力项） |

图：`output/high_shelf.png`。`check_shelf.py` 跑 42 组组合（N=4/8，±12/±24 dB，
ripple 0.5/1/3 dB，f0=1k/5k，另加 N=16 压力项），全部稳定并通过阈值。
High shelf 是最严格的一类：误差随阶数、boost、ripple、f0 增长；典型的
N=8、+12 dB、1 dB、1 kHz 仍在 0.062 dB 内。

## 观察

**Butterworth（默认配置）。** 两条曲线在全频段几乎重合；最大线性偏差出现在
过渡带 5.5 kHz 附近。24 kHz 处虚线抬高到约 −92 dB，是因为映射把数字 Nyquist
精确对应到模拟的 `fs/σ ≈ 18.7 kHz`，而不是 24 kHz：该处数字值与
`A(18.7 kHz)` 在机器精度内一致。两者在 Nyquist 的线性差只有 `2.2e-5`，
相对通带增益可以忽略。

**Chebyshev II / 虚轴零点（默认配置）。** 一般形式可以继续工作。纯虚轴零点
对应 `b = 0` 的二次因子；当零点频率低于 `fs/σ` 时映射给出 `c2 = 1`，数字零点
严格落在单位圆上，阻带零点频率与模拟值的最大偏差 0.43%（5100/6018/9039 Hz
对 5098/6013/9000 Hz）。模拟零点 25.6 kHz 高于 Nyquist，此时 `c2 < 1`，该对
零点落到实轴 `z = −0.2838`，无法再在单位圆上形成零陷——这是采样率限制，
不是映射失效。−80 dB 等波纹阻带与零点结构都保持；零点附近的 dB 尖峰来自
零点频率的微小偏移，两侧电平此时都已低于 −120 dB。

**预畸变 + σ=2。** 截止点精确，Nyquist 匹配点从 24 kHz 移到 23.67 kHz，
已经很接近 Nyquist；但整体幅度误差比默认配置大：Butterworth 的
`≥ −20 dB` 内误差 0.716 dB（默认 0.085 dB），Chebyshev II 为 0.708 dB
（默认 0.007 dB），其零点频率也偏移得更多（8760 Hz 对模拟 9000 Hz）。
它适合“必须精确截止、同时希望 Nyquist 端点尽量靠近”的场合。

**怎么选。** 只追求整体幅度精度：默认的 `σ=π√(2/3)`、无预畸变。必须精确
截止：加预畸变；其中 `σ=π√(2/3) + 预畸变` 的整体误差更小（Butterworth
`7.4e-4`、Chebyshev II `1.3e-3`），但 Nyquist 匹配点仍在 18.72 kHz。
`σ=2 + 预畸变` 是“截止精确 + Nyquist 尽可能接近”的折中，代价是通带 /
过渡带误差增大。

## 限制

- 支持实系数、因果、无 `s = 0` 极点（积分器）的有理函数；LP / HP / BP / BS
  均覆盖。`s = 0` 的零点按 (4.2) 第二行映射；`s = 0` 的极点会显式报错。
- 映射是逐因子的，高阶时误差会累积；主表是 N=8 的实测值，四类型图是 N=4
  原型的实测值。
- 只验证幅度，不涉及相位；论文指出 `σ = 0` 时相位最优，`σ` 在幅度区间时会有
  相位偏差。
- 极点位于右半平面时使用 `φ` 的倒数分支，论文提醒某个截止频率会出现除零；
  实现里对此显式报错。

## 运行

```bash
python qwqdsp/labs/ivantsovy/check_impl.py      # 公式 / 四类型 / 预畸变 / 稳定性自检，≈3 s
python qwqdsp/labs/ivantsovy/check_shelf.py     # shelf 原型 + 42 组离散化自检，≈9 s
python qwqdsp/labs/ivantsovy/exp_magnitude.py   # 主实验 4 张图，≈7 s
python qwqdsp/labs/ivantsovy/exp_four_types.py  # 两族四种类型图（2 张），≈6 s
python qwqdsp/labs/ivantsovy/exp_high_shelf.py  # high shelf 图，≈5 s
python qwqdsp/labs/ivantsovy/run_all.py         # 全部（7 张图）
```

依赖：`numpy`、`scipy`、`matplotlib`。图片写到 `output/`（已 gitignore）。
