# 仓库实现（`resample_iir.hpp` + `elliptic_blep.hpp`）基于论文的哪些工作

## 结论

**只用了论文 §3.1「输入滤波器」那一半 + §3 的 modified impulse-invariant transform。**
另一半（§3.2 输出滤波器的阶跃/矩形保持形式）没有用；BBD 延迟线本身、Juno-60 电路系数也没用。

| 论文 | 仓库 |
|---|---|
| §3 引子：不做额外插值，直接用电路自带低通完成重采样 | 整个思路 |
| 式 (2) `H(s)=Σ r_m/(s−p_m)`（部分分式） | `resample_coeffs.h` 的 `complexPoles/realPoles` + `complexCoeffsDirect/realCoeffsDirect`（由 notebook 对模拟椭圆滤波器做 `scipy.signal.residue` 得到） |
| 式 (3) `h(t)=Σ r_m e^{p_m t}` | `partial_step_poles_[s][i] = exp(s/N · pole)`、`impluse_coeffs_[i]` |
| 式 (4)(5) 把时间拆成 `t=(l_n+d_n)T_s` | `ResampleIIR::Process` 的 `rpos`(整数) + `phase`(分数) |
| 式 (6) `g_in,m(d_n)=T_s·r_m·p̄^{d_n}` | `Add()` 上 `r_m·scale`、`Get(d)` 上 `p̄^d`（两者相乘） |
| 式 (8)(9) 递归 `x_m(l)=p̄_m x_m(l−1)+ū(l)` | `Step()`：`state_[i] *= partial_step_poles_.back()[i]`；`Add(amount)` 注入 |
| 式 (7)(10) `u((l_n+d_n)T_s)=Σ_m g_in,m(d_n)·x_in,m(l_n)` | `Get(samplesInFuture)`：`Σ Re(state_[i]·lerp(p̄_i^{d}))` |
| §3.3 "…may use polynomial approximations or **look-up tables** for the b coefficients" | `kPartialSteps` + 线性插值，正是论文允许的实现选择（论文说误差分析超出其范围，本目录把它量了） |

**没用到**：§3.2 输出滤波器（式 12–25）、`H_0`、BBD 的变速率延迟线、Table 1 的 Juno-60 系数、
§4 的过采样建议、§3.3 的实系数二阶合并（仓库直接留复数极点，不做合并）。

## 相对论文的三处改动（就是"更强"的来源）

1. **滤波器从「电路给定的」变成「可设计的」**：论文的 `H_in` 是电路固定截止；这里用任意椭圆低通
   （`BestCoeffs` 25 阶 / `MedianCoeffs` 31 阶 / `FastCoeffs` 11 阶），
   截止由 `SetCutoffByFpass(min(target_fs, source_fs)/2, source_fs)` 摆到 `min(fs)/2`。
   （NonIntegerSRC 是摆在 `fs/2`，降采样不抗混叠；这里是对的。）
2. **`p̄^{d}` 用 LUT 代替逐样本 `exp`**：论文只提了可以这么做。
3. **同一套冲激不变核被复用到 BLEP**：`Add(amount, samplesInPast)` 与 coeffs 里的
   "direct bandlimited synthesis of a polynomial-segment waveform" 注释说明这是给
   bandlimited 振荡器用的，和 SRC 共用 `partial_step_poles_`。

## 实测（`temp/hp_probe.cpp` / `hp_analyze.py`，clang++ -O2，与 CMake 同工具链）

方法与 §2 一致：C++ 探针只 dump 样本，Python 侧做最小二乘拟合 + 逐点比对设计方程。
期望增益 = `|H_old(j·2f·fpass/min(fs_in,fs_out))|`（`H_old` 由 coeff 表按式 (2) 解析求值）。

### 1. 实测增益 vs 设计方程（`BestCoeffs<float>`，`kPartialSteps=128`）

| 用例 | f (Hz) | 实测 amp | 期望 \|H\| | 残差 rms |
|---|---|---|---|---|
| 48k→96k | 1000 | 0.9955 | 0.9955 | 0.0000 |
| 48k→96k | 10000 | 0.9936 | 0.9935 | 0.0002 |
| 48k→96k | 20000 | 0.9993 | 0.9993 | 0.0005 |
| 48k→96k | 23000 | 0.9948 | 0.9947 | 0.0015 |
| 96k→48k | 20000 | 0.9993 | 0.9993 | 0.0001 |
| 96k→48k | 26000 / 30000 / 40000 | — | — | **0.0002**（音在输出 Nyquist 之上，被滤干净） |
| 48k→44.1k | 20000 | 0.9941 | 0.9941 | 0.0004 |
| 48k→44.1k | 23000 | — | — | 0.0003（`min(fs)/2=22.05k` 之上，滤掉） |

结论：**幅频逐点吻合到 ~1e-4，降采样的抗混叠真的生效**（对照 NonIntegerSRC：30 kHz 折回 18 kHz 幅度 ≈1.0057）。
残差随 f 靠近截止上升，是滤波器的过渡带/镜频残余，和滤波器阶数一致
（`FastCoeffs` 11 阶在 23 kHz 处残差 0.1063，`Best` 25 阶只有 0.0015）。

### 2. 与"精确连续核直接求和"的独立复算（double 精度）

把 `y(t)=Σ_k x[k]·h_a(t−k)`、`h_a(τ)=Σ_m Re(c_m e^{λ_m τ})` 用 numpy 直接算（不用 LUT、不递归），
与 C++ 输出逐点比：

| 实例化 | max\|e\| | rms |
|---|---|---|
| `BestCoeffs<double>` | **6.4e-15** | 3.0e-15 |
| `BestCoeffs<float>` | 1.97e-03 | 4.5e-04 |

float 版的误差**随 n 线性增长**、且拟合出来是正交（相位）分量、增益分量仅 1e-6
→ 是 float32 每步极点相位舍入的相干累积（`Sample=float` 的固有代价，不是算法错）。
想要更好就把 `TSample` 换成 `double`（`EllipticBlep` 全部按 `Sample` 模板化，改模板参数即可）。

### 3. `kPartialSteps` 的实际代价（double 精度，48k→44.1k，20 kHz，只有 LUT 是误差源）

| kPartialSteps | max\|e\| | rms |
|---|---|---|
| 8 | 1.42e-02 | 6.85e-03 |
| 32 | 9.14e-04 | 4.29e-04 |
| 128 | 5.40e-05 | 2.69e-05 |
| 512 | 3.25e-06 | 1.68e-06 |

二阶收敛（每 4 倍 → 误差 /16.6），与线性插值 `e^{dλ}` 的 `~|λ|²/(8N²)` 一致。
**128 档 = -85 dB，远低于 float32 自己的 ~-54 dB**，所以现在这个 128 是够的、不是瓶颈。

### 4. 单位脉冲 → 输出就是连续核在新栅格上的采样（`temp/hp_ir_probe.cpp`）

`x[0]=1`、其余为 0 时，输出应当恰好是 `h_a(n·fs_in/fs_out)`（`h_a(τ)=Σ_m Re(c_m·scale·e^{λ_m τ})`，τ 以输入样本计）：

| 用例 | 分数偏移 `d_n` | max\|C++ − 解析 h_a\| |
|---|---|---|
| 48k→96k | 恒为 {0, 0.5}（正好是 LUT 节点） | **7.2e-16** |
| 48k→44.1k | 0.0884…（落在节点之间，要插值） | 6.5e-06（换成精确 `exp` 后 1e-16 ⇒ 全部来自 LUT） |

即"输出 = 连续冲激响应在新时间栅格上的采样"是字面成立的，映射平移不变
（输入移 1 个样本 → 输出移 `fs_in/fs_out` 个样本）。第二行多出的 6.5e-6 **全部来自 LUT 插值**
——比率非整数时 `d_n` 落在 LUT 节点之间，整数比 1/2 时 `d_n∈{0,0.5}` 正好是节点，所以误差为零。

## 顺带发现（不是本次调查的目标，但记录）

`resample_iir_dynamic.hpp:44` 调用了 `blep_.SetCutoff(cutoff, source_fs_)`，
而 `EllipticBlep` 只有 `SetCutoffByFpass/SetCutoffByFstop/SetScale` —— **没有 `SetCutoff`**。
模板成员只有实例化时才查名，所以头文件被 `fx/fx.hpp` 包含也不报错；
一旦实例化 `ResampleIIRDynamic<...>`（我试了 `BestCoeffs<float>,128`）就编译失败：

```
error: no member named 'SetCutoff' in 'signalsmith::blep::EllipticBlep<...>'
```

即它是**没被实例化过的死代码**（全仓库只有 `fx.hpp` 的 `#include`，无使用点）。
`SetRatio` 里"按 cut/ratio"的意图应写成 `SetCutoffByFpass(min(source_fs, source_fs/ratio)/2, source_fs)`
或直接 `SetScale`。

## 相关

- 原型的**偶数阶修正**（`even_order_modify`）与**直接项**（`d ≠ 0` 的原型怎么处理、怎么实现、
  换截止时怎么缩放）都在 [algorithm.md](algorithm.md) §9；四步数值研究见
  [direct_term_study.py](direct_term_study.py)，图在 `output/`。
- 索引与结论速览见 [README.md](README.md)。
