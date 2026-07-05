
# qwqdsp

存储我曾经使用过的东西（重复造轮子！），为了防止每次都要去仓库找代码复制，将它们全部打包到一起

> [!NOTE]
> 这不是什么正经的DSP库，里面存储了大量未优化和尝试的代码

## 贡献

很乐意有人提交pull request（如果有的话就好）

## Modules

### 依赖 (Dependencies)

- **C++20** 或更高版本。
- **[Eigen3](https://eigen.tuxfamily.org)**（可选）：启用 `QWQDSP_HAVE_EIGEN`。
- **[Intel IPP](https://www.intel.com/content/www/us/en/developer/tools/oneapi/ipp.html)**（可选）：启用 `QWQDSP_HAVE_IPP`，使用 IPP 作为 FFT 后端。
- **[SIMDe](https://github.com/simd-everywhere/simde)**（可选）：启用 `QWQDSP_HAVE_SIMDE`，在非 x86 平台上模拟 SIMD 指令。
- **[raylib](https://www.raylib.com)**（可选）：构建 GUI 测试和 `playing` 可执行程序。

### 核心模块 (Core Modules)

`core`

- `convert.hpp`: 频率/单位转换工具 — `Freq2W`, `Freq2Pitch`, `Pitch2Freq`, `Freq2Mel`, `Mel2Freq`, `Samples2Decay`, `Db2Gain`, `Gain2Db`，以及模拟域频率/带宽/Q 值转换函数。
- `extension_marcos.hpp`: 编译器扩展宏 — `QWQDSP_FORCE_INLINE`, `QWQDSP_AUTO_VECTORLIZE` 等。
- `polymath.hpp`: 多精度数学工具 — 快速 `SinPi`, `SinCycle`, `CosCycle` 等三角函数近似实现。
- `interpolation4.hpp`: 备用插值聚合 — `Lagrange3rd`, `SPPCHIP` 等内联实现。
- `algebraic_waveshaper.hpp`: 代数波形塑形器 — `Naive`, `ADAA`（反锯齿，0.5 采样延迟）, `ADAA_MV`, `ADAA_MV_Compensation`。
- `adsr_envelope.hpp`: ADSR 包络发生器 — `AdsrEnvelope`（Attack/Decay/Sustain/Release 阶段控制）。

### 滤波器模块 (Filter Modules)

`filter`

- `biquad.hpp`: 双二阶滤波器（TDF2 直接 II 型转置）— `Biquad`（`Tick`, `Set`, `Reset`）。
- `biquad_coeff.hpp`: 双二阶系数结构体 — `BiquadCoeff`（`b0,b1,b2,a1,a2`）。
- `lattice_biquad.hpp`: 格型双二阶滤波器 — `LatticeBiquad`（数值稳定性优于直接型）。
- `rbj.hpp`: RBJ 滤波器设计（Audio EQ Cookbook）— `RBJ::Lowpass/Highpass/Bandpass/Peak/Lowshelf/HighShelf/Notch/Allpass`，含数字域 Q 值预畸变。
- `svf.hpp`: Chamberlin 状态变量滤波器 — `SVF`（`MakeLowpass/Highpass/Bandpass/Bell/LowShelf/HighShelf`）。
- `svf_tpt.hpp`: 梯形（TPT）状态变量滤波器 — `SVF_TPT`。
- `svf_tpt_shelf.hpp`: TPT SVF 搁架滤波器 — 低架/高架。
- `one_pole.hpp`: 一阶经典滤波器。
- `one_pole_tpt.hpp`: 一阶 TPT 滤波器。
- `onepole_tpt_shelf.hpp`: 一阶 TPT 搁架滤波器。
- `ladder.hpp`: 梯形滤波器（ladder filter）— 经典 4 极低通。
- `ladder_8pole.hpp`: 8 极梯形滤波器。
- `lattice.hpp`: 格型滤波器 — `Lattice`。
- `window_fir.hpp`: 窗函数法 FIR 设计 — `WindowFIR::Lowpass/Highpass/Bandpass/Bandstop`。
- `fir.hpp`: FIR 滤波器运行时 — `FIRDirect`（批处理）, `FIRTranspose`（转置型，逐采样）, `FirOptimise`（循环缓冲区优化）。
- `int_delay.hpp`: 整数延迟线 — `IntDelay`。
- `iir_hilbert.hpp`: IIR Hilbert 变换器（宽带 90° 移相）。
- `iir_hilbert4.hpp`: 可变阶 IIR Hilbert 变换器。
- `iir_cpx_hilbert.hpp`: 复数 IIR Hilbert 变换器。
- `any_hilbert.hpp`: 任意阶 IIR Hilbert 变换器 — `AnyHilbert`（基于并联全通结构）。
- `formant.hpp`: 共振峰滤波器 — `FormantFilter`。
- `linkwitz_riley.hpp`: Linkwitz-Riley 分频滤波器的简单实现。
- `match_biquad.hpp`: 双二阶匹配滤波器 — `MatchBiquad`。
- `median.hpp`: 中值滤波器 — `MedianFilter`。
- `ota_one_pole.hpp`: OTA 一阶滤波器 — `OTAOnePole`。
- `parallel_allpass.hpp`: 并联全通滤波器。
- `thiran_filter.hpp`: Thiran 分数延迟滤波器 — `ThiranFilter`（IIR 分数延迟近似）。
- `transpose_sallen_key.hpp`: 转置 Sallen-Key 2 极点滤波器。
- `allpass.hpp`: 全通滤波器。
- `analog_responce.hpp`: 模拟响应建模工具。
- `gold_rader.hpp`: Gold-Rader 滤波器。
- `iir_design.hpp`: IIR 滤波器设计工具。
- `iir_design_extra.hpp`: 额外 IIR 设计工具。
- `fast_set_iir_paralle.hpp`: Fast Set IIR 并联实现。
- `fixed/`: 定点数（fixed-point）滤波器实现子模块：
  - `fixed.hpp`: 聚合头。
  - `acc_traits.hpp`: 累加器类型特性。
  - `df1_biquad.hpp`: 定点 DF1 双二阶。
  - `df1_biquad_q2.hpp`: Q2 格式定点 DF1 双二阶。
  - `df1_biquad_split.hpp`: 拆分式定点 DF1 双二阶。
  - `gold_rader.hpp`: 定点 Gold-Rader 滤波器。

### 振荡器模块 (Oscillator Modules)

`oscillator`

- `blep_coeff.hpp`: BLEP 系数 — `namespace blep_coeff`，预计算的黑曼-纳托尔 BLEP 校正系数。
- `polyblep.hpp`: 多项式 BLEP 振荡器 — `PolyBlep`（锯齿波/方波/三角波/PWM）。
- `polyblep_sync.hpp`: 硬同步 BLEP 振荡器 — `PolyBlepSync`。
- `blit.hpp`: BLIT（带宽受限脉冲训练）振荡器 — `Blit`（脉冲/锯齿波/方波/三角波/正弦波）。
- `blit_pwm.hpp`: BLIT PWM 振荡器 — `BlitPWM`。
- `dsf.hpp`: 经典 DSF（离散求和公式）振荡器 — `DSFClassic`。
- `dsf_correct.hpp`: 修正 DSF 振荡器 — `DSFCorrect`。
- `noise.hpp`: 噪声发生器 — 白噪声/粉红噪声。
- `smooth_noise.hpp`: 平滑噪声发生器 — 插值平滑噪声。
- `vic_sine_osc.hpp`: Vicanek 正弦振荡器 — 双采样递归结构。
- `dr_sine_osc.hpp`:数字谐振器 正弦振荡器。
- `coridc_sine_osc.hpp`: CORDIC 正弦振荡器。
- `elliptic_sine_osc.hpp`: 椭圆函数正弦振荡器。
- `mcf_sine_osc.hpp`: MCF（Modified Couple Form/Magic Circle）正弦振荡器。
- `table_sine_v2.hpp`: 查表正弦振荡器 v2 — 波形表 + 线性插值。
- `table_sine_v3.hpp`: 查表正弦振荡器 v3 — 改进版波形表。
- `raw_oscillor.hpp`: 原始振荡器 — 基础波形生成。

### 频谱分析模块 (Spectral Modules)

`spectral`

- `real_fft.hpp`: 实数 FFT 核心接口 — `RealFFT`（`Init`, `FFT`, `IFFT`, `FFTGainPhase`, `IFFTGainPhase`, `Hilbert`, `TimeDomainShift`），自动选择 OOURA 或 IPP 后端。
- `complex_fft.hpp`: 复数 FFT — `ComplexFFT`（复数/实数输入、增益-相位变换、Hilbert），输出格式 `[0, 2π)`；`ShuffleFrequency` 静态方法。
- `oouras_real_fft.hpp`: OOURA 实数 FFT 裸封装 — `OourasRealFFT`（封装 `rdft`，输出 `[re,im]×num_bins`）。
- `oouras_complex_fft.hpp`: OOURA 复数 FFT 裸封装 — `OourasComplexFFT`（封装 `cdft`，输出 `[0, 2π)`）。
- `ipp_real_fft.hpp`: Intel IPP 实数 FFT 封装 — `IppRealFFT`（基于 `ippsFFT{Fwd,Inv}_RToCCS_32f`）。
- `ipp_complex_fft.hpp`: Intel IPP 复数 FFT 封装 — `IppComplexFFT`（基于 `ippsFFT{CToC}_{Fwd,Inv}_32f`）。
- `reassignment.hpp`: 频谱重分配（reassignment）— `Reassignment`（时频细化），`ReassignmentCorrect`（导数窗口精确重分配）。

### 音频效果模块 (FX Modules)

`fx`

- `uniform_convolution.hpp`: 均匀分块卷积 — `UniformConvolution`（基于 FFT 的分块卷积，支持任意长度脉冲响应）。
- `pitch_shifter.hpp`: 音高移位器 — `PitchShifter`（环形缓冲 + 双读取指针，2048 点）。
- `pitch_shifter2.hpp`: 音高移位器 v2。
- `plat_reverb.hpp`: Dattorro 板式混响 — `PlateReverb`（经典 Dattorro 1997 算法）。
- `limiter.hpp`: 限制器 — `SimpleLimiter`（峰值保持 + 前视延迟）。
- `resample.hpp`: 重采样器 — `Resample`。
- `resample_iir.hpp`: IIR 重采样器 — `ResampleIIR`。
- `resample_iir_dynamic.hpp`: 动态 IIR 重采样器 — 可变采样率转换。
- `delay_line.hpp`: 延迟线效果 — `DelayLine`（带反馈、混音控制）。
- `elliptic_blep.hpp`: 椭圆 BLEP 效果 — 基于椭圆滤波器的 BLEP 校正。
- `oversample.hpp`: 过采样器 — `Oversample`（多级半带滤波器过采样/降采样）。
- `polyphase_resample_fir.hpp`: 多相 FIR 重采样 — `PolyphaseDownsamplerFir` / `PolyphaseUpsamplerFir`（多相结构上采样/下采样）。

### SIMD 优化模块 (SIMD Element Modules)

`simd_element`

- `simd_pack.hpp`: 泛型 SIMD 向量模板 — `Pack4Bytes<T,N>`（`Load/Store/Broadcast`、算术运算符、`ReduceAdd`）。
- `simde_wrap4.hpp`: SIMDe SSE4.1 128 位封装 — `V4f`/`V4i`（基于 `simde__m128`）。
- `simde_wrap8.hpp`: SIMDe 256 位封装。
- `biquads.hpp`: SIMD 双二阶滤波器组 — `Biquads<N>`（N 通道并行双二阶滤波）。
- `delay_line_mono.hpp`: 单声道多延迟线 — `DelayLineMono<N,kFastTrick>`（Hermite 插值延迟读取）。
- `delay_line_single.hpp`: 多声道单延迟线 — `DelayLineSingle<N>`。
- `delay_line_stereo.hpp`: 立体声延迟线 — `DelayLineStereo<N,kFastTrick>`。
- `delay_line_multiple.hpp`: 多通道独立延迟线 — `DelayLineMultiple<N,kFastTrick>`。
- `delay_allpass.hpp`: 全通延迟 — `DelayAllpass<N>`（混响基础构建块）。
- `one_pole_tpt.hpp`: SIMD 一阶 TPT 滤波器 — `OnePoleTPT<N>`（低通/高通/高架/全通）。
- `one_pole_tpt_shelf.hpp`: SIMD 一阶 TPT 搁架滤波器 — `OnepoleTPTShelf<N>`（低架/高架/倾斜架）。
- `envelope_follower.hpp`: SIMD 包络跟随器 — `EnevelopeFollower<N>`（峰值保持 + attack/release 平滑）。
- `stereo_iir_hilbert_cpx.hpp`: SIMD 立体声 IIR Hilbert 复数变换 — `StereoIIRHilbertCpx<N>` / `StereoIIRHilbertDeeperCpx<N>`（8/16 级全通级联）。
- `plate_reverb.hpp`: SIMD 板式混响 — `PlateReverb<N>`（Dattorro 算法的 SIMD 优化实现）。
- `algebraic_waveshaper.hpp`: SIMD 代数波形塑形器 — `AlgebraicWaveshaper<N>`（`Naive`, `ADAA`, `ADAA_MV`, `ADAA_MV_Compensation`）。
- `align_allocator.hpp`: 内存对齐分配器 — `AlignedAllocator<T,Alignment>`（跨平台 `_aligned_malloc` / `aligned_alloc`）。

### 窗函数模块 (Window Modules)

`window`

- `blackman.hpp`: Blackman 窗。
- `blackman_harris.hpp`: 4 项 Blackman-Harris 窗（`kSidelobe=-92dB`）。
- `blackman_harris_3term.hpp`: 3 项 Blackman-Harris 窗（`kSidelobe=-71.48dB`，主瓣宽 3）。
- `blackman_nuttall.hpp`: Blackman-Nuttall 窗（`kSidelobe=-98.3dB`）。
- `hamming.hpp`: Hamming 窗（`kSidelobe=-43.8dB`）。
- `hann.hpp`: Hann 窗（`kSidelobe=-31.6dB`）。
- `kaiser.hpp`: Kaiser 窗（beta 参数可调），内置 `cephes::i0/i1`。
- `lanczos.hpp`: Lanczos 窗（`kSidelobe=-26.6dB`）。
- `nuttall.hpp`: Nuttall 窗（`kSidelobe=-93.3dB`）。
- `taylor.hpp`: Taylor 窗（旁瓣电平 + nbars 参数）。
- `helper.hpp`: 窗工具函数（归一化、时间加权、零填充等）。

### 插值模块 (Interpolation Modules)

`interpolation`

- `linear.hpp`: 线性插值。
- `catmull_rom_spline.hpp`: Catmull-Rom 样条插值。
- `makima.hpp`: 修正 Akima 插值。
- `sppchip.hpp`: 单调保形三次 Hermite 插值（SPPCHIP）。

### 分段处理模块 (Segmentation Modules)

`segement`

- `analyze.hpp`: 分析基类。
- `analyze_auto.hpp`: 自动分段分析 — `AnalyzeAuto`（帧处理循环，回调模式）。
- `analyze_synthsis_offline.hpp`: 离线分析综合。
- `analyze_synthsis_online.hpp`: 在线（实时）分析综合。
- `mono_reader.hpp`: 单声道文件读取器 — `MonoReader`。
- `slice.hpp`: 切片工具 — `Slice1D`（一维数据分块迭代）。

### 自适应滤波模块 (Adaptive Filter Modules)

`adaptive`

- `nlms.hpp`: 归一化最小均方（NLMS）自适应滤波器。
- `rls_filter.hpp`: 递归最小二乘（RLS）自适应滤波器。
- `burg_lp.hpp`: Burg 法线性预测（LP）系数估计。
- `lag_buffer.hpp`: 延迟缓冲区 — 自适应滤波器所需的信号延迟管理。

### 基音检测模块 (Pitch Detection Modules)

`pitch`

- `pitch.hpp`: 基音数据结构 — `Pitch`（`pitch_hz`, `non_period_ratio`）。
- `yin.hpp`: YIN 基音检测算法。
- `fast_yin.hpp`: 快速 YIN 基音检测算法。
- `mpm.hpp`: McLeod Pitch Method（MPM）基音检测。
- `helmholtz.hpp`: Helmholtz 基音检测。

### Cephes 数学函数模块 (Cephes Math Functions)

`cephes`

- `bessel.hpp`: 修正 Bessel 函数 I₀(x)、I₁(x)（双精度），切比雪夫级数逼近。
- `besself.hpp`: 修正 Bessel 函数 I₀(x)、I₁(x)（单精度）。
- `elliptic.hpp`: 完全 / 不完全椭圆积分（`ellpk`, `ellik`, `ellpj`）。
- `ellipticf.hpp`: 单精度椭圆积分（`ellpkf`, `ellikf`, `ellpjf`）。

### 杂项工具模块 (Miscellaneous Modules)

`misc`

- `smoother.hpp`: 指数平滑器 — `ExpSmoother`（attack/release 可调）。
- `envelope_follower.hpp`: 包络跟随器。
- `peakfind.hpp`: 峰值查找 — `PeakFind`。
- `ampd_peak.hpp`: AMPD（Automatic Multi-scale Peak Detection）峰值检测算法。
- `crossover.hpp`: 分频器 — `Crossover`。
- `integrator.hpp`: 积分器 — `Integrator`（泄漏积分器）。

## Tests

### 滤波器测试 (Filter Tests)

`tests/nogui/filter/`

- `biquad.cpp`: 比较 `Biquad` 与 `LatticeBiquad` 的脉冲响应（RBJ 低通系数）。
- `fir_design.cpp`: `WindowFIR::Bandstop` 带阻 FIR 设计 + Hamming 窗 + FFT 频响验证。
- `minimum_phase_fir.cpp`: 离散 Hilbert 变换（倒谱法）将线性相位 FIR 转为最小相位。
- `paralle_allpass.cpp`: 原型/设计文档 — 被注释掉的 SIMD 并行全通，`main()` 为空。
- `reverse_iir.cpp`: Vicanek 反向时间 IIR 线性相位补偿。
- `residual_chebyshev.cpp`: Chebyshev Type II 设计（s 域极点 → 双线性变换 → 部分分式 → 双二阶）。

### 音频效果测试 (Audio FX Tests)

`tests/nogui/audiofx/`

- `auto_notch.cpp`: FIR 房间脉冲响应 + FFT 啸叫检测 + 自适应陷波。
- `wiener_filter.cpp`: 4 场景维纳滤波去噪（噪声抑制、声码器等）。
- `convolution.cpp`: `UniformConvolution` 分块卷积正确性。
- `conv_feedback.cpp`: 均匀卷积 + 反馈环路的梳状/延迟/混响效果。
- `conv_feedback2.cpp`: 反馈卷积变体。
- `reverb.cpp`: C 风格 FDN 混响（全通 + 16 路 FDN + 架桥滤波 + 合唱）。
- `freq_shifter.cpp`: IIR Hilbert 复单边带调制，偏移 -150 Hz。
- `resample.cpp`: `ResampleIIR` 扫频重采样至 96 kHz。
- `polyphase.cpp`: 原型/设计文档 — `#if 0` 注释掉的多相抽取器。

### 频谱处理测试 (Spectral Tests)

`tests/nogui/spectral/`

- `fft_interpolation.cpp`: FFT 域频谱缩放实现时域插值（上采样）。
- `ifft_random_phase.cpp`: 随机相位谱 + IFFT 生成全通 FIR。

### 合成测试 (Synth Tests)

`tests/nogui/synth/`

- `dsf.cpp`: `DSFCorrect` 复现 BLIT 论文波形。
- `stupid_resynthsis.cpp`: 频谱重分配 + 峰值检测 + 查表正弦波重合成。
- `wsola.cpp`: WSOLA 时间伸缩（Hann 窗 + 相似性搜索，2.0 倍拉伸）。

### Cephes 数学验证 (Cephes Math Verification)

`tests/nogui/cephes/`

- `cephes.cpp`: Bessel 函数、椭圆积分正确性验证。
- `cephes_ref.inc`: 参考值数据。
- `check_cephes.py`: Python 精度对比。

### GUI 实时测试 (Realtime GUI Tests)

`tests/gui/`

- `interpolations.cpp`: GUI 实时 (raylib) — 拖拽控制点对比 SPPCHIP / MAKIMA / Catmull-Rom / Linear 插值曲线。
- `auto_notch_rt.cpp`: GUI 实时 (raylib + miniaudio) — 实时啸叫检测与自适应陷波。

`tests/gui/spectral/`

- `polyphase_filter_bank_view.cpp`: GUI 实时 (raylib + miniaudio) — 多相分析滤波器组（DFT 调制，M=512）。
- `polyphase_filter_bank_view2.cpp`: GUI 实时 (raylib + miniaudio) — 余弦调制分析滤波器组（CMFB，M=256，DCT-IV）。
- `resonate_bank.cpp`: GUI 实时 (raylib + SIMD) — AVX Float256 加速谐振器组（最多 1024）。
- `spectrum3.cpp`: GUI 实时 (raylib + miniaudio) — 1/12 倍频程平滑频谱分析仪（4096 点 FFT）。
- `spectrum/`: GUI 实时 (raylib + miniaudio) — 25 bin 频谱分析仪（LED 条形图，嵌入式 C 移植）。
- `spectrum2/`: GUI 实时 (raylib + miniaudio) — 81 bin + 32 频带 Mel 滤波器组。

`tests/gui/synth/`

- `formant_filter.cpp`: GUI 实时 (raylib) — 5 级联格型双二阶 + DSF 共振峰合成器。
- `blit.cpp`: GUI 实时 (raylib) — BLIT 振荡器 8 波形实时演示。
- `polyblep.cpp`: GUI 实时 (raylib) — PolyBLEP 振荡器 9 波形（含硬同步）演示。

### 轻量级 RNNoise 噪声抑制 (RNNoise Lite)

`tests/rnnoise_lite/`

- `main.cpp`: WAV 加载 → RNNoise 降噪 → 保存输出。
- `rnnoise.c` / `rnnoise.h`: C 推理引擎（20 频带双二阶 + MFCC + 3 层 GRU + nnom int8）。
- `rnn.hpp`: C++ 推理封装（nnom + MFCC，50% 重叠帧）。
- `mfcc.c` / `mfcc.h`: MFCC 特征提取（FFT + 三角梅尔滤波器组 + DCT）。
- `weights.h`: 3 层 GRU int8 量化权重。
- `denoise_weights.h`: 对应浮点 GRU 权重。
- `equalizer_coeff.h`: 20 频带带通双二阶系数。
- `nnom/`: NNoM 轻量级神经网络运行时（int8 量化，GRU / 全连接 / 张量运算）。

## Notebooks

### Python Jupyter 笔记本

`notebooks/`

- `adaa.ipynb`: ADAA（反导数抗混叠）波形塑形验证 — 使用 SymPy 推导代数波形塑形器的积分和 ADAA 公式（Naive / ADAA-2 / ADAA-MV），并用 Numpy + FFT 频谱分析对比各版本的混叠抑制效果。
- `blep.ipynb`: BLEP 多项式系数优化 — 使用 SymPy 定义分段 BLEP 校正多项式，通过 SciPy 约束优化（端点条件 + 单调性约束）求解 MSE 最优的 4 阶多项式系数。
- `cheby_hilbert.ipynb`: 切比雪夫 II 型 / 椭圆滤波器并联 IIR 设计 — 设计模拟原型滤波器，经双线性变换和部分分式分解后导出 C++ 残差-极点系数代码。
- `fastset_iir.ipynb`: Fast Set IIR 滤波器系数设计 — 通过极点配置优化不同建立时间常数（1e2~1e7 采样）的 2/4/6/8 阶 IIR 滤波器，导出并联形式的 C++ 结构体代码。
- `holters_parker_coeff.ipynb`: Holters-Parker 滤波器系数设计 — 设计椭圆 / 切比雪夫 I 型模拟原型，通过部分分式展开为并联 IIR，生成 C++ 系数代码。
- `tpt_filter.ipynb`: TPT（梯形）滤波器拓扑符号推导 — 使用 SymPy 推导 SVF、转置 Sallen-Key 2 极点、4 极梯形低通 / 高通滤波器的差分方程。
- `ola.ipynb`: 时变 FFT 卷积与 OLA/WOLA 重建实验 — ① 三种步长（R=M, R=M/2, R=M/4）下时变核（FIR 低通→高通突变）的 FFT 卷积重建质量对比；② WOLA 中 sqrt-Hann 窗在 50% 重叠下的完美重建演示；③ Hann 与 sqrt-Hann 窗的频谱对比；④ 75% 重叠下 OLA（Hann）与 WOLA（Hann²）的重建增益常数分析；⑤ 多种窗函数（Rect / Hann / Hamming / Blackman / Blackman-Harris）在其最小完美重建步长下的 OLA 稳态纹波精度验证。
