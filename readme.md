
# qwqdsp

存储我曾经使用过的东西（重复造轮子！），为了防止每次都要去仓库找代码复制，将它们全部打包到一起

> [!NOTE]
> 这不是什么正经的DSP库，里面存储了大量未优化和尝试的代码

> 贡献  
> 很乐意有人提交pull request（如果有的话就好）  

---

## 📦 依赖 (Dependencies)

- **C++20** 或更高版本。
- **[Eigen3](https://eigen.tuxfamily.org)**（可选）：启用 `QWQDSP_HAVE_EIGEN`。
- **[Intel IPP](https://www.intel.com/content/www/us/en/developer/tools/oneapi/ipp.html)**（可选）：启用 `QWQDSP_HAVE_IPP`，使用 IPP 作为 FFT 后端。
- **[SIMDe](https://github.com/simd-everywhere/simde)**（可选）：启用 `QWQDSP_HAVE_SIMDE`，在非 x86 平台上模拟 SIMD 指令。
- **[raylib](https://www.raylib.com)**（可选）：构建 GUI 测试和 `playing` 可执行程序。

---

## 🎮 demos

### 🎛️ 滤波器示例 (Filter Demos)

`tests/nogui/filter/`

| 文件 | 说明 |
|------|------|
| [`biquad.cpp`](tests/nogui/filter/biquad.cpp) | 对比 Biquad 与 LatticeBiquad 的脉冲响应差异 |
| [`fir_design.cpp`](tests/nogui/filter/fir_design.cpp) | WindowFIR 带阻滤波器设计 + Hamming 窗 + FFT 频响验证 |
| [`minimum_phase_fir.cpp`](tests/nogui/filter/minimum_phase_fir.cpp) | 离散 Hilbert 变换将线性相位 FIR 转为最小相位 |
| [`paralle_allpass.cpp`](tests/nogui/filter/paralle_allpass.cpp) | 被注释掉的并行全通 |
| [`reverse_iir.cpp`](tests/nogui/filter/reverse_iir.cpp) | Vicanek 反向时间 IIR 实现线性相位补偿 |
| [`residual_chebyshev.cpp`](tests/nogui/filter/residual_chebyshev.cpp) | Chebyshev II 型从 s 域极点推导到并行双二阶滤波器 |
| [`any_hilbert.cpp`](tests/nogui/filter/any_hilbert.cpp) | AnyHilbert 脉冲响应输出 |
| [`plot_hilbert_ir.py`](tests/nogui/filter/plot_hilbert_ir.py) | AnyHilbert 脉冲响应绘图脚本 |

### 🔊 音频效果示例 (Audio FX Demos)

`tests/nogui/audiofx/`

| 文件 | 说明 |
|------|------|
| [`auto_notch.cpp`](tests/nogui/audiofx/auto_notch.cpp) | FIR 房间脉冲响应 + FFT 啸叫检测 + 自适应陷波抑制 |
| [`wiener_filter.cpp`](tests/nogui/audiofx/wiener_filter.cpp) |维纳滤波去噪和vocoder(?) |
| [`convolution.cpp`](tests/nogui/audiofx/convolution.cpp) | UniformConvolution 分块卷积正确性验证 |
| [`conv_feedback.cpp`](tests/nogui/audiofx/conv_feedback.cpp) | 均匀卷积 + 反馈环路构成的梳状/混响效果 |
| [`conv_feedback2.cpp`](tests/nogui/audiofx/conv_feedback2.cpp) | 反馈卷积的另一种变体实现 |
| [`freq_shifter.cpp`](tests/nogui/audiofx/freq_shifter.cpp) | 频率移动，有无抗混叠对比 |
| [`resample.cpp`](tests/nogui/audiofx/resample.cpp) | 扫频重采样至 96 kHz，测试FIR/IIR/多相 重采样 |
| [`reverb.cpp`](tests/nogui/audiofx/reverb.cpp) | C 风格 FDN 混响 |
| [`oversample.cpp`](tests/nogui/audiofx/oversample.cpp) | 过采样对代数波形塑形器的失真抑制效果对比 |

### 📊 频谱处理示例 (Spectral Demos)

`tests/nogui/spectral/`

| 文件 | 说明 |
|------|------|
| [`fft_interpolation.cpp`](tests/nogui/spectral/fft_interpolation.cpp) | FFT 域频谱缩放实现时域上采样 |
| [`ifft_random_phase.cpp`](tests/nogui/spectral/ifft_random_phase.cpp) | 随机相位谱 + IFFT 生成全通 FIR |
| [`fold.cpp`](tests/nogui/spectral/fold.cpp) | FFT的时间混叠增强通道不相干性 |

### 🎹 合成示例 (Synth Demos)

`tests/nogui/synth/`

| 文件 | 说明 |
|------|------|
| [`dsf.cpp`](tests/nogui/synth/dsf.cpp) | DSFCorrect 复现 BLIT 论文波形 |
| [`stupid_resynthsis.cpp`](tests/nogui/synth/stupid_resynthsis.cpp) | 简单的重新合成测试 |
| [`wsola.cpp`](tests/nogui/synth/wsola.cpp) | WSOLA |

### ✅ Cephes 数学验证 (Cephes Math Verification)

`tests/nogui/cephes/`

| 文件 | 说明 |
|------|------|
| [`cephes.cpp`](tests/nogui/cephes/cephes.cpp) | Bessel 函数和椭圆积分正确性验证 |
| [`cephes_ref.inc`](tests/nogui/cephes/cephes_ref.inc) | 参考值数据 |
| [`check_cephes.py`](tests/nogui/cephes/check_cephes.py) | Python 精度对比脚本 |

### 🖥️ GUI 实时演示 (Realtime GUI Demos)

`tests/gui/`

| 文件 | 说明 |
|------|------|
| [`interpolations.cpp`](tests/gui/interpolations.cpp) | 拖拽控制点实时对比四种插值曲线 |
| [`auto_notch_rt.cpp`](tests/gui/auto_notch_rt.cpp) | 实时啸叫检测与自适应陷波 |
| [`filters.cpp`](tests/gui/filters.cpp) | 多种滤波器（SVF / Ladder / Sallen-Key / OTA） |
| [`filters2.cpp`](tests/gui/filters2.cpp) | 过采样 + 非线性滤波器 |

### 📈 频谱分析 (Spectral)
`tests/gui/spectral/`

| 文件 | 说明 |
|------|------|
| [`polyphase_filter_bank_view.cpp`](tests/gui/spectral/polyphase_filter_bank_view.cpp) | DFT 调制多相分析滤波器组（M=512） |
| [`polyphase_filter_bank_view2.cpp`](tests/gui/spectral/polyphase_filter_bank_view2.cpp) | 余弦调制分析滤波器组（M=256，DCT-IV） |
| [`resonate_bank.cpp`](tests/gui/spectral/resonate_bank.cpp) | AVX Float256 加速谐振器组（最多 1024） |
| [`spectrum3.cpp`](tests/gui/spectral/spectrum3.cpp) | 1/12 倍频程平滑频谱分析仪（4096 点 FFT） |
| [`spectrum/`](tests/gui/spectral/spectrum/) | LED 条形图频谱分析仪（25 bin，嵌入式 C 移植） |
| [`spectrum2/`](tests/gui/spectral/spectrum2/) | Mel 滤波器组频谱分析仪（81 bin + 32 频带） |
| [`spectrum4.cpp`](tests/gui/spectral/spectrum4.cpp) | 粗略的多分辨率频谱分析仪 |
| [`reassignment.cpp`](tests/gui/spectral/reassignment.cpp) | 实时频谱重分配时频图 |

### 🎛️ 合成器 (Synth)
`tests/gui/synth/`

| 文件 | 说明 |
|------|------|
| [`formant_filter.cpp`](tests/gui/synth/formant_filter.cpp) | DSF 共振峰合成器 |
| [`blit.cpp`](tests/gui/synth/blit.cpp) | BLIT 振荡器 8 波形实时演示 |
| [`polyblep.cpp`](tests/gui/synth/polyblep.cpp) | PolyBLEP 振荡器 9 波形（含硬同步）演示 |
| [`noise.cpp`](tests/gui/synth/noise.cpp) | 噪声（白/粉红/布朗/点击）演示 |

### 🧠 轻量级 RNNoise 噪声抑制 (RNNoise Lite)

`tests/rnnoise_lite/`

| 文件 | 说明 |
|------|------|
| [`main.cpp`](tests/rnnoise_lite/main.cpp) | WAV 加载 → RNNoise 降噪 → 保存输出的完整流程 |

---

### 📓 Notebooks

`notebooks/`

| 文件 | 说明 |
|------|------|
| [`adaa.ipynb`](notebooks/adaa.ipynb) | SymPy 推导 ADAA 波形塑形器积分公式 + FFT 混叠抑制对比 |
| [`any_hilbert.ipynb`](notebooks/any_hilbert.ipynb) | 任意 IIR 半带滤波器留数分解为 Hilbert 变换器并联形式 |
| [`blep.ipynb`](notebooks/blep.ipynb) | BLEP相关的东西 |
| [`cheap_elliptic.ipynb`](notebooks/cheap_elliptic.ipynb) | 廉价椭圆滤波器设计，含 SymPy 推导和交互式相频响应 |
| [`colormap.ipynb`](notebooks/colormap.ipynb) | matplotlib colormap |
| [`fastset_iir.ipynb`](notebooks/fastset_iir.ipynb) | Fast Set IIR 极点配置系数优化，导出 C++ 结构体 |
| [`halfband.ipynb`](notebooks/halfband.ipynb) | 半带滤波器的设计与分析 |
| [`holters_parker_coeff.ipynb`](notebooks/holters_parker_coeff.ipynb) | Holters-Parker 重采样技术的 系数生成 |
| [`iir_compare.ipynb`](notebooks/iir_compare.ipynb) | IIR 原型（Butterworth/Chebyshev I-II/Elliptic/Bessel）幅频对比 |
| [`ivantsovy.ipynb`](notebooks/ivantsovy.ipynb) | Yuriy Ivantsov 理想双线性/双二次数字滤波器理论验证 |
| [`ivantsovy_1st_svf.ipynb`](notebooks/ivantsovy_1st_svf.ipynb) | Ivantsov 一阶 SVF 直接计算方法推导 |
| [`ivantsovy_2nd_svf.ipynb`](notebooks/ivantsovy_2nd_svf.ipynb) | Ivantsov 二阶 SVF 多模式输出传递函数推导(**损坏**) |
| [`ivantsovy_state_space.ipynb`](notebooks/ivantsovy_state_space.ipynb) | Ivantsov 状态空间滤波器与其他滤波器架构(TDF2/SVF/LatticeLadder)对比 |
| [`ola.ipynb`](notebooks/ola.ipynb) | 时变 FFT 卷积与 OLA/WOLA 完美重建实验 |
| [`ola_long_convolution.ipynb`](notebooks/ola_long_convolution.ipynb) | 均匀分块 OLA 长脉冲响应卷积 |
| [`reassignment.ipynb`](notebooks/reassignment.ipynb) | 时间-频率重分配方法实现与分析 |
| [`sst.ipynb`](notebooks/sst.ipynb) | ssqueezepy 库的 STFT/CWT 时频分析演示 |
| [`tpt_filter.ipynb`](notebooks/tpt_filter.ipynb) | SymPy 推导 TPT 滤波器（SVF/Sallen-Key/梯形）差分方程 |

---

### 🏗️ 核心模块 (Core Modules)

| 文件 | 说明 |
|------|------|
| [`convert.hpp`](include/qwqdsp/convert.hpp) | 频率/单位转换工具 |
| [`extension_marcos.hpp`](include/qwqdsp/extension_marcos.hpp) | 编译器扩展宏 |
| [`polymath.hpp`](include/qwqdsp/polymath.hpp) | 快速三角函数近似实现 |
| [`interpolation4.hpp`](include/qwqdsp/interpolation4.hpp) | 备用插值算法 — `Interpolation4` |
| [`algebraic_waveshaper.hpp`](include/qwqdsp/algebraic_waveshaper.hpp) | 代数波形塑形器 — `AlgebraicWaveshaper` |
| [`adsr_envelope.hpp`](include/qwqdsp/adsr_envelope.hpp) | ADSR 包络发生器 — `AdsrEnvelope` |

### 🎛️ 滤波器模块 (Filter Modules)

| 文件 | 说明 |
|------|------|
| [`allpass.hpp`](include/qwqdsp/filter/allpass.hpp) | `AllpassOrder1` 一阶全通滤波器<br>`AllpassOrder2` 二阶全通滤波器<br>`AllpassPolyphase` Schroeder 多相全通<br>`AllpassOrder1Polyphase` 一阶可调延迟全通<br>`AllpassOrder2Polyphase` 二阶可调延迟全通 |
| [`analog_responce.hpp`](include/qwqdsp/filter/analog_responce.hpp) | `AnalogResponce` 模拟频响计算工具 |
| [`any_hilbert.hpp`](include/qwqdsp/filter/any_hilbert.hpp) | `AnyHilbert` 任意并联 Hilbert 变换器 |
| [`biquad.hpp`](include/qwqdsp/filter/biquad.hpp) | `Biquad` TDF2 双二阶滤波器 |
| [`biquad_coeff.hpp`](include/qwqdsp/filter/biquad_coeff.hpp) | `BiquadCoeff` 单精度双二阶系数<br>`DoubleBiquadCoeff` 双精度双二阶系数 |
| [`fast_set_iir_paralle.hpp`](include/qwqdsp/filter/fast_set_iir_paralle.hpp) | `FastSetIirParalle` fast set iir实现 |
| [`fir.hpp`](include/qwqdsp/filter/fir.hpp) | `FIRDirect` 批处理 FIR<br>`FIRTranspose` 逐采样转置 FIR<br>`FirOptimise` 循环缓冲区优化 FIR |
| [`formant.hpp`](include/qwqdsp/filter/formant.hpp) | `Formant` 元音共振峰查询结构体 |
| [`gold_rader.hpp`](include/qwqdsp/filter/gold_rader.hpp) | `GoldRader` Gold-Rader 滤波器 |
| [`iir_cpx_hilbert.hpp`](include/qwqdsp/filter/iir_cpx_hilbert.hpp) | `IIRHilbertCpx` 8 级复数 IIR Hilbert<br>`IIRHilbertDeeperCpx` 16 级复数 IIR Hilbert |
| [`iir_design.hpp`](include/qwqdsp/filter/iir_design.hpp) | `IIRDesign` IIR 原型滤波器设计工具 |
| [`iir_design_extra.hpp`](include/qwqdsp/filter/iir_design_extra.hpp) | `IIRDesignExtra` 额外 IIR 设计方法 |
| [`iir_hilbert.hpp`](include/qwqdsp/filter/iir_hilbert.hpp) | `IIRHilbert` 8 级实数 IIR Hilbert<br>`IIRHilbertDeeper` 16 级实数 IIR Hilbert |
| [`iir_hilbert4.hpp`](include/qwqdsp/filter/iir_hilbert4.hpp) | `IIRHilbertFull` 可变阶 IIR Hilbert 变换器 |
| [`int_delay.hpp`](include/qwqdsp/filter/int_delay.hpp) | `IntDelay` 整数采样延迟线 |
| [`ladder.hpp`](include/qwqdsp/filter/ladder.hpp) | `Ladder` 4 极梯形滤波器 |
| [`ladder_8pole.hpp`](include/qwqdsp/filter/ladder_8pole.hpp) | `Ladder8Pole` 8 极梯形滤波器 |
| [`lattice.hpp`](include/qwqdsp/filter/lattice.hpp) | `LatticeZero` 零点格型节<br>`LatticeZeroPolyphase` 多相零点格型节<br>`LatticePole` 极点格型节<br> `LatticePolePolyphase` 多相极点格型节 |
| [`lattice_biquad.hpp`](include/qwqdsp/filter/lattice_biquad.hpp) | `LatticeBiquad` 格型双二阶滤波器 |
| [`linkwitz_riley.hpp`](include/qwqdsp/filter/linkwitz_riley.hpp) | `LinkwitzRiley2` 二阶 Linkwitz-Riley<br>`LinkwitzRiley4` 四阶 Linkwitz-Riley<br>`LinkwitzRiley6` 六阶 Linkwitz-Riley<br>`LinkwitzRiley8` 八阶 Linkwitz-Riley |
| [`match_biquad.hpp`](include/qwqdsp/filter/match_biquad.hpp) | `MatchBiquad` 双二阶匹配设计工具 |
| [`median.hpp`](include/qwqdsp/filter/median.hpp) | `MedianDynamic` 动态大小中值滤波器<br>`Median` 编译期大小中值滤波器 |
| [`one_pole.hpp`](include/qwqdsp/filter/one_pole.hpp) | `OnePoleFilter` 一阶经典滤波器 |
| [`one_pole_tpt.hpp`](include/qwqdsp/filter/one_pole_tpt.hpp) | `OnePoleTPT` 一阶 TPT 滤波器 |
| [`onepole_tpt_shelf.hpp`](include/qwqdsp/filter/onepole_tpt_shelf.hpp) | `OnepoleTPTShelf` 一阶 TPT 搁架滤波器 |
| [`ota_one_pole.hpp`](include/qwqdsp/filter/ota_one_pole.hpp) | `OTAOnePole` OTA 一阶滤波器 |
| [`parallel_allpass.hpp`](include/qwqdsp/filter/parallel_allpass.hpp) | `ParallelAllpass` 双路并联全通滤波器 |
| [`rbj.hpp`](include/qwqdsp/filter/rbj.hpp) | `RBJ` Audio EQ Cookbook 系数设计 |
| [`svf.hpp`](include/qwqdsp/filter/svf.hpp) | `SVF` A.Simpler 状态变量滤波器 |
| [`svf_tpt.hpp`](include/qwqdsp/filter/svf_tpt.hpp) | `SvfTPT` TPT状态变量滤波器 |
| [`svf_tpt_shelf.hpp`](include/qwqdsp/filter/svf_tpt_shelf.hpp) | `SvfTPTShelf` TPT SVF 搁架滤波器 |
| [`thiran_filter.hpp`](include/qwqdsp/filter/thiran_filter.hpp) | `ThiranFilter` Thiran 分数延迟滤波器 |
| [`transpose_sallen_key.hpp`](include/qwqdsp/filter/transpose_sallen_key.hpp) | `TransposeSallenKey` 转置 Sallen-Key 滤波器 |
| [`window_fir.hpp`](include/qwqdsp/filter/window_fir.hpp) | `WindowFIR` 窗函数 FIR 设计 |
| [`fixed/acc_traits.hpp`](include/qwqdsp/filter/fixed/acc_traits.hpp) | `AccTypeTrait` 定点累加器类型特性 |
| [`fixed/df1_biquad.hpp`](include/qwqdsp/filter/fixed/df1_biquad.hpp) | `DF1_Biquad` 定点 DF1 双二阶 |
| [`fixed/df1_biquad_q2.hpp`](include/qwqdsp/filter/fixed/df1_biquad_q2.hpp) | `DF1_Biquad2` Q2 格式定点 DF1 双二阶 |
| [`fixed/df1_biquad_split.hpp`](include/qwqdsp/filter/fixed/df1_biquad_split.hpp) | `DF1_BiquadSplit` 拆分式定点 DF1 双二阶 |
| [`fixed/gold_rader.hpp`](include/qwqdsp/filter/fixed/gold_rader.hpp) | `GoldRader` 定点 Gold-Rader 滤波器 |

### 〰️ 振荡器模块 (Oscillator Modules)

| 文件 | 说明 |
|------|------|
| [`blep_coeff.hpp`](include/qwqdsp/oscillator/blep_coeff.hpp) | `Triangle` 三角窗 BLEP 校正系数<br>`Hann` Hann 窗 BLEP 校正系数<br>`BSpline` B 样条 BLEP 校正系数<br>`BlackmanNutall` Blackman-Nuttall BLEP 校正系数<br>`BlackmanNutallApprox` 近似 Blackman-Nuttall BLEP 校正系数 |
| [`blit.hpp`](include/qwqdsp/oscillator/blit.hpp) | `Blit` BLIT 带宽受限脉冲振荡器 |
| [`blit_pwm.hpp`](include/qwqdsp/oscillator/blit_pwm.hpp) | `BlitPWM` BLIT PWM 振荡器 |
| [`coridc_sine_osc.hpp`](include/qwqdsp/oscillator/coridc_sine_osc.hpp) | `CordicSineOsc` CORDIC 正弦振荡器 |
| [`dr_sine_osc.hpp`](include/qwqdsp/oscillator/dr_sine_osc.hpp) | `DROscFull` 数字谐振器正交正弦振荡器 |
| [`dsf.hpp`](include/qwqdsp/oscillator/dsf.hpp) | `DSFClassic` 经典离散求和公式振荡器<br>`DSFComplexFactor` 复数因子版本 |
| [`dsf_correct.hpp`](include/qwqdsp/oscillator/dsf_correct.hpp) | `DSFCorrect` 修正 DSF 振荡器<br>`DSFCorrectComplex` 复数修正版本 |
| [`elliptic_sine_osc.hpp`](include/qwqdsp/oscillator/elliptic_sine_osc.hpp) | `EllipticSineOsc` 椭圆正弦振荡器 |
| [`mcf_sine_osc.hpp`](include/qwqdsp/oscillator/mcf_sine_osc.hpp) | `MCFSineOsc` MCF 正弦振荡器<br>`FullMCFSineOsc` 完整的正交 MCF 振荡器 |
| [`noise.hpp`](include/qwqdsp/oscillator/noise.hpp) | `WhiteNoise` 白噪声<br>`PinkNoise` 粉红噪声<br>`BrownNoise` 布朗噪声<br>`Clicks` 点击噪声 |
| [`polyblep.hpp`](include/qwqdsp/oscillator/polyblep.hpp) | `PolyBlep` 多项式 BLEP 振荡器 |
| [`polyblep_sync.hpp`](include/qwqdsp/oscillator/polyblep_sync.hpp) | `PolyBlepSync` 硬同步 BLEP 振荡器 |
| [`raw_oscillor.hpp`](include/qwqdsp/oscillator/raw_oscillor.hpp) | `RawOscillor` 基础波形生成振荡器 |
| [`smooth_noise.hpp`](include/qwqdsp/oscillator/smooth_noise.hpp) | `SmoothNoise` 插值平滑噪声发生器 |
| [`table_sine_v2.hpp`](include/qwqdsp/oscillator/table_sine_v2.hpp) | `TableSineV2` 线性插值查表正弦振荡器 |
| [`table_sine_v3.hpp`](include/qwqdsp/oscillator/table_sine_v3.hpp) | `TableSineV3` 抛物线插值查表正弦振荡器 |
| [`vic_sine_osc.hpp`](include/qwqdsp/oscillator/vic_sine_osc.hpp) | `VicSineOsc` Vicanek 正交正弦振荡器 |

### 📊 频谱分析模块 (Spectral Modules)

| 文件 | 说明 |
|------|------|
| [`complex_fft.hpp`](include/qwqdsp/spectral/complex_fft.hpp) | `ComplexFFT` 复数 FFT 接口，自动选择 OOURA/IPP 后端 |
| [`ipp_complex_fft.hpp`](include/qwqdsp/spectral/ipp_complex_fft.hpp) | `IppComplexFFT` Intel IPP 复数 FFT 封装 |
| [`ipp_real_fft.hpp`](include/qwqdsp/spectral/ipp_real_fft.hpp) | `IppRealFFT` Intel IPP 实数 FFT 封装 |
| [`oouras_complex_fft.hpp`](include/qwqdsp/spectral/oouras_complex_fft.hpp) | `OourasComplexFFT` OOURA 复数 FFT 裸封装 |
| [`oouras_real_fft.hpp`](include/qwqdsp/spectral/oouras_real_fft.hpp) | `OourasRealFFT` OOURA 实数 FFT 裸封装 |
| [`real_fft.hpp`](include/qwqdsp/spectral/real_fft.hpp) | `RealFFT` 实数 FFT 核心接口，自动选择 OOURA/IPP 后端 |
| [`reassignment.hpp`](include/qwqdsp/spectral/reassignment.hpp) | `Reassignment` 频谱重分配时频细化<br>`ReassignmentCorrect` 导数窗口精确频谱重分配 |

### 🎚️ 音频效果模块 (FX Modules)

| 文件 | 说明 |
|------|------|
| [`delay_line.hpp`](include/qwqdsp/fx/delay_line.hpp) | `DelayLine` 可调延迟线效果 |
| [`elliptic_blep.hpp`](include/qwqdsp/fx/elliptic_blep.hpp) | `EllipticBlep` 椭圆滤波器 BLEP 校正 |
| [`limiter.hpp`](include/qwqdsp/fx/limiter.hpp) | `SimpleLimiter` 峰值保持前视限制器 |
| [`oversample.hpp`](include/qwqdsp/fx/oversample.hpp) | `Oversample` 多级半带过采样/降采样器 |
| [`pitch_shifter.hpp`](include/qwqdsp/fx/pitch_shifter.hpp) | `PitchShifter` 环形缓冲音高移位器 |
| [`pitch_shifter2.hpp`](include/qwqdsp/fx/pitch_shifter2.hpp) | `PhaseVocoder` 相位声码器音高移位 |
| [`plat_reverb.hpp`](include/qwqdsp/fx/plat_reverb.hpp) | `PlateReverb` Dattorro 板式混响 |
| [`polyphase_resample_fir.hpp`](include/qwqdsp/fx/polyphase_resample_fir.hpp) | `PolyphaseDownsamplerFir` 多相下采样 FIR<br>`PolyphaseUpsamplerFir` 多相上采样 FIR |
| [`resample.hpp`](include/qwqdsp/fx/resample.hpp) | `Resample` 重采样器 |
| [`resample_iir.hpp`](include/qwqdsp/fx/resample_iir.hpp) | `ResampleIIR` IIR 重采样器 |
| [`resample_iir_dynamic.hpp`](include/qwqdsp/fx/resample_iir_dynamic.hpp) | `ResampleIIRDynamic` 可变采样率 IIR 重采样 |
| [`uniform_convolution.hpp`](include/qwqdsp/fx/uniform_convolution.hpp) | `UniformConvolution` FFT 均匀分块卷积 |

### ⚡ SIMD 优化模块(已弃用) (SIMD Element Modules)

| 文件 | 说明 |
|------|------|
| [`align_allocator.hpp`](include/qwqdsp/simd_element/align_allocator.hpp) | `AlignedAllocator` 内存对齐分配器 |
| [`algebraic_waveshaper.hpp`](include/qwqdsp/simd_element/algebraic_waveshaper.hpp) | `AlgebraicWaveshaper` SIMD 代数波形塑形器 |
| [`biquads.hpp`](include/qwqdsp/simd_element/biquads.hpp) | `Biquads` SIMD 并行双二阶滤波器组 |
| [`delay_allpass.hpp`](include/qwqdsp/simd_element/delay_allpass.hpp) | `DelayAllpass` SIMD 全通延迟线 |
| [`delay_line_mono.hpp`](include/qwqdsp/simd_element/delay_line_mono.hpp) | `DelayLineMono` SIMD 单声道多延迟线 |
| [`delay_line_multiple.hpp`](include/qwqdsp/simd_element/delay_line_multiple.hpp) | `DelayLineMultiple` SIMD 多通道独立延迟线 |
| [`delay_line_single.hpp`](include/qwqdsp/simd_element/delay_line_single.hpp) | `DelayLineSingle` SIMD 多声道单延迟线 |
| [`delay_line_stereo.hpp`](include/qwqdsp/simd_element/delay_line_stereo.hpp) | `DelayLineStereo` SIMD 立体声延迟线 |
| [`envelope_follower.hpp`](include/qwqdsp/simd_element/envelope_follower.hpp) | `EnevelopeFollower` SIMD 包络跟随器 |
| [`one_pole_tpt.hpp`](include/qwqdsp/simd_element/one_pole_tpt.hpp) | `OnePoleTPT` SIMD 一阶 TPT 滤波器 |
| [`one_pole_tpt_shelf.hpp`](include/qwqdsp/simd_element/one_pole_tpt_shelf.hpp) | `OnepoleTPTShelf` SIMD 一阶 TPT 搁架滤波器 |
| [`plate_reverb.hpp`](include/qwqdsp/simd_element/plate_reverb.hpp) | `PlateReverb` SIMD Dattorro 板式混响 |
| [`simd_pack.hpp`](include/qwqdsp/simd_element/simd_pack.hpp) | `Pack4Bytes` 泛型 SIMD 向量模板<br>`PackOps` SIMD 运算包装 |
| [`simde_wrap4.hpp`](include/qwqdsp/simd_element/simde_wrap4.hpp) | `V4f` SIMDe SSE4.1 128 位浮点向量<br>`V4i` SIMDe SSE4.1 128 位整数向量 |
| [`simde_wrap8.hpp`](include/qwqdsp/simd_element/simde_wrap8.hpp) | `V8f` SIMDe 256 位浮点向量<br>`V8i` SIMDe 256 位整数向量 |
| [`stereo_iir_hilbert_cpx.hpp`](include/qwqdsp/simd_element/stereo_iir_hilbert_cpx.hpp) | `StereoIIRHilbertCpx` 8 级 SIMD 立体声 IIR Hilbert<br>`StereoIIRHilbertDeeperCpx` 16 级 SIMD 立体声 IIR Hilbert |

### 🪟 窗函数模块 (Window Modules)

| 文件 | 说明 |
|------|------|
| [`blackman.hpp`](include/qwqdsp/window/blackman.hpp) | `Blackman` Blackman 窗（旁瓣 -58dB） |
| [`blackman_harris.hpp`](include/qwqdsp/window/blackman_harris.hpp) | `BlackmanHarris` 4 项 Blackman-Harris 窗（旁瓣 -92dB） |
| [`blackman_harris_3term.hpp`](include/qwqdsp/window/blackman_harris_3term.hpp) | `BlackmanHarrisThreeTerm` 3 项 Blackman-Harris 窗（旁瓣 -71.48dB） |
| [`blackman_nuttall.hpp`](include/qwqdsp/window/blackman_nuttall.hpp) | `BlackmanNuttall` Blackman-Nuttall 窗（旁瓣 -98.3dB） |
| [`hamming.hpp`](include/qwqdsp/window/hamming.hpp) | `Hamming` Hamming 窗（旁瓣 -43.8dB） |
| [`hann.hpp`](include/qwqdsp/window/hann.hpp) | `Hann` Hann 窗（旁瓣 -31.6dB） |
| [`kaiser.hpp`](include/qwqdsp/window/kaiser.hpp) | `Kaiser` Kaiser 窗（beta 参数可调） |
| [`lanczos.hpp`](include/qwqdsp/window/lanczos.hpp) | `Lanczos` Lanczos 窗（旁瓣 -26.6dB） |
| [`nuttall.hpp`](include/qwqdsp/window/nuttall.hpp) | `Nuttall` Nuttall 窗（旁瓣 -93.3dB） |
| [`taylor.hpp`](include/qwqdsp/window/taylor.hpp) | `Taylor` Taylor 窗（旁瓣电平 + nbars 参数） |
| [`helper.hpp`](include/qwqdsp/window/helper.hpp) | `Helper` 窗函数工具（归一化、时间加权、零填充） |

### 📐 插值模块 (Interpolation Modules)

| 文件 | 说明 |
|------|------|
| [`catmull_rom_spline.hpp`](include/qwqdsp/interpolation/catmull_rom_spline.hpp) | `CatmullRomSpline` Catmull-Rom 样条插值 |
| [`linear.hpp`](include/qwqdsp/interpolation/linear.hpp) | `Linear` 线性插值 |
| [`makima.hpp`](include/qwqdsp/interpolation/makima.hpp) | `Makima` 修正 Akima 插值 |
| [`sppchip.hpp`](include/qwqdsp/interpolation/sppchip.hpp) | `SPPCHIP` 单调保形三次 Hermite 插值 |

### ✂️ 分段处理模块 (Segmentation Modules)

| 文件 | 说明 |
|------|------|
| [`analyze.hpp`](include/qwqdsp/segement/analyze.hpp) | `Analyze` 分块分析基类 |
| [`analyze_auto.hpp`](include/qwqdsp/segement/analyze_auto.hpp) | `AnalyzeAuto` 自动帧处理循环 |
| [`analyze_synthsis_offline.hpp`](include/qwqdsp/segement/analyze_synthsis_offline.hpp) | `AnalyzeSynthsisOffline` 离线分析综合 |
| [`analyze_synthsis_online.hpp`](include/qwqdsp/segement/analyze_synthsis_online.hpp) | `AnalyzeSynthsisOnline` 在线分析综合 |
| [`mono_reader.hpp`](include/qwqdsp/segement/mono_reader.hpp) | `MonoReader` 单声道文件读取器 |
| [`slice.hpp`](include/qwqdsp/segement/slice.hpp) | `Slice1D` 一维数据分块迭代<br>`Slice2D` 二维数据分块迭代 |

### 🔄 自适应滤波模块 (Adaptive Filter Modules)

| 文件 | 说明 |
|------|------|
| [`burg_lp.hpp`](include/qwqdsp/adaptive/burg_lp.hpp) | `BurgLP` Burg 法线性预测系数估计 |
| [`lag_buffer.hpp`](include/qwqdsp/adaptive/lag_buffer.hpp) | `LagBuffer` 动态延迟缓冲区<br>`LagBufferStatic` 静态延迟缓冲区 |
| [`nlms.hpp`](include/qwqdsp/adaptive/nlms.hpp) | `NLMS` 归一化最小均方自适应滤波器 |
| [`rls_filter.hpp`](include/qwqdsp/adaptive/rls_filter.hpp) | `RLSFIlter` 递归最小二乘自适应滤波器 |

### 🎵 基音检测模块(没一个好用的) (Pitch Detection Modules)

| 文件 | 说明 |
|------|------|
| [`fast_yin.hpp`](include/qwqdsp/pitch/fast_yin.hpp) | `FastYin` 快速 YIN 基音检测 |
| [`helmholtz.hpp`](include/qwqdsp/pitch/helmholtz.hpp) | `Helmholtz` Helmholtz 基音检测 |
| [`mpm.hpp`](include/qwqdsp/pitch/mpm.hpp) | `MPM` McLeod Pitch Method 基音检测 |
| [`pitch.hpp`](include/qwqdsp/pitch/pitch.hpp) | `Pitch` 基音数据结构（pitch_hz / non_period_ratio） |
| [`yin.hpp`](include/qwqdsp/pitch/yin.hpp) | `Yin` YIN 基音检测 |

### 🔢 Cephes 数学函数模块 (Cephes Math Functions)

| 文件 | 说明 |
|------|------|
| [`bessel.hpp`](include/qwqdsp/cephes/bessel.hpp) | `Bessel` 修正 Bessel 函数 I₀ / I₁（双精度） |
| [`besself.hpp`](include/qwqdsp/cephes/besself.hpp) | `Besself` 修正 Bessel 函数 I₀ / I₁（单精度） |
| [`elliptic.hpp`](include/qwqdsp/cephes/elliptic.hpp) | `Elliptic` 完全/不完全椭圆积分（双精度） |
| [`ellipticf.hpp`](include/qwqdsp/cephes/ellipticf.hpp) | `Ellipticf` 完全/不完全椭圆积分（单精度） |

### 🧩 杂项工具模块 (Miscellaneous Modules)

| 文件 | 说明 |
|------|------|
| [`ampd_peak.hpp`](include/qwqdsp/misc/ampd_peak.hpp) | `AMPDPeakFinding` 自动多尺度峰值检测 |
| [`crossover.hpp`](include/qwqdsp/misc/crossover.hpp) | `CrossoverGain` 增益渐变分频<br>`CrossoverPower` 功率渐变分频 |
| [`envelope_follower.hpp`](include/qwqdsp/misc/envelope_follower.hpp) | `EnevelopeFollower` 包络跟随器 |
| [`integrator.hpp`](include/qwqdsp/misc/integrator.hpp) | `IntegratorNaive` 朴素积分器<br>`IntegratorNaiveLeak` 泄漏积分器<br>`IntegratorTrapezoidal` 梯形积分器 |
| [`peakfind.hpp`](include/qwqdsp/misc/peakfind.hpp) | `PeakFinding` 峰值查找 |
| [`smoother.hpp`](include/qwqdsp/misc/smoother.hpp) | `ExpSmoother` 指数平滑器<br>`ConstantTimeSmoother` 恒定时长平滑器 |
