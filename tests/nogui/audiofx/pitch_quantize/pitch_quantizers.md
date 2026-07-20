# 音高量化器实验版本

本目录中的六个离线测试程序用于比较不同的相位声码器音高量化策略。所有版本都使用 Hann 分析/综合窗、FFT、重叠相加和同一套调性/MIDI 查找逻辑，并输出单声道 WAV。版本 1、2 使用 `SweepWav()`；版本 3 至 6 使用 `drumloop.wav`，跨组比较前应先统一输入素材。

`pitch_quatizer.cpp` 的文件名沿用历史拼写，其中缺少 `n`；它是第一个版本。

| 版本 | 源文件 | 频率/谱形映射 | 相位策略 | 输出文件 | 听感（鼓循环） | 听感（人声） |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | `pitch_quatizer.cpp` | 每个原 bin 留在原位置，仅改变该 bin 的瞬时频率 | 朴素 phase vocoder 递推 | `pitch_quatizer.wav` | 瞬态丢失 | 中频明显量化 |
| 2 | `pitch_quantizer2.cpp` | 每个源 bin 独立映射到最近 MIDI bin，碰撞保留最大幅度 | 目标谱峰递推 + Puckette identity phase lock | `pitch_quantizer2.wav` | 瞬态丢失 | 听起来很流畅 |
| 3 | `pitch_quantizer3.cpp` | 局部谱峰量化，峰域随峰整体平移 | PGHI 源相位差的 phase lock | `pitch_quantizer3.wav` | 瞬态还行 | 低频有 auto tune 的感觉 |
| 4 | `pitch_quantizer4.cpp` | 与版本 3 相同 | 原始分析相位差的 phase lock | `pitch_quantizer4_phase_lock.wav` | 瞬态还行 | 低频有 auto tune 的感觉 |
| 6 | `pitch_quantizer6.cpp` | 与版本 3 相同 | 纯 PGHI 源绝对相位投影，无 phase lock | `pitch_quantizer6_pure_pghi.wav` | 瞬态丢失 + 大量 music noise | 低中频调制感 |

> 版本 5 与版本 3 代码完全相同，已删除。

## 共同处理流程

1. 估计瞬时频率：版本 1 使用相邻 STFT 帧相位差，版本 2 至 6 使用当前帧和一采样延迟帧的 FFT 相位差。
2. 将峰或 bin 的频率查找至最近的、属于当前调性的 MIDI 音符。
3. 重建半谱并 IFFT。
4. 加 Hann 综合窗并 OLA，最后按数值测得的 OLA 峰值归一化。

各测试默认使用 C Major，`SetKey()` 可选择大调、自然小调、和声小调、旋律小调或半音阶。

## 版本 1：朴素 Phase Vocoder

`pitch_quatizer.cpp` 在原始频率 bin 中保留幅度，只调整每个 bin 的综合相位增量，使其向最近 MIDI 频率移动。

- 优点：结构最简单，适合作为基准。
- 限制：修正量受 `+-fs / (2 * hop)` 限制；每个 bin 独立传播相位，容易出现相位涂抹和金属感。

## 版本 2：逐 bin 映射

`pitch_quantizer2.cpp` 把每个源 bin 直接搬到最近 MIDI 频率对应的目标 bin。多个源 bin 落入同一目标时，仅保留幅度最大的源 bin。

- 目标谱峰按目标频率递推；其他目标 bin 锁定到其最近目标峰的线性 FFT 相位偏移。
- 优点：可任意移动目标 bin，不受朴素 PV 修正范围限制。
- 限制：宽带峰域和瞬态会在碰撞中收缩为少数孤立 bin，难以保留局部谱形与相对相位。

## 峰域整体映射

版本 3 至 6 共享此映射。先检测局部谱峰 `p`，将该峰量化到目标 bin `q`；属于该峰域的源 bin `k` 映射为：

```text
target_bin = q + (k - p)
```

峰中心完成音高量化，峰域内部保留频谱宽度与邻接关系。若两个峰域在目标频谱重叠，仍保留幅度较大的源 bin；落出 DC/Nyquist 边界的 bin 会被丢弃。

这一步改善了版本 1 因孤立 bin 导致的频谱破碎问题，但**仅靠峰域映射不足以保留瞬态**——版本 6 虽有相同的映射，仍丢失瞬态。保留瞬态的关键在于 version 3–5 使用的**目标峰相位递推 + phase lock** 的组合。

## 版本 3：峰域映射 + PGHI Phase Lock 原型

`pitch_quantizer3.cpp` 是峰域映射、PGHI 与 phase lock 的工作原型。PGHI 在源频谱域产生相位场，再由目标峰锚定峰域内的相对相位。

目标峰仍按量化后的目标频率递推，因此 PGHI 不直接决定输出的绝对峰相位。

**鼓循环瞬态保留较好；人声听感与版本 4、5 相同。**

## 版本 4：峰域映射 + 原始分析相位 Phase Lock

`pitch_quantizer4.cpp` 是不含 PGHI 的基线。目标峰按 MIDI 目标频率递推；峰域 bin 使用同一分析帧中的相对相位：

```text
phase(target_bin) = phase(target_peak)
                  + wrap(analysis_phase(source_bin)
                       - analysis_phase(source_peak))
```

此版本用于判断峰域映射和 phase lock 本身的效果。

**鼓循环瞬态保留较好；人声听感与版本 3 相同。**

> 原版本 5 与版本 3 代码完全相同，已删除。版本 6 序号不变。

## 版本 6：峰域映射 + 纯 PGHI

`pitch_quantizer6.cpp` 保留峰域映射和源域 PGHI，但不使用目标峰相位递推或 phase lock。每个保留的目标 bin 直接投影源 PGHI 绝对相位：

```text
phase(target_bin) = pghi_phase(source_bin) + frequency_shift(target_bin)

frequency_shift[n] = frequency_shift[n - 1]
                   + (omega_target - omega_source) * hop
```

`frequency_shift` 在目标 bin 非活动后清零。它只把源相位场逐帧平移到量化后的目标频率，不使用目标谱峰锚定或 phase lock。该版本用于隔离“PGHI 绝对相位投影”与“目标峰递推 + phase lock”的听感差异。
**实测瞬态丢失且伴随大量 music noise，与版本 1、2 的瞬态丢失表现不同。**
## 推荐对比顺序

对鼓循环或带瞬态的素材，建议先比较：

1. **版本 2 与版本 4**，检查峰域整体映射带来的改善（版本 2 瞬态丢失 vs 版本 4/3/5 保留瞬态）。
2. **版本 4 与版本 6**，比较 phase lock + 目标峰递推与纯 PGHI 绝对相位投影（两者瞬态表现迥异）。
3. **版本 3/4/5** 三者鼓循环听感一致（瞬态还行），人声听感也一致（低频 auto tune 感），无需在同一组内反复对比。

比较时请保持输入、调性、FFT 大小、hop 和归一化方式一致，并以相同响度试听。
