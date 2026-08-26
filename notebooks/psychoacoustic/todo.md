这份心理声学全景学习路线图涵盖了从**解剖生理学、心理物理学、数学建模到工业级 DSP 算法实现**的完整知识体系。

路线图采用 Markdown To-Do List 形式，并针对每个模块标注了**核心理论**与对应的**实际应用程序 / 工程落地场景**。

---

## 阶段一：听觉生理与心理物理学基础（Perceptual Fundamentals）

* [ ] ### 1.1 耳蜗力学与听觉感知基石
* [ ] **物理与生理结构**：外耳/中耳传输特性、耳蜗基底膜（Basilar Membrane）行波理论、毛细胞换能机制。
* [ ] **频域-位置映射**：听觉滤波器组（Auditory Filterbank）、 Place Theory（位置理论）与 Phase-locking（相位锁定）。
* [ ] **应用程序/工程例子**：
* **人工耳蜗（Cochlear Implants）信号处理器**：基于基底膜频域映射的带通滤波与电刺激编码。
* **WebAudio API / DAWs**：基于听觉特性的对数频域分析仪（Logarithmic Spectrum Analyzer）。




* [ ] ### 1.2 临界频带与频域标尺 (Critical Bands & Scales)
* [ ] **频带模型**：Bark Scale、ERB Scale (Equivalent Rectangular Bandwidth)、Mel Scale 的导出与数学公式。
* [ ] **临界带宽测定**：通过带通噪声掩蔽纯音实验理解临界带宽（Critical Bandwidth）。
* [ ] **应用程序/工程例子**：
* **MFCC（梅尔频率倒谱系数）提取器**：语音识别（ASR）与音乐特征检索（MIR）的前端特征工程。
* **感知 EQ（Perceptual Equalizer）**：按 Bark Scale 划频段的动态均衡器设计。




* [ ] ### 1.3 强度感知与等响度 (Loudness Perception)
* [ ] **感知曲线**：Fletcher-Munson & ISO 226 等响度曲线（Equal-Loudness Contours）。
* [ ] **响度单位与量化**：dB SPL、Phon（响度级）与 Sone（主观响度）的推导与转换。
* [ ] **应用程序/工程例子**：
* **夜间模式 / 动态等响度补偿（Dynamic Loudness Compensation）**：如 Apple HomePod、Yamaha YPAO 在低音量播放时动态提升低频与高频。
* **广播级响度限制器（Loudness Limiter）**：电视与流媒体广播的峰值与响度控制。





---

## 阶段二：感知掩蔽机制与感知编码（Masking & Audio Coding）

* [ ] ### 2.1 频域与时域掩蔽 (Frequency & Temporal Masking)
* [ ] **频域掩蔽**：单频掩蔽（Tone-masking-Tone）、噪声掩蔽（Tone-masking-Noise）、扩散函数（Spreading Function）。
* [ ] **时域掩蔽**：前掩蔽（Pre-masking, 5~20ms）与后掩蔽（Post-masking, 50~200ms）机制。
* [ ] **绝对听阈（ATH）**：Quiet Threshold 曲线在算法中的表达。
* [ ] **应用程序/工程例子**：
* **Psychoacoustic Dither（心理声学抖动）**：如 iZotope Ozone 导出 16-bit CD 音频时，将量化噪声推入高频掩蔽区/高听阈区。
* **降噪与语音增强（Spectral Subtraction with Masking Threshold）**：只滤除高于掩蔽阈值的噪声，消除“音乐噪声（Musical Noise）”。




* [ ] ### 2.2 心理声学模型与音频压缩 (Psychoacoustic Models in Codecs)
* [ ] **MPEG Psychoacoustic Model 1 & 2**：音调/非音调分量提取、掩蔽曲线叠加、SMR（信掩比）计算。
* [ ] **感知加权与比特分配**：利用 SMR 和 NMR（噪声掩蔽比）指导霍夫曼编码与量化阶梯。
* [ ] **预冲激控制**：长短窗动态切换（Window Switching）解决前掩蔽时域失真。
* [ ] **应用程序/工程例子**：
* **经典与现代感知编码器**：MP3 (LAME)、AAC (Fraunhofer)、Opus (Xiph.Org)、Dolby Digital (AC-3)。
* **蓝牙音频传输协议**：LDAC、aptX Adaptive 等低延迟高压缩比编解码器。





---

## 阶段三：空间听觉与 3D 音频（Spatial Hearing & Binaural Audio）

* [ ] ### 3.1 双耳定位线索 (Binaural Localization Cues)
* [ ] **双重理论（Duplex Theory）**：低频依赖 ITD（双耳时间差，<1.5kHz），高频依赖 ILD（双耳声级差，>1.5kHz）。
* [ ] **头影效应与衍射**：球形头模型与实际人体解剖结构的声学散射。
* [ ] **应用程序/工程例子**：
* **立体声拓宽插件（Stereo Widener / Panner）**：利用微秒级延迟（Haas 效应）与高频 EQ 差值模拟左右空间感。




* [ ] ### 3.2 HRTF 与三维空间渲染 (HRTF & Binaural Rendering)
* [ ] **HRTF/HRIR**：头部相关传输函数，耳廓（Pinna）对高频的谱谷（Spectral Cues）与垂直面（Elevation）定位。
* [ ] **SOFA 文件格式**：Spatially Oriented Format for Acoustics 标准与其解析。
* [ ] **动态追踪**：头部追踪（Head Tracking）下的双耳实时卷积算法。
* [ ] **应用程序/工程例子**：
* **空间音频渲染引擎**：Apple Spatial Audio、Sony 360 Reality Audio、Dolby Atmos Binaural Renderer。
* **游戏空间音频中间件**：Wwise (Spatial Audio)、FMOD、Steam Audio、Meta Audio SDK。




* [ ] ### 3.3 房间声学与多声道重放 (Room Perception & Multi-channel)
* [ ]  precedence effect（哈斯效应/优先效应）：直达声与早期反射声的融合时间窗（1~30ms）。
* [ ] **去相关性与 ASW**：IACC（双耳互相关系数）与空间宽广感（Apparent Source Width）。
* [ ] **应用程序/工程例子**：
* **虚拟房间建模插件**：Waves Abbey Road Studio、Plugin Alliance DearVR PRO。
* **汽车音响声场重构**：Bose / Harman Kardon 的车载虚拟环绕声算法。





---

## 阶段四：高级心理声学品质与主客观评价（Sound Quality & Measurement）

* [ ] ### 4.1 高级心理声学描述符 (Psychoacoustic Descriptors)
* [ ] **粗糙度（Roughness, Asper）**：20~300 Hz 振幅调制产生的刺耳感。
* [ ] **波动度（Fluctuation Strength, Vacil）**：<20 Hz 慢速调制的起伏感。
* [ ] **锐度（Sharpness, Acum）**：高频能量比例与频谱质心分析。
* [ ] **音调感（Tonalness / Tonality）**：声音中正弦分量与随机噪声的比例。
* [ ] **应用程序/工程例子**：
* **NVH 汽车声学评估软件**：Head Acoustics ArtemiS SUITE、Siemens Simcenter Testlab。
* **消费电子声学调校**：智能音箱、耳机生产线的自动声学质检与 Sound Profile 调校。




* [ ] ### 4.2 客观主观化测量标准 (Standardized Perceptual Metrics)
* [ ] **ITU-R BS.1770-4**：K-Weighting 滤波器与 LUFS (Loudness Units Full Scale) 计算。
* [ ] **PEAQ (ITU-R BS.1387)**：音频质量主观感知测量标准（ODG 得分预测）。
* [ ] **POLQA / PESQ**：语音质量感知评估标准（Telecommunication Speech Quality）。
* [ ] **应用程序/工程例子**：
* **流媒体响度标准化系统**：Spotify, Apple Music, YouTube 的自动响度 Normalize 算法。
* **音频编解码算法自动测试 pipeline**：自动化评估压缩算法后的 ODG/PESQ 分数。





---

## 阶段五：心理声学信号处理应用（Psychoacoustic DSP Processing）

* [ ] ### 5.1 虚拟听感与硬件极限突破 (Virtual Psychoacoustic DSP)
* [ ] **Missing Fundamental（缺失基频效应）**：利用高阶非线性谐波生成（Harmonic Generator）塑造“虚拟低音（Virtual Bass）”。
* [ ] **高频重建（SBR, Spectral Band Replication）**：利用低频信息与控制参数伪造高频包络。
* [ ] **应用程序/工程例子**：
* **微型扬声器低频增强算法**：Waves MaxxBass、Dirac Panorama Sound。
* **低码率编码优化**：HE-AAC 中的 SBR 模块、Bluetooth LC3 编码增强。




* [ ] ### 5.2 听觉场景分析与智能分离 (Auditory Scene Analysis - ASA)
* [ ] **鸡尾酒会效应与听觉流（Auditory Stream）**：大脑按照谐波关系、共同起伏（Common Fate）、频域连续性进行声源分离机制。
* [ ] **应用程序/工程例子**：
* **智能降噪与 AI 提音**：Krisp.ai、iZotope RX Spectral Repair（结合传统心理声学约束与神经网络）。
* **助听器 DSP 引擎**：Oticon、Phonak 助听器中的动态声场感知与人声增强算法。





---

## 🛠 推荐工具与实战演练环境

| 类别 | 工具/库 | 用途说明 |
| --- | --- | --- |
| **Python 分析库** | `librosa`, `scipy.signal`, `pyroomacoustics` | 提取 MFCC、计算 Spectrogram、模拟房间脉冲响应 |
| **心理声学量化库** | `MOSQITO` (Python), `SQiLab` | 开源心理声学指标计算库（Loudness, Sharpness, Roughness） |
| **空间音频开发** | `SOFA API`, `VST3 / JUCE SDK` | 解析 SOFA 格式，在 C++ JUCE 中编写 HRTF 双耳卷积插件 |
| **标准参考源码** | `LAME MP3 Source Code`, `Opus Codec` | 深入研读工业级 C 语言实现的心理声学模型与比特分配逻辑 |
