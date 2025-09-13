
# qwqdsp

存储我曾经使用过的东西（重复造轮子！），为了防止每次都要去仓库找代码复制，将它们全部打包到一起

## 贡献

很乐意有人提交pull request（如果有的话就好）

## 物品

- 滤波器
  - 定点滤波器
    - 直接1型二次量化补偿biquad
    - 直接1型一次量化补偿biquad
    - 直接1型ba分离biquad
    - gold-rader结构
  - 全通滤波器
    - Schroeder全通滤波器
    - 一阶（多相）全通滤波器
    - 二阶（多相）全通滤波器
    - Thiran全通滤波器
  - FIR
    - 直接FIR
    - 转置FIR
    - 加窗FIR系数设计
  - 共振峰参数
  - 转置II型biquad
  - gold-rader结构
  - iir希尔伯特滤波器
  - iir复数希尔伯特滤波器
  - iir滤波器设计（无奇数）
    - 巴特沃斯
    - 切比雪夫I型
      - 偶数阶段修正
    - 切比雪夫II型
      - 偶数阶段修正
    - 椭圆滤波器（通带是等波纹但是振幅在某些情况下有点偏差）
  - iir滤波器设计额外
    - 指定截止衰减巴特沃斯
    - 指定截止切比雪夫
  - 格型（多相）滤波器
  - linkwitz-riley分频滤波器
  - 中位数滤波器
  - RBJ-biquad系数计算
  - SVF滤波器
- fx
  - 延迟线
  - 重采样器
    - 凯泽尔窗重采样
    - elliptic blep IIR重采样（非线性相位，速度很快，归一化音量，约-100dB抑制（见elliptic-blep仓库滤波器设计））
    - elliptic blep IIR动态重采样（可能是为采样器实现的）
  - uniform分块卷积
- 插值
  - catmull-rom-spline
  - linear
  - makima
  - sppchip
- 其他（misc）
  - 峰值查找
    - 自适应峰值查找
    - 另一个不知道哪里来的峰值查找（matlab原版）
  - 积分器
    - 普通积分器
    - 梯形积分器
    - 泄露积分器
  - 交叉淡化
    - 恒定增益
    - 恒定功率
  - 参数平滑
    - 恒定时间
    - 恒定变化
    - 一阶指数滤波
- 振荡器
  - blit（宽带限制脉冲串）
    - blit
    - 带有pwm的blit
    - blit经典波形（噪声大）
  - 正弦发生器
    - cordic-正交
    - 数字谐振器（digital resonator）
    - 椭圆-伪正交
    - 修正耦合结构（magic circle）-正交
    - 查表-正交
    - Vicanek/Levine-正交
  - dsf（离散求和公式）
    - 复数dsf
    - 实数dsf
      - 普通dsf
      - 梳状滤波dsf
  - blep
    - 2/4阶blep
    - 经典波形blep
    - hardsync（无三角形）
  - 噪声
    - 白色
    - 粉色
    - 棕色
- 音高跟踪
  - yin
  - fft加速yin
- 切片处理
  - 单声道多声道切片
  - 分析切片
  - 实时分析合成切片
  - 离线分析合成切片
  - 手动分析切片
- 光谱处理
  - 实数fft
  - 复数fft
  - 时频重分频
  - nc（比矩形窗更窄，超过-50dB的旁瓣抑制）<=论文说的
- 窗口函数（普通和导数以及加时间）
  - blackman
  - hamming
  - hann
  - kaiser
  - lanczos
  - taylor
  - 窗函数的相关信息（主瓣宽度，旁瓣抑制等）
- 各种转换函数
- 适用于等距的插值器
- 适用于普通的插值器
- 多项式近似（基于remez）

## TODO

- [ ] DSF
  - [ ] 限制DSF末端的增益
  - [ ] 添加DSF增益的斜率

- [ ] 滤波器
  - [ ] 留数法部分分式分解-级联转并联
  - [ ] SIMD IIR
  - [ ] 椭圆滤波器改进
  - [ ] 奇数阶IIR设计
  - [ ] FIR设计
  - [ ] 等波纹FIR设计
  - [ ] 最小二乘FIR设计
  - [ ] FARROW插值滤波器
  - [ ] FIR希尔伯特变换
  - [ ] comformal变换

- [ ] 振荡器抗混叠
  - [ ] DPW
  - [ ] 优化BLIT的噪声
  - [ ] BLEP三角波硬同步
  - [ ] 优化BLEP代码的分支

- [ ] B-spline
