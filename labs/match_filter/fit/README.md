# 任意幅频曲线拟合为数字 biquad 级联

输入覆盖 DC 到 Nyquist 的频率–dB 锚点，PCHIP 插值得到连续目标，也可以直接给出返回有限 dB 的函数。输出是可直接用于数字滤波的 SOS 系数，每节的分子、分母系数自由求解。

## 方法

令 `t = cos(ω)`。实系数二阶节的功率是两个关于 `t` 的二次多项式之比，级联 `K` 节后是 `2K/2K` 型有理函数。去掉分子分母的公用比例后，用 `4K+1` 个节点即可构造齐次线性插值系统。

分子和分母在全带须同号且无零点。通过 `t = (z + 1/z)/2` 将每个根映射成一对互为倒数的根，取单位圆内的根构造稳定、最小相位的 SOS，再校准整体增益。

BRASIL 移动频率插值节点，使各段峰值等高。实验在其后用峰值对节点位置的数值雅可比修正节点。波纹仍不齐时，对 dB 误差上限做二分，每一步在频率网格上解线性不等式 `P_target · 10^(-E/10) · D ≤ N ≤ P_target · 10^(E/10) · D`。取误差过零点重启节点搜索，比较实际 SOS 的峰值误差。低频端用保持有理次数的分式线性坐标变换改善条件数。[BRASIL 原论文](https://www.ricam.oeaw.ac.at/files/reports/20/rep20-37.pdf)

## 实测

48 kHz 下测试四条手工给定的曲线。独立的密集线性与对数频率网格用于复算峰值误差及符号交替，所有结果的极点半径均小于 1。

各格为全带峰值误差（dB）和近峰值的正负交替数。

| 目标曲线 | 1 节 | 2 节 | 3 节 |
|---|---:|---:|---:|
| 宽缓轮廓 | 0.6167 · 6/6 | 0.1487 · 10/10 | 0.0892 · 14/14 |
| 交替峰谷 | 3.7606 · 6/6 | 2.4298 · 10/10 | 0.3732 · 14/14 |
| 低频细节（40–320 Hz） | 3.0045 · 6/6 | 1.5381 · 10/10 | 0.5119 · 14/14 |
| 窄凹口（10.1 kHz，−22 dB） | 4.5646 · 1/6 | 1.9996 · 2/10 | 0.8799 · 14/14 |

对于无次数亏损的 `2K/2K` 逼近，等波纹检查需要 `4K+2` 个等幅、正负交替的极值。以峰值的 99% 为门槛，12 组中有 10 组达到所需个数。低频细节的 1 节解刚好达到门槛，波纹高度差为 0.0300 dB。窄凹口的 1、2 节仍未通过检查。网格上的近似交替是数值证据。实验只拟合有限的 dB 幅度，精确零点需要额外结构约束。

## 工业方法对照

MathWorks 和 NI 都提供任意幅度 IIR 的 least‑p 设计，其中较大的 `p` 逼近峰值误差最小化。两家官方工具箱仅作为方法参考，下面的实测使用独立的 SciPy 对照。[MathWorks Least Pth-Norm Optimal IIR Filter Design](https://www.mathworks.com/help/dsp/ug/least-pth-norm-optimal-iir-filter-design.html) · [NI Digital Filter Design Toolkit 手册](https://download.ni.com/support/manuals/371353a.pdf)

对照在同一 dB 目标上直接优化稳定 SOS 的自由系数。采用多初值最小二乘，再逐步将误差范数从 `p=2` 提高到 `128`，按拟合网格上的峰值选解。三种结果共用独立的密集验证网格。RMS 按线性频率积分，作为峰值之外的第二项指标。

| 目标（3 节） | 最小二乘峰值 | Least‑p 峰值 | BRASIL + 网格峰值 |
|---|---:|---:|---:|
| 宽缓轮廓 | 0.269 dB | 0.120 dB | 0.089 dB |
| 交替峰谷 | 0.877 dB | 0.380 dB | 0.373 dB |
| 低频细节 | 1.560 dB | 0.522 dB | 0.512 dB |
| 窄凹口 | 5.626 dB | 0.972 dB | 0.880 dB |

全 12 组中，BRASIL + 网格有 11 组峰值低于 least‑p。例外是 2 节窄凹口，least‑p 为 1.972 dB，BRASIL + 网格为 2.000 dB。与 least‑p 和 BRASIL + 网格分别比较，最小二乘都有 11 组 RMS 更低。本机单次对照的逐组耗时比中位数为 3.37，BRASIL + 网格耗时更长。所有结果的极点均在单位圆内。

```powershell
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
python qwqdsp/labs/match_filter/fit/compare_lpnorm.py
```

结果写入 `output/industry_comparison.csv`。使用 `--target` 或 `--sos` 时，文件名会包含筛选条件。

有 MATLAB DSP System Toolbox 时，可用下列命令测宽缓轮廓的 3 节官方实现。传入线性幅度，并用目标幅度的倒数给相对误差加权。验证仍使用 dB 峰值。[iirlpnorm 接口](https://www.mathworks.com/help/dsp/ref/iirlpnorm.html)

```matlab
f = unique([linspace(0,1,769), logspace(-6,0,257)]);
hz = [0 400 1500 3400 6800 12300 19400 24000];
db = [2 1 3 7 3 -2 -4 -3];
a = 10.^(pchip(hz,db,24000*f)/20);
[b,d] = iirlpnorm(6,6,f,[0 1],a,1./a,[2 128]);
v = unique([linspace(0,1,65538), logspace(-7,0,8193)]);
h = freqz(b,d,pi*v(:));
reference = pchip(hz,db,24000*v);
e = 20*log10(abs(h(:))) - reference(:);
peak_db = max(abs(e))
```

## 分式 Remez

分式 Remez 直接在参考极值上交替施加 `+E` 和 `−E` dB。每轮解出波纹高度 `E` 与对应的功率有理函数，谱分解后根据 SOS 的实际误差交换极值。初值取 BRASIL + 网格结果，交替下界仍有明显间隙时也尝试 least‑p 初值。

四条目标、1–3 节共 12 组均得到稳定的近等波纹解。10.1 kHz 窄凹口的 1、2、3 节峰值误差分别为 4.4372、1.8982、0.8783 dB，其中 2 节需要 least‑p 初值。

```powershell
python qwqdsp/labs/match_filter/fit/remez_experiment.py
```

结果写入 `output/remez.csv`。使用 `--target` 时，文件名会包含目标编号。

用 `from remez_fit import FitRemez` 调用求解器。接口为 `FitRemez(target_db, sections)`，返回 SOS 数组。目标函数接收 DC 至 Nyquist 归一化频率，并返回 dB。

## 运行

依赖 NumPy、SciPy，画图另需 Matplotlib。在工作区根目录执行

```powershell
python qwqdsp/labs/match_filter/fit/experiment.py --plot
```

图保存为 `output/curves.png`。

在实验目录中可用以下方式拟合自己的锚点。SOS 每行为 `b0, b1, b2, 1, a1, a2`。

```python
from fit_curve import FitCurve, SampledCurve

target_db = SampledCurve([0, 400, 1500, 3400, 6800, 12300, 19400, 24000],
                          [2, 1, 3, 7, 3, -2, -4, -3], 48000)
result = FitCurve(target_db, sections=2)
if result is not None:
    sos = result.sos
```
