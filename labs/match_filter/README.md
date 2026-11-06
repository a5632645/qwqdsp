# match_filter — 滤波器匹配方法

各方法一个子目录，代码、结论与运行方式见各自的 README。

| 目录 | 内容 |
|---|---|
| `best_fit/` | Vicanek 功率函数匹配后，以差分进化优化二阶低通/高通；二阶高通保留 DC 二阶零点与 fc，放弃 Nyquist 匹配，并比较冲击不变基线。 |
| `best_pole_fit_polynominal/` | 将一阶低通逐点 minimax 分母工程化为归一化截止频率上的分段 Chebyshev 多项式，并生成 C++ 头文件。 |
| `ivantsovy/` | Ivantsov《On the Ideal Bilinear and Biquadratic Digital Filter》§4 的数值验证，以及由它导出的**双点匹配**解（截止/中点 + Nyquist 同时精确）、该变换的**精确频率映射**。15 张图。 |
| `vicanek/` | Vicanek `match_biquad.hpp` 的移植与对照：脉冲响应不变极点 + 分子拟合，**没有全局映射**、每种类型一条曲线。 |
| `best_pole_fit/` | Vicanek 功率函数三点匹配（0 / fc / Nyquist）+ 等波纹极点；一阶低通全带峰值 0.094–0.120 dB，一阶高通 1.5e−7–0.213 dB，高 shelf（±6/±12 dB）2.3e−7–0.148 dB，二阶低通（Q=0.5–2）全带 0.07–0.70 dB（≥−60 dB 频段 0.004–0.67），二阶高通（Vicanek 两点，DC+fc）≤0.477 dB（23.9k），同结构 II 1.24–3.27；fc≤4k ≤0.065，二阶带通（归一化峰值 1）全带 0.05–0.97（23.9k 0.05–0.23），II 0.90–5.98；二阶 notch（DC+fc 两点）≤0.103（fc≤4k），23.9k 1.64–1.77，II 34.8–41.3；二阶 peaking（+12/+6/−6/−12 dB × Q=0.5–2，三点）fc≤4k ≤0.30，23.9k 0.04–0.14，同结构 II（DC+fc 两点）0.41–3.57、cut 档最高 12.2；二阶高 shelf（+12/+6/−6/−12 dB × Q=0.5–2，0/fc+斜率）fc≤4k ≤0.107，23.9k 0.07–2.99（+12/Q=2 2.99），II 0.30–13.7（cut 档 11.8–13.7）。 |
| `constraint_fit/` | 一/二阶低通、高通、带通、peaking、high shelf：结构约束 + 自由节点插值。一阶低通 **0.026–0.080 dB**、高通 **1.1e-7–4.8e-2**；二阶低通 **0.38–0.45**、高通 **0–0.063**、带通 **0.003–0.48**、peaking **0–0.18**、high shelf **0–0.43**。硬约束只保留用途必需项，其余交给自由节点。 |

```bash
python qwqdsp/labs/match_filter/ivantsovy/run_all.py           # Ivantsov 侧全部自检 + 15 张图
python qwqdsp/labs/match_filter/vicanek/exp_vicanek.py         # Vicanek 对照图 + 表
python qwqdsp/labs/match_filter/best_fit/lowpass_1st.py        # 三点匹配 + 差分进化（5 个 fc）
python qwqdsp/labs/match_filter/best_fit/highpass_1st.py       # 一阶高通：三点匹配 + 差分进化（5 个 fc）
python qwqdsp/labs/match_filter/best_fit/highpass_2nd.py       # 二阶高通：DC/fc 匹配 + 极点差分进化（4 个 Q × 5 个 fc）
python qwqdsp/labs/match_filter/best_pole_fit_polynominal/lowpass_1st.py # 一阶低通：minimax 极点多项式拟合 + C++ 头文件
python qwqdsp/labs/match_filter/best_pole_fit/lowpass_1st.py   # 一阶低通：三点匹配 + 最佳极点
python qwqdsp/labs/match_filter/best_pole_fit/highpass_1st.py  # 一阶高通：三点匹配 + 最佳极点
python qwqdsp/labs/match_filter/best_pole_fit/highshelf_1st.py # 高 shelf：4 个 boost × 5 个 fc
python qwqdsp/labs/match_filter/best_pole_fit/lowpass_2nd.py   # 二阶低通：4 个 Q × 5 个 fc
python qwqdsp/labs/match_filter/best_pole_fit/highpass_2nd.py  # 二阶高通：Vicanek 两点 + 极点优化（vs II）
python qwqdsp/labs/match_filter/best_pole_fit/bandpass_2nd.py  # 二阶带通：三点匹配 + 极点优化（vs II）
python qwqdsp/labs/match_filter/best_pole_fit/notch_2nd.py     # 二阶 notch：两点（DC+缺口）+ 极点优化（vs II）
python qwqdsp/labs/match_filter/best_pole_fit/peaking_2nd.py   # 二阶 peaking：boost×Q 三点 + 极点优化（vs II）
python qwqdsp/labs/match_filter/best_pole_fit/highshelf_2nd.py # 二阶高 shelf：boost×Q 三点 + 极点优化（vs II）
python qwqdsp/labs/match_filter/constraint_fit/lowpass_1st_interp.py # 一阶低通：5 点有理插值 + BRASIL
python qwqdsp/labs/match_filter/constraint_fit/highpass_1st_interp.py # 一阶高通：4 点有理插值 + BRASIL
python qwqdsp/labs/match_filter/constraint_fit/lowpass_2nd_interp.py # 二阶低通：Q=0.5–2，4 点插值
python qwqdsp/labs/match_filter/constraint_fit/highpass_2nd_interp.py # 二阶高通：Q=0.5–2，3 点插值
python qwqdsp/labs/match_filter/constraint_fit/bandpass_2nd_interp.py # 二阶带通：Q=0.5–2，4 点插值
python qwqdsp/labs/match_filter/constraint_fit/peaking_2nd_interp.py  # 二阶 peaking：boost×Q
python qwqdsp/labs/match_filter/constraint_fit/highshelf_2nd_interp.py # 二阶 high shelf：boost×Q
python qwqdsp/labs/match_filter/constraint_fit/highshelf_1st_interp.py # 一阶 high shelf：boost
python qwqdsp/labs/match_filter/constraint_fit/lowpass_1st.py        # 一阶低通：硬约束枚举（负结果）
```

各实验都以 `fs = 48 kHz`、截止频率 20 / 117.6 / 691.4 / 4065 / 23900 Hz（对数 5 点）
为基准；图的 `output/` 各自独立、已 gitignore。
