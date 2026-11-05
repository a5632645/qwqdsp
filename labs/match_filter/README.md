# match_filter — 滤波器匹配方法

各方法一个子目录，代码、结论与运行方式见各自的 README。

| 目录 | 内容 |
|---|---|
| `ivantsovy/` | Ivantsov《On the Ideal Bilinear and Biquadratic Digital Filter》§4 的数值验证，以及由它导出的**双点匹配**解（截止/中点 + Nyquist 同时精确）、该变换的**精确频率映射**。15 张图。 |
| `vicanek/` | Vicanek `match_biquad.hpp` 的移植与对照：脉冲响应不变极点 + 分子拟合，**没有全局映射**、每种类型一条曲线。 |

```bash
python qwqdsp/labs/match_filter/ivantsovy/run_all.py     # Ivantsov 侧全部自检 + 15 张图
python qwqdsp/labs/match_filter/vicanek/exp_vicanek.py   # Vicanek 对照图 + 表
```

两侧都以 `fs = 48 kHz`、截止频率 20 / 117.6 / 691.4 / 4065 / 23900 Hz（对数 5 点）
为基准；图的 `output/` 各自独立、已 gitignore。
