# vicanek — `match_biquad` 的移植与对照

`match_biquad.hpp`（Vicanek：脉冲响应不变极点 + 分子拟合）的 Python 移植与
和 Ivantsov 双点解的对照。映射定义与 Ivantsov 那边一致：给定设计，找模拟频率
`y` 使 `|Ha(y)| = |Hd(x)|`；参考实现见 `../ivantsovy/README.md`。

与 Ivantsov 的根本区别：**它不是频率变换**。极点用脉冲响应不变
（`a1 = -2 exp(-ζwc) cos(√(1-ζ²)wc)`、`a2 = exp(-2ζwc)`，角度就是 `wc`），零点
再由 `|H|²` 在 2–3 个频点上的线性方程解出（`|H|²` 是
`φ = (1-s, s, 4s(1-s))`、`s = sin²(w/2)` 的二次型，所以三个频点就是一个 3×3
线性方程组）。于是：

- **没有全局映射**，每种类型各有一条曲线；
- 低通/高通只匹配 DC 与 `wc`（Nyquist 不匹配），映射右端停在 0.84–0.85·Nyquist；
- shelf 在三点解可行时匹配 DC / Nyquist / `wc`，右端才回到 1；极点频率
  `wc·√G` 越过 Nyquist 时（如 +12 dB、`fc > 0.354·fs`）退化为重分布匹配点，
  放弃 Nyquist 精确；
- 映射几乎与 `fc` 无关（极点角度就是 `wc`，拟合是局部的），所以高频档不会像
  Ivantsov 那样塌陷。

## 映射曲线

| fc / (fs/2) | Vicanek 低通 | Vicanek 高通 | Vicanek shelf +12 | Ivantsov（任意类型） |
|---|---|---|---|---|
| 0.2 | 0.8380 | 0.8526 | 1.0000 | 1.0000 |
| 0.9 | 0.8404 | 0.8550 | 1.0000 | 1.0000 |

（表内是 x → Nyquist 处映射到的模拟频率 / (fs/2)。）

图：`output/vicanek_mapping.png`（fc = 0.2 / 0.9·fs/2 两张映射曲线；右图可见
Vicanek 的三条曲线随 fc 几乎不动，而 Ivantsov 那条塌了）。

## 精度对照

全带 `max|ΔdB|`（fs = 48 kHz，Q = 1/√2；低通/高通各用自身原型，shelf 用同一
模拟原型比较）：

| fc | V-LP | I-LP | V-HP | I-HP | V-shelf | I-shelf |
|---|---|---|---|---|---|---|
| 20 Hz | 1.05 | 0.000 | 1.93 | 0.000 | 0.000 | 0.000 |
| 117.6 Hz | 3.07 | 0.000 | 1.93 | 0.000 | 0.000 | 0.000 |
| 691.4 Hz | 3.07 | 0.000 | 1.93 | 0.005 | 0.000 | 0.003 |
| 4065 Hz | 3.05 | 0.003 | 1.88 | 0.164 | 0.010 | 0.087 |
| 23.9 kHz | 0.81 | 2.50 | 2.13 | 37.99 | 1.29 | 4.52 |

读法：Vicanek 的误差与 `fc` 基本无关（低通稳定在 ~3 dB、高通 ~1.9 dB），所以
23.9 kHz 档它明显更好；4 kHz 以下 Ivantsov 几乎精确，而 Vicanek 的低通仍差 3 dB。
shelf 两边都是三点匹配，Vicanek 在中高频更省（0.010 对 0.087、1.29 对 4.52），
低频持平。

## 文件

- `vicanek.py`：`match_biquad.hpp` 的幅度部分移植（`impulse_invariant` /
  `get_phi` / `fit_numerator` / `lowpass` / `highpass` / `highshelf` / `lowshelf` /
  `shelf_roots`）。不变量已核对：`|H(0)|`、`|H(wc)| = Q`、shelf 的 DC 与 `wc`
  精确，Nyquist 仅在三点解可行时精确。
- `exp_vicanek.py`：对照图与上表。

## 运行

```bash
python qwqdsp/labs/match_filter/vicanek/exp_vicanek.py   # 对照图 + 两张表，≈3 s
```
