# 数字 Thiran 低通：指定截止频率的系数公式

问题：固定阶数 $n$，要求 $-3\mathrm{dB}$ 截止频率为 $\omega_c$，不管群延迟 —— 求系数。

## 文献结论：没有现成公式

| 出处 | 给的是什么 |
|---|---|
| Thiran 1971（原文在 IEEE 付费墙后，未能直接核读） | 按**延迟**参数给出全极点/全通的闭式系数 |
| Laakso 等的综述（Aalto 章节，转述 Thiran 1971） | 明确写：*"A drawback of Thiran's design technique is that the magnitude response of the all-pole lowpass filter cannot be controlled"* |
| Fettweis 1972（maximally flat delay 数字滤波器） | 仍是延迟参数化 |
| Zhang 2008（maximally flat IIR + flat group delay） | 按（群延迟 τ₀, 平坦度阶数 K/L, 阶数）列线性方程组，同样不是按截止频率 |

也就是说：**该族只有一个自由参数**（$\beta=2\delta$），它同时决定直流群延迟（$=\delta$）、
幅度响应（截止频率）和直流增益，三者不能各自指定。"指定截止频率"因此等于在这一个参数上解一个标量方程，
系数公式本身不变。

## 推导

数字全极点 Thiran 低通（含直流增益）：

$$H(z)=\frac{1}{A(z)},\quad
A(z)=\sum_{k=0}^{n}a_kz^{-k},\quad
a_k=(-1)^k\binom{n}{k}\prod_{i=0}^{n}\frac{\beta+i}{\beta+k+i},\quad a_0=1$$

$-3\mathrm{dB}$ 的定义（相对直流增益归一化）：

$$|H(e^{j\omega_c})|=\frac{|H(1)|}{\sqrt2}
\iff |A(e^{j\omega_c})|^2=2\,A(1)^2$$

展开成 $\beta$ 的方程（唯一的未知量）：

$$F(\beta,\omega_c)=\sum_{k,l}a_k(\beta)a_l(\beta)\cos\bigl((k-l)\omega_c\bigr)-2\Bigl(\sum_k a_k(\beta)\Bigr)^2=0$$

通分后是 $\beta$ 的 $2n$ 次多项式方程。令 $u=1-\cos\omega_c=2\sin^2(\omega_c/2)$：

- $n=1$： $u\beta^2+2u\beta-2=0$，闭式
  $$\boxed{\;\beta=-1+\sqrt{1+\frac{2}{u}}\;}$$
- $n=2$： $u^2\beta^4+8u^2\beta^3+u(6+19u)\beta^2+12u(2+u)\beta-36=0$（四次，取正根）
- 一般 $n$：$2n$ 次，无实用闭式 → 数值求根

**唯一性**（数值求根能成立的前提）：$\beta\to0^+$ 时 $A(z)\to1$，$F\to-1<0$；
$\beta\to\infty$ 时 $a_k\to(-1)^k\binom{n}{k}$、$A(1)\to0$，$F\to+\infty>0$。
实测 $n=1\ldots10$、25 个 $\omega_c$ 在 $\beta\in[10^{-4},10^{6}]$ 上 $F$ **最多一次符号变化**
（少数档位的根落在网格之外），$n\le5$ 时每次都是 1 次 ⇒ 根唯一。所以：从 0 向外扩括号 + Brent 必收敛。

## 判定

- 闭式（$n=1$、$n=2$ 四次根）与数值解一致到 $10^{-14}$（5 个截止频率各测）。
- 设计扫描（$n=1,2,3,4,6,8,10$ × $\omega_c\in[0.03\pi,0.97\pi]$ 共 133 组）：
  解出的 $\beta$ 代入后实测 $-3\mathrm{dB}$ 频率与目标最大差 $1.8\times10^{-15}$（$n\le4$）、
  $1.3\times10^{-9}$（$n=10$）。
- 精度上限由 $A(1)=\sum_k a_k$ 这个交替和的消位决定：条件数 $\sum|a_k|/|A(1)|$ 随 $\beta$（即随截止频率变小）迅速增大 ——
  $n=10$、$\omega_c=0.3\pi$ 时仅 $35$，$0.1\pi$ 时 $1.1\times10^{5}$，$0.05\pi$ 时 $7\times10^{7}$（可达精度 ≈ $\varepsilon\cdot$cond）。
  很窄的截止频率（大 $\beta$）在这一族 + 双精度下不可用。

## 代价：截止频率与增益绑定

同一个 $\beta$ 也定下直流增益 $1/A(1)$，不能各自指定：

| $n$ ＼ $\omega_c/\pi$ | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 | 0.9 |
|---|---|---|---|---|---|---|---|
| 1 | 6.89 | 3.74 | 2.19 | 1.71 | 1.37 | 1.25 | 1.21 |
| 4 | 461 | 43.6 | 6.38 | 2.85 | 1.52 | 1.24 | 1.16 |
| 10 | 2.2e5 | 961 | 18.6 | 4.57 | 1.74 | 1.27 | 1.13 |

（表里同时是直流群延迟 $\delta=\beta/2$：$\omega_c/\pi=0.1$、$n=10$ 时 $\delta=7.45$、增益 $961$。）

要把增益归一，只需把分子取 $A(1)$；但要再独立指定延迟/增益，就得离开这一族
（例如对模拟原型做频率缩放后再双线性变换，代价是数字侧多出 $n$ 个 $z=-1$ 零点）。

## 扫描图

`thiran_cutoff_scan.py` 用上面的设计器出两张图（也在 `output/`），每张两栏（幅度按直流归一化 + 群延迟）：

| 图 | 扫描 |
|---|---|
| `thiran_cutoff_scan_wc.png` | 固定 $n=6$，扫 $\omega_c/\pi=0.05,0.08,0.12,0.18,0.25,0.35,0.5,0.7,0.9$ |
| `thiran_cutoff_scan_n.png` | 固定 $\omega_c/\pi=0.2$，扫 $n=1,2,3,4,5,6,8,10,12$ |

- 幅度栏在目标 $\omega_c$ 处画点（实测 $-3\mathrm{dB}$）和虚线：18 组全部落在目标上（误差 $\le7\times10^{-15}$），
  可以直接看出"截止频率被钉住"这件事。
- 幅度栏的直流增益标注量化了代价：$n=6$ 时 $\omega_c/\pi=0.05$ 要 $+73.6\mathrm{dB}$、$0.9$ 只要 $+1.2\mathrm{dB}$；
  $\omega_c/\pi=0.2$ 时 $n=1\to12$ 从 $+6.8\mathrm{dB}$ 涨到 $+27.1\mathrm{dB}$。
- 群延迟栏的直流值就是 $\delta=\beta/2$（$\omega_c$ 越小、$n$ 越大，$\delta$ 越大），
  低频段一段平坦后掉到负值（该族没有零点，故不会出现尖峰）。

## 运行

```sh
python labs/thiran_filter/thiran_cutoff_design.py     # 打印闭式对照、唯一性、设计扫描、设计表
python labs/thiran_filter/thiran_cutoff_scan.py       # 两张扫描图 -> output/
```
