# tf2ca — 数字 IIR 低通分解为两条全通链之和

本目录探索：把一个**奇数阶**数字 IIR 低通滤波器 `H(z) = B(z)/A(z)` 写成两条
全通滤波器链之和

```
H(z)  = 1/2 * ( A0(z) + A1(z) )        LP
Hc(z) = 1/2 * ( A0(z) - A1(z) )        HP（功率互补：|H|^2 + |Hc|^2 = 1）
```

其中每条链都是一串 1 阶 / 2 阶全通节点。这就是 MATLAB `tf2ca` / `ca2tf` 做的事，
也是 EMQF（Elliptic Minimal Q-Factor）滤波器「一个滤波器 + 一个加法」实现的基础。

## 结论速览

- **可分解条件**：阶数 `n` 为奇数 + 分子 `B` 为镜像对称（`b[k] == b[n-k]`）+ 归一化 `H(1)=1`。
  `butter` / `cheby1` / `cheby2` / `ellip` 天然满足后两条，所以**这四族奇数阶低通都能分解**。
- **极点分组规则**（四族经典设计的工程规则）：把 `A` 的所有极点按**模长 |p| 升序**排列
  （1 个实极点 + `(n-1)/2` 个共轭对），**逐个节点交替**分到两条链。
  实极点所在链的阶数为奇数。
- **不是按幅角交替**：对 `cheby2` 有反例（见下文），网上/ Lyons 文中「pole interlacing
  property / 按角度交替」的说法在四族统一意义下不成立，按模长才统一成立。
- **通用精确判据**：`D0 = A / gcd(A, B+B_c)`、`D1 = A / gcd(A, B-B_c)`，
  其中 `B_c = sqrt(B^2 - A*A^#)` 是**反镜像**分子（见下文推导）。四族经典设计与
  上面的模长规则完全一致。

## 运行

```bash
# 默认 ellip 7 阶、cutoff=0.3（与 EMQF 文档里 MATLAB tf2ca 的例子同参）
python qwqdsp/labs/tf2ca/run.py

# 指定滤波器
python qwqdsp/labs/tf2ca/run.py --kind butter --order 7 --cutoff 0.3
python qwqdsp/labs/tf2ca/run.py --kind cheby2 --order 9 --cutoff 0.25 --stop 60

# 批量校验表（四族 × 3/5/7/9/11 阶 × 多截止频率）
python qwqdsp/labs/tf2ca/run.py --check
```

输出图片在 `output/`（本地生成）。`run.py` 会打印每条全通链的 1 阶/2 阶节点系数，
可直接对照 `include/qwqdsp/filter/allpass.hpp` 里的 `AllpassOrder1` / `AllpassOrder2` 实现。

## 数学推导

设两条链为全通：`A0 = N0/D0`、`A1 = N1/D1`，其中 `N_i` 是 `D_i` 的**系数反序**
（即镜像多项式 `N_i(z) = z^{-d_i} D_i(1/z)`，保证 |N/D| ≡ 1）。代入并通分：

```
H = 1/2 (N0/D0 + N1/D1) = (N0*D1 + N1*D0) / (2*D0*D1)
```

于是必须有

```
A = D0 * D1
B = (N0*D1 + N1*D0) / 2
```

`B` 是「镜像多项式的对称部分」，因此**必然是镜像对称**的。反过来，互补高通分子

```
B_c = (N0*D1 - N1*D0) / 2
```

是**反镜像**的（`c[k] == -c[n-k]`，奇数阶故 `B_c(1)=0`）。两式相加/相减得

```
B + B_c = N0*D1        B - B_c = N1*D0
```

相乘并注意 `N0*N1 = D0^# * D1^# = A^#`（`A^#` 为 A 的镜像）：

```
B^2 - B_c^2 = A * A^#        =>        B_c^2 = B^2 - A*A^#
```

### 极点归哪条链

对 `A` 的每个（单）极点 `p`：`p ∈ D0 ⟺ (B+B_c)(p)=0`... 反过来由

`B + B_c = N0*D1`、`B - B_c = N1*D0`

可直接得到判据

```
p ∈ D0  <=>  B_c(p) = +B(p)
p ∈ D1  <=>  B_c(p) = -B(p)
```

（因为 `B_c² = B² - A·A^#`，而 `A(p)=0`，所以 `B_c(p) = ±B(p)` 恒成立。）
等价地 `D0 = A / gcd(A, B+B_c)`、`D1 = A / gcd(A, B-B_c)`。

### B_c 的快速求法（无需开根/找根）

`B_c` 系数满足反镜像 `c[k] = -c[n-k]`；记 `R = B^2 - A*A^#`。对 `k ≤ (n-1)/2`，
方程 `(c*c)[k] = R[k]` 是下三角的，可逐项递推：

```
c[0] = sqrt(R[0])
c[k] = ( R[k] - Σ_{i=1}^{k-1} c[i]*c[k-i] ) / (2*c[0])
```

其余系数由反镜像补出，且 `k > (n-1)/2` 的方程自动成立。见 `tf2ca.bc_antisymmetric_sqrt`。

## 为什么工程上「按模长交替」能用

对 `butter/cheby1/cheby2/ellip` 这四族，正确的极点分组恰好等于「按 |p| 升序交替」。
两组独立证据：

1. **MATLAB `tf2ca` 发布样例**（EMQF 文档，`ellip(7,2,40,0.3)`）：
   本目录 `validate_reference()` 复现出与官方完全相同的结果（连 `tf2sos` 的分段都一致）：
   ```
   官方 d0 = [1 -2.5163 3.3183 -2.2130 0.7457]   d1 = [1 -1.9550 1.8366 -0.6983]
   本目录   D1 = [1 -2.5164 3.3197 -2.2151 0.7469]  D0 = [1 -1.9549 1.8352 -0.6971]
   ```
2. **EMQF `apellip_du.m`**：它按 `beta = |pole|^2` 排序后**交替**把二阶节点分给 `p0`/`p1`，
   与「按 |p| 交替」等价。

### 反例：不能按幅角交替

`cheby2` 9 阶、cutoff=0.3 时，唯一正确的分组是（按极点幅角列出节点）

```
幅角:      0.00   0.82   0.87   0.96   1.03
归属:      D0     D0     D1     D0     D1      <- 不是交替
```

而按 |p| 排序后归属恰好交替。类似地 `cheby2` 7 阶、`ellip(5,0.7)` 也都不满足幅角交替。
所以 Lyons 文中「pole interlacing property（按角度交错）」只对部分情形成立，
**统一规则是模长**。

## 与仓库现有代码的关系

`include/qwqdsp/filter/parallel_allpass.hpp`：

- 四族（butter/cheby1/cheby2/ellip）都在此实现，系数由「奇数阶模拟原型 + 双线性变换」
  直接生成（`detail::*Prototype`），极点按 |p| 升序交替分配到两条链，**不求根**。
- 该文件原先的 `BuildChebyshev1` 有符号错误（实极点用了 `+sinh(A)`，应取 `-sinh(A)`），
  所有极点被映射到单位圆外；因为 `N = 镜像(D)` 使 `|N/D|` 在单位圆上恒为 1，所以
  幅度谱正确但系统不稳定——与旧注释「幅度谱对但是不稳定」吻合。**此错误已修复**。
- 数值验证见 `tests/filter/parallel_allpass.cpp`（全通平坦度、功率互补、
  与各族解析幅频逐点对比、椭圆阻带边沿与传输零点个数、退化规格返回 false）。

## 数值注意

`B_c^2 = B^2 - A*A^#` 是两个近似多项式的差，**窄带 / 高阶**时会严重消位（条件数问题），
`split_by_bc` 的精度随之下降。四族经典设计建议直接用 `split_by_radius`（只依赖极点，
自 `--check` 表可见绝大多数组合频响偏差 < 1e-7）。另外用 `np.roots`/`np.poly` 反复
往返也会在极点聚集时损失精度；真正要落地（如 C++）时应像 EMQF 那样**直接由设计公式
生成节点系数**，而不是先求根再回代。

## 文件

| 文件 | 说明 |
|------|------|
| `tf2ca.py` | 核心：极点分组(`split_by_radius`/`split_by_bc`)、`B_c` 递推、还原、频响校验、节点拆分 |
| `run.py`   | 驱动：设计→分解→校验→打印节点→对照 MATLAB 样例→出图；`--check` 批量校验 |
| `output/`  | 生成的图片（本地，未纳入版本管理） |

## 参考

- MATLAB `tf2ca` / `ca2tf`（Coupled-Allpass decomposition）
- Lj. D. Milic, M. D. Lutovac, *Efficient Algorithm for the Design of High-Speed Elliptic IIR Filters*, AEU 57(4), 2003 —— EMQF / `apellip_du`
- https://github.com/vadkudr/EMQFfilters 及 http://vadkudr.org/Algorithms/EMQFdemo/EMQFdemo.html
- f. harris, *A Most Efficient Digital Filter: The Two-Path Recursive All-Pass Filter*（researchgate 278320928）
- R. Lyons, *Reducing IIR Filter Computational Workload*（仓库 `parallel_allpass.hpp` 引用之一）
