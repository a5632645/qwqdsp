# tf2ca — 数字 IIR 低通分解为两条全通链之和

奇数阶数字低通 `H(z) = B(z)/A(z)` 写成两条全通链之和：

```
H(z)  = 1/2 * ( A0(z) + A1(z) )        LP
Hc(z) = 1/2 * ( A0(z) - A1(z) )        HP，功率互补：|H|^2 + |Hc|^2 = 1
```

每条链是一串 1 阶 / 2 阶全通节点。这是 MATLAB `tf2ca` / `ca2tf` 的结构，也是
EMQF（Elliptic Minimal Q-Factor）滤波器“一个滤波器 + 一个加法”实现的基础。

## 结论

1. **奇数阶实系数可分解条件**：`n = deg A` 为奇数、`B` 镜像对称（`b[k] = b[n-k]`）、
   归一化 `H(1) = 1`。`butter` / `cheby1` / `cheby2` / `ellip` 的奇数阶低通满足后两条。
2. **极点分组**：把 `A` 的极点按 `|p|` 升序排列（1 个实极点 + `(n-1)/2` 个共轭对），
   逐个节点交替分给两条链；实极点所在链的阶数为奇数。对四族经典设计，这与按
   **s 平面极点角 `arg(s)`** 升序交替给出同一分区（允许两支互换）。
3. **数字域 `arg(p)` 不是等价规则**：`cheby2` 整族会出现与 `|p|` 不同的分区。
4. **精确判据**：`D0 = A / gcd(A, B+B_c)`、`D1 = A / gcd(A, B-B_c)`，其中
   `B_c^2 = B^2 - A*A^#` 是反镜像互补分子；`B_c(p) = ±B(p)` 决定极点归属。
5. **偶数阶**：实系数两路全通做不到（Nyquist 奇偶性）；放开复系数后
   `H = alpha*A1 + conj(alpha)*A2`，`D2 = conj(D1)`。四族极点互异的设计可归一成
   `H = 1/2(A1 + A2)`（两条复全通链）。对实输入只需算一条复链再取实部。
6. **重极点**：双链是否存在取决于 `B`。Linkwitz–Riley `LR_{2n} = (BW_n)^2` 不满足
   双链判据，正确结构是三条链（权重 `1/4, 1/2, 1/4`）。
7. **解析滤波器**：半带低通经频率旋转 `z -> -jz` 得到解析滤波器，落地形式是
   “两组 AP 链 + 延迟单元”。偶阶用复权重形式，奇阶用两组 `z^-2` 实全通链 + `z^-1`
   延迟；两者共用同一套半带设计与闭式 `Wn`。

## 数学

设 `A0 = N0/D0`、`A1 = N1/D1`，其中 `N_i = D_i^#`（系数反序，保证 `|N_i/D_i| ≡ 1`）。
代入并通分：

```
H = (N0*D1 + N1*D0) / (2*D0*D1)   =>   A = D0*D1,   B = (N0*D1 + N1*D0)/2
```

互补高通的分子为 `B_c = (N0*D1 - N1*D0)/2`，它是反镜像多项式（奇数阶故 `B_c(1) = 0`）。
两式相加、相减并相乘：

```
B + B_c = N0*D1,    B - B_c = N1*D0
B^2 - B_c^2 = N0*N1*D0*D1 = A*A^#      =>      B_c^2 = B^2 - A*A^#
```

### 极点归属

对 `A` 的每个单极点 `p`（`A(p) = 0`，故 `B_c(p)^2 = B(p)^2`）：

```
p ∈ D0  <=>  B_c(p) = +B(p)
p ∈ D1  <=>  B_c(p) = -B(p)
```

等价地 `D0 = A / gcd(A, B+B_c)`、`D1 = A / gcd(A, B-B_c)`。

### `B_c` 的递推

`B_c` 的系数满足 `c[k] = -c[n-k]`。记 `R = B^2 - A*A^#`，对 `k ≤ (n-1)/2` 的方程
`(c*c)[k] = R[k]` 是下三角的，可逐项递推：

```
c[0] = sqrt(R[0])
c[k] = ( R[k] - Σ_{i=1}^{k-1} c[i]*c[k-i] ) / (2*c[0])
```

其余系数由反镜像补出，`k > (n-1)/2` 的方程自动成立。

### 分组规则

对四族经典设计，正确的极点分区是 `|p|` 升序交替。两条独立依据：

- MATLAB `tf2ca` 的发布样例 `ellip(7,2,40,0.3)`：本目录得到
  `D0 = [1, -1.9549, 1.8352, -0.6971]`、`D1 = [1, -2.5164, 3.3197, -2.2151, 0.7469]`，
  与官方 `d0 = [1, -2.5163, 3.3183, -2.2130, 0.7457]`、`d1 = [1, -1.9550, 1.8366, -0.6983]`
  在四位小数内一致（两支互换）。
- EMQF 的 `apellip_du` 按 `beta = |p|^2` 排序后交替分配二阶节点，等价于按 `|p|` 交替。

#### `arg(s)` 与 `|p|` 等价

把数字极点经双线性反变换映回 s 平面（`s = 2*fs*(p-1)/(p+1)`，`fs = 2`；正实比例不改变角度），
按 `arg(s)` 升序交替，与按 `|p|` 升序交替给出同一分区。四族 × 奇数阶 3…71 ×
cutoff `1e-4…0.499999` × `rp 0.001…3 dB` × `rs 20…160 dB` 共 6020 例，分区分歧 0 例。
两者的排序裕度不同：`|p|` 的最小相邻间隔约 `8.3e-14`（`ellip` 51 阶 @ `1e-4`），
`arg(s)` 的最小间隔约 `1.3e-9` 度（`ellip` 41 阶 @ `0.3`）。

#### 数字域 `arg(p)` 不成立

`cheby2` 9 阶、cutoff `0.3` 的节点：

| 模长 | `arg(p)` | `arg(s)` | 模长链 | `arg(p)` 链 |
|---|---|---|---|---|
| 0.10043 | 0.000° | 180.000° | 0 | 0 |
| 0.35145 | 59.093° | 145.469° | 1 | 0 |
| 0.60934 | 54.920° | 122.227° | 0 | 1 |
| 0.79359 | 49.797° | 106.982° | 1 | 0 |
| 0.93400 | 47.104° | 95.329° | 0 | 1 |

`cheby2` 的极点角度挤在一起且随 `|p|` 非单调，按 `arg(p)` 交替会打乱分区：同一例子
`max|ΔH| = 1.004`，而 `|p|` 与 `arg(s)` 都在 `2e-13` 量级。对 2700 例做阳性对照，
`arg(p)` 失败 196 例（`butter` 0/675、`cheby1` 0/675、`cheby2` 188/675、`ellip` 8/675）。

#### 文献对应

Gazsi 1985 §2.2（*Explicit Formulas for Lattice Wave Digital Filters*）在 **ψ 平面**
（双线性变换后的参考滤波器复频率）上陈述“`g1` 与 `g2` 的零点在左半平面交替”，
并明确覆盖 inverse Chebyshev；ψ 与这里的 `s` 只差正实因子，角度相同。Vaidyanathan 1993
§3.6 把同一性质重述为 **z 平面**极点角（`θ0 < θ1 < …`，引用 Gazsi），该读法对
`cheby2` 整族不成立。工程实现使用 `|p|`（等价于 `arg(s)`）。

## 适用范围与数值边界

- `B_c^2 = B^2 - A*A^#` 是两个近似多项式的差，窄带 / 高阶时严重消位；`split_by_bc`
  的精度随之下降。四族经典设计用 `split_by_radius`（只依赖极点）。批量校验中，
  四族 × 奇数阶 3/5/7/9/11 × cutoff 0.1/0.2/0.35/0.5/0.7 的最差频响偏差为
  `3.8e-2`（`ellip` 11 阶 @ 0.1），其余组合大多低于 `1e-7`。
- 用 `np.roots` / `np.poly` 反复往返会在极点聚集时损失精度。需要落地到定点/C++ 时，
  应从设计公式直接生成节点系数，不要先求根再回代。
- **重极点**：`A` 有重根时，`split_by_radius` / `split_by_bc` 的“每个极点只属于一条链”
  前提失效；双链是否存在取决于 `B`，不能一概判为不存在。LR 是一类明确不可双链的设计，
  见 [even_order_repeated.md](even_order_repeated.md)。重极点下“逐极点 `alpha` 一致性”
  判据会假阳性，必须用系数空间残差判分。
- **偶数阶**：四族极点互异设计可复系数双链分解；窄带高阶时系数表示的条件数会把
  `B` 还原残差抬到 `1e-3` 量级（`cheby1` 8 阶 @ 0.1），这与奇数阶 `B_c` 的消位同源，
  不是分组错误。非四族设计未做系统验证。
- **解析滤波器**：偶阶半带的自由参数是 `(N, rs)`。闭式 `Wp` 本身满足边沿对称，但
  scipy `ellip` 在过渡带过窄时会失去半带不变量（`|H(0.5)|` 偏离 `1/sqrt(2)`、功率
  互补误差变大）。当前可靠范围：`N ≤ 16` 的常规规格，`N=20` 需要过渡带 `≥ 8e-3`，
  `N ≥ 24` 未找到可靠点。奇阶半带 `N ≤ 17` 时残差 `1e-12…1e-9`，`N ≥ 21` 同样受此
  限制。更高阶需要更高精度的椭圆设计。

## 主题文档

| 主题 | 文档 | 脚本 |
|---|---|---|
| 奇数阶实系数双链 | 本文 | `odd_order.py` |
| 偶数阶复系数双链（极点互异） | [even_order.md](even_order.md) | `even_order.py` |
| 偶数阶 + 重极点（LR 三链） | [even_order_repeated.md](even_order_repeated.md) | `even_order_repeated.py` |
| 偶阶半带 → 解析滤波器 | [hilbert_even.md](hilbert_even.md) | `hilbert_even.py` |
| 奇阶半带 → 解析滤波器 | [hilbert_odd.md](hilbert_odd.md) | `hilbert_odd.py` |

## 运行

```bash
# 奇数阶：默认 ellip 7 阶 @ 0.3；--check 跑四族 × 奇数阶 × 多截止频率
python qwqdsp/labs/tf2ca/odd_order.py
python qwqdsp/labs/tf2ca/odd_order.py --kind cheby2 --order 9 --cutoff 0.25 --stop 60
python qwqdsp/labs/tf2ca/odd_order.py --check

# 偶数阶复双链（极点互异）；--exact 用 60 位精度重算同一分组
python qwqdsp/labs/tf2ca/even_order.py --kind ellip --order 6 --cutoff 0.3
python qwqdsp/labs/tf2ca/even_order.py --check

# 偶数阶 + 重极点（LR）：双链不可，三链可
python qwqdsp/labs/tf2ca/even_order_repeated.py --half-order 2 --cutoff 0.25
python qwqdsp/labs/tf2ca/even_order_repeated.py --check

# 解析滤波器
python qwqdsp/labs/tf2ca/hilbert_even.py --order 16 --rs 60
python qwqdsp/labs/tf2ca/hilbert_even.py --header
python qwqdsp/labs/tf2ca/hilbert_odd.py --N 15 --rs 60
python qwqdsp/labs/tf2ca/hilbert_odd.py --check

# 分组规则：三种排序对照；--stress 跑 6020 例扫描
python qwqdsp/labs/tf2ca/split_rule.py
python qwqdsp/labs/tf2ca/split_rule.py --stress --controls
```

输出图片在 `output/`（本地生成，未纳入版本管理）。`odd_order.py` 会打印每条全通链的
1 阶/2 阶节点系数，可直接对照 `qwqdsp/include/qwqdsp/filter/allpass.hpp` 的
`AllpassOrder1` / `AllpassOrder2`。

## 文件

| 文件 | 说明 |
|---|---|
| `tf2ca.py` | 共享核心：节点化、模长分组、`B_c` 递推与精确判据、复全通响应、频响与误差 |
| `odd_order.py` | 奇数阶实双链：设计、分组、节点打印、频响校验、MATLAB 样例对照、出图 |
| `even_order.py` | 偶数阶复双链（极点互异）：共轭分组枚举、最小二乘、逐极点判据、单链形式、实系数对照 |
| `even_order.md` | 偶数阶复双链的推导、数据与边界 |
| `even_order_repeated.py` | 重极点 / LR：谱因子穷举、系数空间残差、三链构造与校验 |
| `even_order_repeated.md` | 重极点的双链判据、LR 的不可分解性与三链结构 |
| `hilbert_even.py` | 偶阶半带设计（闭式 `Wn`）、复全通分解、频率旋转、头文件系数评估 |
| `hilbert_even.md` | 偶阶解析滤波器的构造、归一化、质量与边界 |
| `hilbert_odd.py` | 奇阶半带：两组 `z^-2` 实全通链 + 延迟，旋转后对应头文件结构 |
| `hilbert_odd.md` | 奇阶解析滤波器的度数论证、延迟极点、DC 相位奇偶律 |
| `split_rule.py` | 分组规则对照与压力测试（模长 / `arg(s)` / `arg(p)`） |
| `output/` | 生成的图片（本地，未纳入版本管理） |

## 参考

- MATLAB `tf2ca` / `ca2tf`（Coupled-Allpass decomposition）
- L. Gazsi, "Explicit Formulas for Lattice Wave Digital Filters", IEEE Trans. Circuits and
  Systems, vol. 32, no. 1, pp. 68–88, Jan. 1985.
- P. P. Vaidyanathan, *Multirate Systems and Filter Banks*, Prentice Hall, 1993, §3.6.
- Lj. D. Milic, M. D. Lutovac, "Efficient Algorithm for the Design of High-Speed Elliptic
  IIR Filters", AEU 57(4), 2003（EMQF / `apellip_du`）
- https://github.com/vadkudr/EMQFfilters
- f. harris, "A Most Efficient Digital Filter: The Two-Path Recursive All-Pass Filter"
- R. Lyons, "Reducing IIR Filter Computational Workload"
