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
- **按"角度"交替要看是哪个角度**：**数字域** `arg(p)` 对 `cheby2` 有反例（见下文）；
  但把极点映回 **s 平面**后按 `arg(s)` 交替，在测试范围内与模长规则完全等价
  （**仅为数值测试支持的猜想，没有文献证据**，见下文「猜想」一节）。
- **统一规则是模长 |p|**：按 |p| 升序交替对四族经典设计统一成立。
- **通用精确判据**：`D0 = A / gcd(A, B+B_c)`、`D1 = A / gcd(A, B-B_c)`，
  其中 `B_c = sqrt(B^2 - A*A^#)` 是**反镜像**分子（见下文推导）。四族经典设计与
  上面的模长规则完全一致。
- **解析滤波器（Hilbert）**：把半带低通做**频率旋转** `z -> -jz` 就得到解析滤波器，落地形式正是
  `qwqdsp/include/qwqdsp/filter/iir_hilbert.hpp` 的「两组 AP 链 + 延迟单元」。
  偶阶走复权重 `alpha` 形式（`hilbert_even.py`/`.md`），奇阶走"两组 `z^-2` 实全通链 + 延迟"
  （`hilbert_odd.py`/`.md`）。两者共用本目录的同一套半带设计与闭式 `Wn`。
- **偶数阶**：实系数两路全通做不到（Nyquist 奇偶性），但放开**复系数全通 + 复权重**后
  偶数阶同样能分解（`H = alpha*A1 + conj(alpha)*A2`，只算一条复链再取实部）。
  推导、数据与未验证边界见 **[even_order.md](even_order.md)**，脚本 `even_order.py`。

## 运行

```bash
# 默认 ellip 7 阶、cutoff=0.3（与 EMQF 文档里 MATLAB tf2ca 的例子同参）
python qwqdsp/labs/tf2ca/odd_order.py

# 指定滤波器
python qwqdsp/labs/tf2ca/odd_order.py --kind butter --order 7 --cutoff 0.3
python qwqdsp/labs/tf2ca/odd_order.py --kind cheby2 --order 9 --cutoff 0.25 --stop 60

# 批量校验表（四族 × 3/5/7/9/11 阶 × 多截止频率）
python qwqdsp/labs/tf2ca/odd_order.py --check
```

（`odd_order.py` 只跑**奇数阶**；偶数阶的**复系数**分解见 `even_order.py`，
结论摘要见 [even_order.md](even_order.md)。）

输出图片在 `output/`（本地生成）。`odd_order.py` 会打印每条全通链的 1 阶/2 阶节点系数，
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
所以 Lyons 文中「pole interlacing property（按角度交错）」**在数字域角度下**只对部分情形成立；
统一规则是模长。（s 平面角另见下节。）

### 猜想：文献里的"角度"可能指 s 平面上的极点角

把数字极点经**双线性反变换**映回 s 平面（`s = 2*fs*(p-1)/(p+1)`，本目录取 `fs = 2`），
再按 `arg(s)` 升序交替分配节点，与本库「按 |p| 升序交替」在测试范围内**给出完全相同的分区**：

| 脚本 | 覆盖 | 出问题 |
|---|---|---|
| `arg_s_plane.py` | 4 族 × 奇数阶 3/5/7/9/11 × cutoff 0.1…0.5（100 例） | **0/100**（分区与 `\|p\|` 全同，残差 ~1e-13） |
| `arg_s_stress.py` A–E | 4 族 × 奇数阶 3…71 × cutoff 1e-4…0.499999（E 组再细扫 400 个点）× rp 0.001–3 dB × rs 20–160 dB（6020 例） | **0/6020** |

同一套件下**两条规则都没失败过**（分区 0 分歧），差别只在**排序裕度**：`|p|` 相邻间隔最小到
**8.2e-14**（绝对量，A 组 ellip 31/0.0001），`arg(s)` 最小到 **1.3e-9 度**——真要有规则先崩，
先崩的会是 `|p|`。

同一套件的阳性对照：数字域 `arg(p)` 失败 196/2700，`arg(p-1)`（丢掉 `arg(p+1)` 项的残缺写法）
失败 245/2700；两者失效区互不重叠，都是 O(1) 残差量级，而 `arg(s)` 在排序间隔小到 1.3e-9 度的
用例里依然正确。这说明「按角度交替」的正确读法很可能是 **s 平面角**，而 `|p|` 规则是它的等价形式。

> ⚠️ **这只是数值测试支持的猜想**：我们**没有检索或核对文献原文**，"文献里说的按角度交替 =
> s 平面极点角"这一步**没有文献证据**。若要写进正式说明，需要先找到原文出处并逐字确认。

### 复测：`arg(p)` 交替在 `cheby2` 上确实失败（2026-10-03 重跑）

四族 × 奇数阶 3/5/7/9/11 × cutoff 0.1/0.2/0.35/0.5（共 80 例），比较两种分法：
**按 `|p|` 升序交替** vs **按 `arg(p)` 升序交替**（角度取共轭对上半平面根的 `|arg|`）。

| 规则 | 失败例数（残差 > 1e-8） | 备注 |
|---|---|---|
| 按 `\|p\|` 交替 | **16 / 80** | 失败者都是窄带高阶（复现实验的固有条件数限制，见下） |
| 按 `arg(p)` 交替 | **20 / 80** | 失败集合 = 上表 16 例 **∪ 4 个 cheby2 例** |

> 注：本表数字来自早期的多项式管线（`np.roots` + `B/A` 多项式求值）。改用 scipy 的 zpk
> 直接取极点、`freqz_zpk` 求响应后，同一批例子里两种规则的残差都回到 1e-13 量级——
> "16/80 条件数失败"是管线误差，不是规则问题（见 `arg_s_plane.py` / `arg_s_stress.py`）。

那 4 个额外失败（`|p|` 规则在同样例子上是机器精度）：

| 例 | `\|p\|` 残差 | `arg` 残差 |
|---|---|---|
| cheby2 n=7 cutoff=0.1 | 2.2e-10 | **6.5e-02** |
| cheby2 n=7 cutoff=0.2 | 5.6e-12 | **2.3e-01** |
| cheby2 n=9 cutoff=0.2 | 6.4e-12 | **2.9e-01** |
| cheby2 n=9 cutoff=0.35 | 9.3e-13 | **1.7e+00** |

失败的原因可以直接看出来 —— `cheby2` 的极点**角度几乎挤在一起**，按角度排序后"交替"不再等价于正确分组：

```
cheby2 n=9, cutoff=0.2 的节点（按 |p| 排序）：
  节点0: |p|=0.3147  arg=  0.000°
  节点1: |p|=0.4820  arg= 33.221°
  节点2: |p|=0.7002  arg= 35.264°
  节点3: |p|=0.8485  arg= 32.701°
  节点4: |p|=0.9531  arg= 31.046°     <- 角度不是单调的，按角度交替会被打乱
```

所以本目录的结论（**统一规则是按 `|p|`**；Lyons 文中"pole interlacing（按角度交错）"只对部分情形成立）
经复测仍然成立，不建议改成按角度交替。

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
| `odd_order.py` | 奇数阶（实系数）：设计→分解→校验→打印节点→对照 MATLAB 样例→出图；`--check` 批量校验 |
| `even_order.py` | 偶数阶（复系数复权重）：共轭分组枚举、最小二乘、逐极点判据、单链形式、实系数对照 |
| `even_order.md` | 偶数阶分解的研究总结（结论 / 推导 / 数据 / 未验证边界） |
| `hilbert_even.py` | 偶阶「解析滤波器 / Hilbert」构造：半带 → 复全通分解 → 旋转；`--header` 给头文件固定系数打分 |
| `hilbert_even.md` | 上面这条路线的研究总结（含"让数字极点在虚轴上"的三条设计条件与 `Wn` 闭式解） |
| `hilbert_odd.py` | 奇阶「解析滤波器 / Hilbert」构造：半带 → 两组 `z^-2` AP 链 + 延迟 → 旋转（= `iir_hilbert.hpp` 的形式） |
| `hilbert_odd.md` | 奇阶路线的研究总结（含"延迟单元 = 链 A1 上 `z=0` 的极点"的极点计数论证） |
| `audio_developer_guide.md` | 面向普通音频开发者的并行全通与解析滤波器总结与实战避坑指南 |
| `output/`  | 生成的图片（本地，未纳入版本管理） |

## 参考

- MATLAB `tf2ca` / `ca2tf`（Coupled-Allpass decomposition）
- Lj. D. Milic, M. D. Lutovac, *Efficient Algorithm for the Design of High-Speed Elliptic IIR Filters*, AEU 57(4), 2003 —— EMQF / `apellip_du`
- https://github.com/vadkudr/EMQFfilters 及 http://vadkudr.org/Algorithms/EMQFdemo/EMQFdemo.html
- f. harris, *A Most Efficient Digital Filter: The Two-Path Recursive All-Pass Filter*（researchgate 278320928）
- R. Lyons, *Reducing IIR Filter Computational Workload*（仓库 `parallel_allpass.hpp` 引用之一）
