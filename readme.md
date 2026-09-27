
# qwqdsp

存储我曾经使用过的东西（重复造轮子！），为了防止每次都要去仓库找代码复制，将它们全部打包到一起

> [!NOTE]
> 这不是什么正经的DSP库，里面存储了大量未优化和尝试的代码

> 贡献  
> 很乐意有人提交pull request（如果有的话就好）  

---

## 📦 依赖 (Dependencies)

- **C++20** 或更高版本。

| 依赖 | 开关 | 说明 |
|------|------|------|
| **Eigen3** | `QWQDSP_HAVE_EIGEN`（由外部定义） | 启用 `rls_filter`、`swift_f0` 等（外部提供 Eigen 并负责 include） |
| **Intel IPP** | `QWQDSP_HAVE_IPP`（由外部定义） | 替换 Ooura FFT 为 IPP 后端（外部提供 IPP 并链接） |
| **Apple Accelerate** | `QWQDSP_HAVE_ACCELERATE`（由外部定义） | macOS 上替换 Ooura FFT 为 vDSP 后端（外部链接 Accelerate） |
| **SIMDe** | `QWQDSP_HAVE_SIMDE`（由外部定义） | 非 x86 平台模拟 SIMD（头文件直接 include `<x86/avx2.h>` / `<x86/sse4.1.h>`） |
| **raylib** | `QWQDSP_USE_RAYLIB=ON` | 构建 raylib 依赖的 GUI example（`example/gui`）与 `labs/playing` 实验；raylib 由外部提供 |

---
