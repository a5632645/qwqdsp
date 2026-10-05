#pragma once
#include <array>
#include <complex>
#include <cstddef>
#include <span>

namespace qwqdsp_filter {
// ------------------------------------------------------------
// 半带全通和多项式上/下采样器
//
// 结构来源：奇数阶低通可写成 H = 1/2[A0(z) + z^-1 A1(z)]（两条全通链，极点按 |p| 升序逐节交替分组）；
// 半带的极点成 ±jρ 对，于是两条链都只是 ζ = z^-2 的函数 —— **分支可以跑在半速率**，
// 这是 IIR 也能「多相」的原因（直接型的伴随矩阵 A^L 在高阶下病态，实测误差量级 1e-1）。
//
// 2^K 倍超采样 = K 级级联半带（每级 2×）。这里三级复用同一组系数：
// 限制通带的只有第一级（它的半带边沿 = 该级速率 fs/4 = 输入采样率/2），
// 后两级的交叉点（48/96 kHz）远高于信号带，只压各自补零产生的镜像。
// 实测三级复合响应：通带 [0, 0.45·fs_in] 起伏 < -190 dB，阻带 [>= 0.527·fs_in] 最坏 -100 dB。
//
// 解析（单边）路径 = 半带插值 ∘ 解析全通 H_an = 1/2[A0(-ζ) + j·z^-1·A1(-ζ)]
// （频率旋转 z -> -jz 在 ζ 上的作用；延迟支路接**虚部**槽）。
// 注意解析级跑在最高速率上，其单边性只在 |f| 大于约 (1-2·Wp)·fs/2（@384 kHz 约 5.25 kHz）时成立。
// ------------------------------------------------------------

// ----- 椭圆半带 N=19 / 阻带 100 dB / Wp=0.4726712 的两条 ζ 全通链节参数 -----
// 节为 (ζ^-1 + a)/(1 + a·ζ^-1)，ζ^-1 = z^-2；a = ρ²。
inline constexpr std::array<float, 5> kHalfbandN19Chain0{
    3.993399361e-02f, 2.953322197e-01f, 5.952248425e-01f, 8.146834120e-01f, 9.656342782e-01f,
};
inline constexpr std::array<float, 4> kHalfbandN19Chain1{
    1.478918562e-01f, 4.516440986e-01f, 7.163745215e-01f, 8.951804120e-01f,
};
inline constexpr std::size_t kHalfbandChain0Len = kHalfbandN19Chain0.size();
inline constexpr std::size_t kHalfbandChain1Len = kHalfbandN19Chain1.size();

namespace detail {
/**
 * @brief 低速率一阶全通链（延迟 z^-2 = ζ^-1）
 * @tparam T 采样类型
 * @tparam N 节数
 *
 * `rotated` 表示 A(−ζ)（解析/单边路径）：节 H = (ζ^-1 + a)/(1 + a·ζ^-1)，ζ → −ζ 后等价递推
 * `out = a·in − s`（未旋转为 `out = a·in + s`）。两式都对应同一个 `in = y ∓ a·s`。
 */
template <class T, std::size_t N>
class ZetaAllpassChain {
public:
    void Init(const std::array<T, N>& a, bool rotated) noexcept {
        a_ = a;
        rotated_ = rotated;
        s_.fill(T{});
    }

    void Reset() noexcept { s_.fill(T{}); }

    [[nodiscard]] T Tick(T x) noexcept {
        T y = x;
        for (std::size_t i = 0; i < N; ++i) {
            const T a = a_[i];
            const T s = s_[i];
            if (rotated_) {
                const T in = y + a * s;
                y = a * in - s;
                s_[i] = in;
            }
            else {
                const T in = y - a * s;
                y = a * in + s;
                s_[i] = in;
            }
        }
        return y;
    }

private:
    std::array<T, N> a_{};
    std::array<T, N> s_{};
    bool rotated_{};
};

/// 一级 2× 的公共部件：两条链 + 「偶相位取 A0、奇相位取 A1」
template <class T>
class HalfbandStage {
public:
    void Reset() noexcept {
        a0_.Reset();
        a1_.Reset();
    }

    /// 一个输入样本 → 两个输出样本（上采样用）
    void Tick2(T x, T* out) noexcept {
        out[0] = a0_.Tick(x);
        out[1] = a1_.Tick(x);
    }

    /// 分别取两条链（抽取用：偶相位走 A0、奇相位走 A1）
    [[nodiscard]] T TickA0(T x) noexcept { return a0_.Tick(x); }
    [[nodiscard]] T TickA1(T x) noexcept { return a1_.Tick(x); }

    /// 初始化两条链（rotated = 解析路径的 A(−ζ)）
    void Init(bool rotated) noexcept {
        a0_.Init(kHalfbandN19Chain0, rotated);
        a1_.Init(kHalfbandN19Chain1, rotated);
    }

private:
    ZetaAllpassChain<T, kHalfbandChain0Len> a0_{};
    ZetaAllpassChain<T, kHalfbandChain1Len> a1_{};
};
} // namespace detail

/**
 * @brief 2^K 倍半带全通和多相**上采样**（实数）
 * @tparam T 采样类型
 * @tparam K 级数（输出倍率 = 2^K）
 *
 * 每级的两条链都跑在该级输入速率上（级 k 处理点数 = 2^k）。
 */
template <class T = float, std::size_t K = 3>
class HalfbandAllpassUpsampler {
public:
    static constexpr std::size_t kFactor = std::size_t{1} << K;

    HalfbandAllpassUpsampler() noexcept {
        for (auto& s : stages_) {
            s.Init(false);
        }
    }

    void Reset() noexcept {
        for (auto& s : stages_) {
            s.Reset();
        }
    }

    /**
     * @brief 推进一个输入样本，写出 kFactor 个上采样样本
     * @param x 输入样本
     * @param out 输出缓冲，恰好 kFactor 个
     */
    void Tick(T x, std::span<T, kFactor> out) noexcept {
        // 必须**按时间顺序（升序）**处理：链是有状态的，倒序会把输出分到错误相位上。
        // 乒乓缓冲，避免原地展开时覆盖尚未读取的样本。
        T* src = ping_.data();
        T* dst = pong_.data();
        src[0] = x;
        std::size_t n = 1;
        for (std::size_t k = 0; k < K; ++k) {
            for (std::size_t i = 0; i < n; ++i) {
                T pair[2];
                stages_[k].Tick2(src[i], pair);
                dst[2 * i] = pair[0];
                dst[2 * i + 1] = pair[1];
            }
            T* tmp = src;
            src = dst;
            dst = tmp;
            n *= 2;
        }
        for (std::size_t i = 0; i < kFactor; ++i) {
            out[i] = src[i];
        }
    }

private:
    /// 展开用的临时缓冲（最大 2^K 个）
    std::array<T, kFactor> ping_{};
    std::array<T, kFactor> pong_{};
    std::array<detail::HalfbandStage<T>, K> stages_{};
};

/**
 * @brief 2^K 倍半带全通和多相**解析（单边）上采样**：输出复数样本
 * @tparam T 采样类型
 * @tparam K 级数
 *
 * 用同一个半带全通和做到带限，再接解析全通（H = 1/2[A0(−ζ) + j·z^-1·A1(−ζ)]）得到单边信号。
 * 两支路全通和只能给出「对称半带」或「平坦单边」两种形态，做不出「单边且带限」，
 * 所以必须**级联**：先半带抗镜像，再解析化。
 */
template <class T = float, std::size_t K = 3>
class HalfbandAllpassAnalyticUpsampler {
public:
    static constexpr std::size_t kFactor = std::size_t{1} << K;

    HalfbandAllpassAnalyticUpsampler() noexcept {
        for (auto& s : stages_) {
            s.Init(false);
        }
        ve_re_.Init(kHalfbandN19Chain0, true); ///< 偶相位 → 解析信号实部
        vo_re_.Init(kHalfbandN19Chain0, true); ///< 奇相位 → 实部（上一拍）
        ve_im_.Init(kHalfbandN19Chain1, true); ///< 偶相位 → 虚部（下一拍用）
        vo_im_.Init(kHalfbandN19Chain1, true); ///< 奇相位 → 虚部（上一拍）
    }

    void Reset() noexcept {
        for (auto& s : stages_) {
            s.Reset();
        }
        ve_re_.Reset();
        vo_re_.Reset();
        ve_im_.Reset();
        vo_im_.Reset();
        im_prev_ = T{};
    }

    /**
     * @brief 推进一个输入样本，写出 kFactor 个解析样本
     * @param x 输入样本（实数）
     * @param out 输出缓冲，恰好 kFactor 个复数
     */
    void Tick(T x, std::span<std::complex<T>, kFactor> out) noexcept {
        // 同实数版：升序 + 乒乓（链有状态）
        T* src = ping_.data();
        T* dst = pong_.data();
        src[0] = x;
        std::size_t n = 1;
        for (std::size_t k = 0; k < K; ++k) {
            for (std::size_t i = 0; i < n; ++i) {
                T pair[2];
                stages_[k].Tick2(src[i], pair);
                dst[2 * i] = pair[0];
                dst[2 * i + 1] = pair[1];
            }
            T* tmp = src;
            src = dst;
            dst = tmp;
            n *= 2;
        }
        const T* up = src;
        // 解析级：偶/奇各半，实部来自 A0、虚部来自 A1（延迟一个偶相位样本接奇相位）
        for (std::size_t m = 0; m < kFactor / 2; ++m) {
            const T ve = up[2 * m];
            const T vo = up[2 * m + 1];
            const T re_e = ve_re_.Tick(ve);
            const T re_o = vo_re_.Tick(vo);
            const T im_e = ve_im_.Tick(ve);
            const T im_o = vo_im_.Tick(vo);
            out[2 * m] = std::complex<T>{re_e, im_prev_};
            out[2 * m + 1] = std::complex<T>{re_o, im_e};
            im_prev_ = im_o;
        }
    }

private:
    std::array<T, kFactor> ping_{};
    std::array<T, kFactor> pong_{};
    std::array<detail::HalfbandStage<T>, K> stages_{};
    detail::ZetaAllpassChain<T, kHalfbandChain0Len> ve_re_{};
    detail::ZetaAllpassChain<T, kHalfbandChain0Len> vo_re_{};
    detail::ZetaAllpassChain<T, kHalfbandChain1Len> ve_im_{};
    detail::ZetaAllpassChain<T, kHalfbandChain1Len> vo_im_{};
    T im_prev_{};
};

/**
 * @brief 2^K 倍半带全通和多相**抽取**（实数抗混叠）
 * @tparam T 采样类型
 * @tparam K 级数
 *
 * 与上采样同一组系数：每级 `y = 1/2[A0(偶相位) + A1(奇相位)延迟一拍]`。
 */
template <class T = float, std::size_t K = 3>
class HalfbandAllpassDecimator {
public:
    static constexpr std::size_t kFactor = std::size_t{1} << K;

    HalfbandAllpassDecimator() noexcept {
        for (auto& s : stages_) {
            s.Init(false);
        }
    }

    void Reset() noexcept {
        for (auto& s : stages_) {
            s.Reset();
        }
        prev_.fill(T{});
    }

    /**
     * @brief 吃 kFactor 个输入样本，吐 1 个
     * @param in 输入缓冲，恰好 kFactor 个
     * @return 抽取后的样本
     */
    [[nodiscard]] T Tick(std::span<const T, kFactor> in) noexcept {
        for (std::size_t i = 0; i < kFactor; ++i) {
            buffer_[i] = in[i];
        }
        std::size_t n = kFactor;
        for (std::size_t k = 0; k < K; ++k) {
            const std::size_t half = n / 2;
            for (std::size_t m = 0; m < half; ++m) {
                const T ce = stages_[k].TickA0(buffer_[2 * m]);
                const T co = stages_[k].TickA1(buffer_[2 * m + 1]);
                buffer_[m] = T{0.5} * (ce + prev_[k]);
                prev_[k] = co;
            }
            n = half;
        }
        return buffer_[0];
    }

private:
    std::array<T, kFactor> buffer_{};
    std::array<T, kFactor> prev_{};
    std::array<detail::HalfbandStage<T>, K> stages_{};
};
} // namespace qwqdsp_filter
