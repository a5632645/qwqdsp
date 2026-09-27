#pragma once

#ifdef QWQDSP_HAVE_EIGEN

#include <Eigen/Dense>
#include <cmath>

namespace qwqdsp_adaptive {
/**
 * @brief 递推最小二乘（RLS）自适应滤波器，带指数遗忘因子
 *
 * 每步用回归向量 `φ = [source, x[n-1], x[n-2], …, x[n-ORDER+1]]` 去预测 `target`，
 * 以遗忘因子 λ 最小化加权平方误差 `Σ λ^(n-i)·(target_i − wᵀφ_i)²`，即标准 RLS：
 *
 *     k = P φ / (λ + φᵀ P φ)
 *     w ← w + k · (target − wᵀφ)
 *     P ← (I − k φᵀ) P / λ
 *
 * 采样顺序：`Tick` **先把 source 移进延迟线再做预测**，所以 `φ[0]` 就是本次的 source。
 *
 * 两种典型用法：
 * - **一步预测 / LPC**：`Tick(x[n], x[n+1])`——source 与 target 取同一段信号、错开一个采样。
 *   收敛后 `GetCoeff()` 拿到的 w 就是预测系数 `x̂[m] = Σ_j w[j-1]·x[m-j]`，对应的全极点模型
 *   `H(z) = 1 / (1 − Σ_j w[j-1]·z^-j)` 与 `Filter()` 的内部结构是同一个滤波器。
 * - **自适应辨识**：source 给参考输入、target 给期望输出，收敛后的 w 即被辨识系统的冲激响应。
 *
 * @warning source 与 target 给**同一个值**是退化用法：最小二乘存在零误差解 `w = e0`
 *          （把当前样本直接抄成输出），RLS 必然收敛到它——实测残差 −60 dB、w0≈0.96、‖w[1:]‖≈0.12，
 *          系数里不含任何信号谱信息。要得到全极点模型，target 必须给"后一个样本"。
 *
 * @warning P0（`Init` 的 `identity` 参数）是**初始逆相关矩阵**，越大表示"w=0 的先验"越弱、收敛越快。
 *          类默认 0.01 与信号量级无关，对常见量级偏小到几乎学不动：order 20、λ=1.0 在 RMS≈0.11 的
 *          语音段上跑 1.4 万点，尾部残差只有 −18.7 dB，而同一套回归量/目标直接解法方程是 −49.3 dB；
 *          传 `Init(1e6)` 后两者才对得上（−48.9 dB）。
 *          qwqfixme: 默认值应改成按数据量级归一（如 P0·E[x²]）或直接调大。
 *
 * λ 决定有效记忆 `1/(1-λ)` 个采样，按信号是否平稳来选：同一段 0.3 s 语音、order 20 实测
 * λ=0.999（≈1000 点）尾部残差 −49.8 dB；λ 拉到 0.9999 / 0.99999（≈1e4 / 1e5 点）跨过非平稳段后
 * 变成整段的折中，只剩 −21.2 / −19.0 dB。
 *
 * @note 整个类受 `QWQDSP_HAVE_EIGEN` 保护；ORDER 是编译期常量（模板参数），矩阵全部定长，无堆分配。
 * @note 类名里的 `FIlter` 是历史笔误（大写 I）；改名会破坏库外调用者，故保留。
 */
template <int ORDER>
class RLSFIlter {
public:
    using vec = Eigen::Matrix<double, ORDER, 1>;     // ORDER×1 列向量
    using mat = Eigen::Matrix<double, ORDER, ORDER>; // ORDER×ORDER 方阵

    /**
     * @brief 设定初始逆相关矩阵的对角值并复位
     * @param identity P0 = identity·I。想要"几乎没有先验、直接奔 LS 解"就给 1e6 量级；
     *                 默认 0.01 在常见信号量级下太强、收敛极慢（见类注释里的实测）
     */
    void Init(double identity = 0.01) noexcept {
        identity_val_ = identity;
        Reset();
    }

    /**
     * @brief 复位状态：`P = identity·I`、`w = 0`、输入/输出延迟线清零
     * @note 不改动 P0 与遗忘因子
     */
    void Reset() noexcept {
        identity_.setIdentity();
        p_ = identity_ * identity_val_;
        w_.setZero();
        latch_.setZero();
        iir_latch_.setZero();
    }

    /**
     * @brief 送一个采样：更新系数，并返回本次预测值
     * @param source 本次回归向量的第 0 项（当前输入样本）
     * @param target 期望输出
     * @return 预测值 `wᵀφ`，用的是**更新之前**的 w（先验预测）
     * @note 会改写内部 `w` / `P`；`source` 同时被移进延迟线，所以同值喂 source/target 会退化（见类注释）
     */
    double Tick(double source, double target) noexcept {
        for (int i = ORDER - 1; i > 0; --i) {
            latch_[i] = latch_[i - 1];
        }
        latch_[0] = source;

        // 先算预测与误差，再更新 P / w（a priori 形式）
        double const pred = w_.transpose() * latch_;
        err_ = target - pred;

        double const lamda_inv = 1.0 / forget_;
        p2_.noalias() = p_ * lamda_inv;
        k_.noalias() = (p_ * latch_) / (forget_ + (latch_.transpose() * p_ * latch_).value());
        p_.noalias() = (identity_ - k_ * latch_.transpose()) * p2_;
        w_.noalias() += k_ * err_;

        return pred;
    }

    /**
     * @brief 用当前系数做一次全极点合成
     *
     * 回归量取**输出**自己的延迟线（`iir_latch_`，与 `Tick` 的输入延迟线刻意分开），所以结构就是
     * `yy = wᵀ·y[n-1..n-ORDER] + gain·x`，即 `H(z) = 1 / (1 − Σ_j w[j-1]·z^-j)`——
     * 与用 `GetCoeff()` 的系数算出来的包络是同一个滤波器。
     *
     * @param x 激励
     * @return 本次合成输出
     * @note `gain = sqrt(err² + 1e-18)` 取的是**最近一次 `Tick` 的预测残差**当激励增益，
     *       所以要先跑 `Tick`（让 `err_` 有值）再用 `Filter`
     */
    double Filter(double x) noexcept {
        // 激励增益：取近期预测残差的量级
        double gain = std::sqrt(err_ * err_ + 1e-18f);

        double yy = w_.transpose() * iir_latch_ + gain * x;
        for (int i = ORDER - 1; i > 0; --i) {
            iir_latch_[i] = iir_latch_[i - 1];
        }
        iir_latch_[0] = yy;
        return yy;
    }

    /**
     * @brief 设置遗忘因子
     * @param forget 指数遗忘因子 λ，取 (0, 1]；有效记忆 ≈ `1/(1-λ)` 个采样
     */
    void SetForgetParam(float forget) noexcept {
        forget_ = forget;
    }

    /**
     * @brief 取当前预测系数
     * @return w，长度 ORDER，与 `Tick` 内部抽头一一对应：`pred = w[0]*source + w[1]*x[n-1] + ...`
     * @note 返回的是内部存储的引用，下一次 `Tick`/`Reset` 会改写它
     */
    vec const& GetCoeff() const noexcept {
        return w_;
    }
private:
    double forget_ = 0.999; // 遗忘因子 λ
    double err_{};          // 最近一次预测误差 target − pred；Filter() 拿它当激励增益
    double identity_val_{}; // P0 的对角值（Init 参数）
    mat p_;                 // 逆相关矩阵 P
    mat p2_;                // 每步临时量 P/λ
    mat identity_;          // 单位阵（预先建好，避免每步构造）
    vec w_;                 // 预测系数
    vec latch_;             // 输入延迟线（Tick 用），latch_[0] = 本次 source
    vec k_;                 // RLS 增益向量

    vec iir_latch_;         // 输出延迟线（Filter 用），与 latch_ 刻意分开
};
} // namespace qwqdsp_adaptive

#endif
