#pragma once
#include "qwqdsp/segement/slice.hpp"
#include <algorithm>
#include <cstddef>
#include <span>
#include <vector>

namespace qwqdsp_filter {
/**
 * @tparam kBatchSize 最大一次性处理多少采样，每隔这个数量才会执行一次移动采样，越大平均移动越少，同时内存增加
 */
template <size_t kBatchSize>
class FIRDirect {
public:
    void Reset() noexcept {
        std::fill(latch_.begin(), latch_.end(), 0.0f);
    }

    /**
     * @brief 滤波器的参数是h(0)...h(n-1)排列，且coeff的大小需要手动分配
     * @tparam Func void(std::vector<float>& coeff)
     */
    template <class Func>
    void SetCoeff(Func&& func) {
        func(coeff_);
        std::reverse(coeff_.begin(), coeff_.end());
        size_t const size = coeff_.size() + kBatchSize - 1;
        size_t const old_size = latch_.size();
        if (old_size < size) {
            latch_.resize(size);
            for (size_t i = 0; i < old_size; ++i) {
                size_t const target = size - i - 1;
                size_t const source = old_size - i - 1;
                latch_[target] = latch_[source];
            }
            std::fill_n(latch_.begin(), size - old_size, 0.0f);
        }
    }

    void Process(std::span<float> x) noexcept {
        qwqdsp_segement::Slice1D slice{x};
        size_t wpos = 0;
        while (!slice.IsEnd()) {
            // 首先向左移动缓冲区m个采样，复制m个采样到末尾，然后移动系数相加
            auto block = slice.GetSome(kBatchSize);
            for (size_t i = 0; i < latch_.size() - block.size(); ++i) {
                latch_[i] = latch_[i + block.size()];
            }
            std::copy(block.begin(), block.end(), latch_.end() - block.size());

            for (size_t i = 0; i < block.size(); ++i) {
                float sum = 0.0f;
                for (size_t j = 0; j < coeff_.size(); ++j) {
                    sum += coeff_[j] * latch_[i + j];
                }
                x[wpos++] = sum;
            }
        }
    }

    void Process(std::span<float> in, std::span<float> out) noexcept {
        qwqdsp_segement::Slice1D slice{in};
        size_t wpos = 0;
        while (!slice.IsEnd()) {
            // 首先向左移动缓冲区m个采样，复制m个采样到末尾，然后移动系数相加
            auto block = slice.GetSome(kBatchSize);
            for (size_t i = 0; i < latch_.size() - block.size(); ++i) {
                latch_[i] = latch_[i + block.size()];
            }
            std::copy(block.begin(), block.end(), latch_.end() - block.size());

            for (size_t i = 0; i < block.size(); ++i) {
                float sum = 0.0f;
                for (size_t j = 0; j < coeff_.size(); ++j) {
                    sum += coeff_[j] * latch_[i + j];
                }
                out[wpos++] = sum;
            }
        }
    }
private:
    std::vector<float> coeff_;
    std::vector<float> latch_;
};

class FIRTranspose {
public:
    void Reset() noexcept {
        std::fill(latch_.begin(), latch_.end(), 0.0f);
    }

    template <class Func>
        requires std::invocable<Func, std::vector<float>&>
    void SetCoeff(Func&& func) {
        func(coeff_);
        if (latch_.size() < coeff_.size() - 1) {
            latch_.resize(coeff_.size() - 1);
        }
    }

    void SetCoeff(std::span<const float> coeff) {
        coeff_.assign(coeff.begin(), coeff.end());
        latch_.resize(coeff.size() - 1);
    }

    float Tick(float x) noexcept {
        float const y = coeff_.front() * x + latch_.front();
        for (size_t i = 0; i < latch_.size() - 1; ++i) {
            latch_[i] = latch_[i + 1] + coeff_[i + 1] * x;
        }
        latch_.back() = coeff_.back() * x;
        return y;
    }

    void Process(std::span<float> x) noexcept {
        for (auto& s : x) {
            s = Tick(s);
        }
    }

    std::span<const float> GetCoeff() const noexcept {
        return coeff_;
    }
private:
    std::vector<float> coeff_;
    std::vector<float> latch_;
};

/**
 * @brief 基于循环缓冲区和系数双倍复制的优化 FIR 滤波器
 *
 * 该类通过将系数数组翻倍并反转排列，将延迟线中的循环索引访问
 * 转换为线性访问，避免了每采样计算索引取模的开销。适合高性能
 * 实时音频处理场景。
 *
 * @see FIRDirect, FIRTranspose
 */
class FirOptimise {
public:
    /**
     * @brief 重置滤波器状态
     *
     * 将延迟线清零，写位置归零
     */
    void Reset() noexcept {
        std::fill(lags_.begin(), lags_.end(), 0.0f);
        wpos_ = 0;
    }

    /**
     * @brief 设置滤波器系数
     * @param coeff 系数数组 h[0..N-1]，h[0] 对应最新输入采样
     *
     * 将系数翻倍并反转存储，使 Process 中的内积计算无需取模。
     * 同时初始化延迟线缓冲区大小。
     */
    void SetCoeff(std::span<const float> coeff) {
        coeff_.resize(coeff.size() * 2);
        auto it = std::reverse_copy(coeff.begin(), coeff.end(), coeff_.begin());
        std::reverse_copy(coeff.begin(), coeff.end(), it);

        lags_.resize(coeff.size());
        lag_size_ = static_cast<int>(lags_.size());
    }

    /**
     * @brief 处理单个采样点
     * @param x 输入采样
     * @return 滤波后的输出采样
     *
     * 先将新采样写入循环缓冲区，然后利用预排列的双倍系数数组
     * 计算内积，避免延迟线移动或取模运算。
     */
    float Process(float x) noexcept {
        lags_[wpos_] = x;
        ++wpos_;

        float y = 0.0f;
        int begin = lag_size_ - wpos_;
        for (int i = 0; i < lag_size_; ++i) {
            y += lags_[i] * coeff_[begin + i];
        }

        if (wpos_ >= lag_size_) {
            wpos_ = 0;
        }

        return y;
    }
private:
    std::vector<float> coeff_;
    std::vector<float> lags_;
    int wpos_{};
    int lag_size_{};
};
} // namespace qwqdsp_filter
