#pragma once
#include <span>
#include <vector>

namespace qwqdsp::adaptive {
template<class Sample, bool kLattice>
class BurgLP {
public:
    void Init(size_t block_len, size_t coeff_len) {
        b1_.resize(block_len);
        b2_.resize(block_len);
        aa_.resize(coeff_len);
    }

    float Process(std::span<float> x, std::span<float> a) {
        size_t const N = x.size();
        size_t const M = a.size();

        // (3)

        double gain = 0.0;
        for(size_t j = 0; j < N; j++)
        {
            gain += x[j] * x[j];
        }

        gain /= N;
        if(gain <= 0)
        {
            return 0.0;    // warning empty
        }

        // (9)

        b1_[0] = x[0];
        b2_[N - 2] = x[N - 1];
        for(size_t j = 1; j < N - 1; j++)
        {
            b1_[j] = b2_[j - 1] = x[j];
        }

        for(size_t i = 0; ; ++i)
        {
            // (7)

            double num = 0.0, denum = 0.0;
            for(size_t j = 0; j < N - i - 1; j++)
            {
                num += b1_[j] * b2_[j];
                denum += b1_[j] * b1_[j] + b2_[j] * b2_[j];
            }

            if(denum <= 0)
            {
                return 0.0;    // warning ill-conditioned
            }

            a[i] = 2.0 * num / denum;

            // (10)

            gain *= 1.0 - a[i] * a[i];

            // (5)

            if constexpr (!kLattice) {
                for(size_t j = 0; j+1 <= i; j++)
                {
                    a[j] = aa_[j] - a[i] * aa_[i - j - 1];
                }
            }

            if(i == M-1) break;

            // (8)  Watch out: i -> i+1

            for(size_t j = 0; j <= i; j++)
            {
                aa_[j] = a[j];
            }

            for(size_t j = 0; j <= N - i - 2; j++)
            {
                b1_[j] -= aa_[i] * b2_[j];
                b2_[j] = b2_[j + 1] - aa_[i] * b1_[j + 1];
            }
        }

        return gain;
    }
private:
    std::vector<float> b1_;
    std::vector<float> b2_;
    std::vector<float> aa_;
};
}