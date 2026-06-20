// https://vicanek.de/articles/ReverseIIR.pdf
// 此文件实现了 IIR滤波器 的 线性相位补偿滤波器

#include <array>
#include <complex>
#include <qwqdsp/filter/rbj.hpp>
#include <qwqdsp/filter/biquad.hpp>
#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/convert.hpp>

template <int NDelay>
class Delay {
public:
    static constexpr int kMask = NDelay - 1;

    void Reset() noexcept {
        wpos_ = 0;
        buffer_.fill(0.0f);
    }

    float Get() noexcept {
        return buffer_[wpos_];
    }

    void Push(float x) noexcept {
        buffer_[wpos_] = x;
        ++wpos_;
        wpos_ &= kMask;
    }
private:
    std::array<float, NDelay> buffer_;
    int wpos_{};
};

template <int NState> requires (NState >= 0)
class ReverseState {
public:
    void Reset() noexcept {
        s1_.Reset();
        s2_.Reset();
    }

    void Set(float a1, float a2) noexcept {
        float re_c = -a1 / 2.0f;
        float c2 = a2;
        for (int i = 0; i < NState; ++i) {
            re_c = 2 * re_c * re_c - c2;
            c2 = c2 * c2;
        }
        a1_ = 2 * re_c;
        a2_ = c2;
    }

    float Tick(float x) noexcept {
        float y = a2_ * x + s2_.Get() + a1_ * s1_.Get();
        s2_.Push(s1_.Get());
        s1_.Push(x);
        return y;
    }
private:
    static constexpr int kDelay = 1 << NState;
    Delay<kDelay> s1_;
    Delay<kDelay> s2_;
    float a2_{};
    float a1_{};
};

class ReverseFilter {
public:
    void Reset() noexcept {
        biquad_.Reset();
        s0_.Reset();
        s1_.Reset();
        s2_.Reset();
        s3_.Reset();
        s4_.Reset();
        s5_.Reset();
    }

    void Set(qwqdsp_filter::BiquadCoeff c) noexcept {
        biquad_.Set(c);
        s0_.Set(c.a1, c.a2);
        s1_.Set(c.a1, c.a2);
        s2_.Set(c.a1, c.a2);
        s3_.Set(c.a1, c.a2);
        s4_.Set(c.a1, c.a2);
        s5_.Set(c.a1, c.a2);

        c.a1 = 0;
        c.a2 = 0;
        biquad2_.Set(c);
    }

    float Tick(float x) noexcept {
        x = s0_.Tick(x);
        x = s1_.Tick(x);
        x = s2_.Tick(x);
        x = s3_.Tick(x);
        x = s4_.Tick(x);
        x = s5_.Tick(x);
        x = biquad_.Tick(x);
        x = biquad2_.Tick(x);
        return x;
    }
private:
    qwqdsp_filter::Biquad biquad_;
    qwqdsp_filter::Biquad biquad2_;
    ReverseState<0> s0_;
    ReverseState<1> s1_;
    ReverseState<2> s2_;
    ReverseState<3> s3_;
    ReverseState<4> s4_;
    ReverseState<5> s5_;
};

int main() {
    ReverseFilter f;
    f.Reset();
    
    qwqdsp_filter::RBJ rbj;
    rbj.Lowpass(std::numbers::pi_v<float> / 2, std::numbers::sqrt2_v<float> / 2);
    auto coeff = rbj.ToBiquadCoeff();
    f.Set(coeff);

    float ir[1024]{1.0f};
    for (int i = 0; i < 1024; ++i) {
        ir[i] = f.Tick(ir[i]);
    }

    qwqdsp_spectral::RealFFT fft;
    fft.Init(1024);

    float gains[513];
    float phases[513];
    fft.FFTGainPhase(ir, gains, phases);
    
    for (int i = 0; i < 513; ++i) {
        gains[i] = qwqdsp::convert::Gain2Db<-100.0f>(gains[i]);
    }
}
