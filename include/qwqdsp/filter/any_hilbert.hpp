#pragma once

#include <array>
#include <complex>
#include <concepts>
#include <cstddef>
#include <tuple>
#include <utility>

namespace qwqdsp_filter::hilbert_coeff {

/// @brief 希尔伯特滤波器的一个二阶节系数
/// @note 由「半带低通 → 共轭极点对合并」得到：a = 2Re[r]、b = 2Re[r·p̄]、c = 2Re[p]、d = |p|²，
///       其中 (p, r) 是半带低通在 z⁻¹ 域并行分解的一对共轭极点与留数
///       （合并恒等式：(a - b·z⁻¹) / (1 - c·z⁻¹ + d·z⁻²)）；
///       ±90° 频移由实现在实部/虚部之间的耦合给出，不在这些系数里
struct Section {
    double a;
    double b;
    double c;
    double d;
};

/// @brief 系数集接口：kState 个二阶节 + 直接项
template <class T>
concept CAnyHilbertCoeff =
    requires {
        { T::kState } -> std::convertible_to<std::size_t>;
        { T::kDirect } -> std::convertible_to<double>;
        { T::kSections[0].a } -> std::convertible_to<double>;
        { T::kSections[0].b } -> std::convertible_to<double>;
        { T::kSections[0].c } -> std::convertible_to<double>;
        { T::kSections[0].d } -> std::convertible_to<double>;
    } && (std::tuple_size_v<decltype(T::kSections)> == static_cast<std::size_t>(T::kState));

/// @brief 32 阶 ellip 半带设计（通带波纹 0.1 dB / 阻带 150 dB）导出的希尔伯特滤波器系数
/// @note 由 tools/any_hilbert_gen.py 生成；合并二阶节与一阶复数并行形式的响应偏差 1.55e-13
struct Ellip32 {
    static constexpr std::size_t kState = 16;
    static constexpr double kDirect = 0.0016725282372669745;
    static constexpr std::array<Section, kState> kSections{{
        {1.5117691950e-01, -1.6150960039e-01, 1.3052689240e+00, 4.5547132049e-01},
        {-4.0029448985e-01, 3.1285144583e-01, 1.1135166236e+00, 5.3524281193e-01},
        {4.9774027954e-01, -2.3897509472e-01, 8.3791368666e-01, 6.4989752517e-01},
        {-4.4633141833e-01, 9.8686295334e-02, 5.8031118445e-01, 7.5706439789e-01},
        {3.2473980264e-01, 1.0576186212e-02, 3.8228975485e-01, 8.3944569600e-01},
        {-1.9996267924e-01, -6.5895779752e-02, 2.4442612941e-01, 8.9680198411e-01},
        {1.0336160540e-01, 7.8355485382e-02, 1.5330819938e-01, 9.3471366069e-01},
        {-4.0935386859e-02, -6.6911037420e-02, 9.4761364871e-02, 9.5907883314e-01},
        {6.8327132448e-03, 4.7011218501e-02, 5.7737148988e-02, 9.7449591774e-01},
        {7.9749493092e-03, -2.7827970196e-02, 3.4545671410e-02, 9.8416742751e-01},
        {-1.1561746436e-02, 1.3290298175e-02, 2.0113644317e-02, 9.9020946329e-01},
        {9.7259201265e-03, -4.0972773302e-03, 1.1189489279e-02, 9.9398370257e-01},
        {-6.0695498296e-03, -5.6086795657e-04, 5.7275727638e-03, 9.9635568114e-01},
        {2.6117239910e-03, 1.9432727706e-03, 2.4630086794e-03, 9.9787510426e-01},
        {-4.4278853955e-04, -1.3746116533e-03, 6.3520964134e-04, 9.9889639359e-01},
        {-1.1149996062e-04, 3.8518024932e-04, -1.8305928139e-04, 9.9965877371e-01},
    }};
};

} // namespace qwqdsp_filter::hilbert_coeff

namespace qwqdsp_filter {
/// @brief 由并行复数二阶节构成的希尔伯特滤波器
/// @tparam T 采样类型
/// @tparam TCoeff 系数集，kState / kSections / kDirect 全部从它取得（编译期常量，不占实例内存）
template <typename T, hilbert_coeff::CAnyHilbertCoeff TCoeff>
class AnyHilbert {
public:
    static constexpr std::size_t kState = TCoeff::kState;

    void Reset() noexcept {
        for (auto& s : states_) {
            s = State{};
        }
    }

    /// @brief 实输入一个样点，返回实部/虚部两路输出
    std::complex<T> Tick(T x) noexcept {
        std::complex<T> y = x * static_cast<T>(TCoeff::kDirect);
        [&]<std::size_t... I>(std::index_sequence<I...>) {
            (TickReal<I>(x, y), ...);
        }(std::make_index_sequence<kState>{});
        return y;
    }

    /// @brief 复数输入一个样点，返回复数输出
    std::complex<T> Tick(std::complex<T> x) noexcept {
        std::complex<T> y = x * static_cast<T>(TCoeff::kDirect);
        [&]<std::size_t... I>(std::index_sequence<I...>) {
            (TickComplex<I>(x, y), ...);
        }(std::make_index_sequence<kState>{});
        return y;
    }
private:
    /// @brief 单个二阶节的状态；系数来自模板，不占实例内存
    struct State {
        T sre1{};
        T sre2{};
        T sim1{};
        T sim2{};
    };

    /// @brief 实输入的单个二阶节
    template <std::size_t I>
    void TickReal(T x, std::complex<T>& y) noexcept {
        constexpr hilbert_coeff::Section kSection = TCoeff::kSections[I];
        T const a = static_cast<T>(kSection.a);
        T const b = static_cast<T>(kSection.b);
        T const c = static_cast<T>(kSection.c);
        T const d = static_cast<T>(kSection.d);

        State& s = states_[I];
        T const yre = x * a + s.sre1;
        T const yim = s.sim1;
        s.sre1 = -c * yim + s.sre2;
        s.sim1 = -b * x + yre * c + s.sim2;
        s.sre2 = d * yre;
        s.sim2 = d * yim;
        y += std::complex<T>{yre, yim};
    }

    /// @brief 复数输入的单个二阶节
    template <std::size_t I>
    void TickComplex(std::complex<T> x, std::complex<T>& y) noexcept {
        constexpr hilbert_coeff::Section kSection = TCoeff::kSections[I];
        T const a = static_cast<T>(kSection.a);
        T const b = static_cast<T>(kSection.b);
        T const c = static_cast<T>(kSection.c);
        T const d = static_cast<T>(kSection.d);

        T const xre = x.real();
        T const xim = x.imag();
        State& s = states_[I];
        T const yre = xre * a + s.sre1;
        T const yim = xim * a + s.sim1;
        s.sre1 = b * xim - c * yim + s.sre2;
        s.sim1 = -b * xre + yre * c + s.sim2;
        s.sre2 = d * yre;
        s.sim2 = d * yim;
        y += std::complex<T>{yre, yim};
    }

    std::array<State, kState> states_{};
};
} // namespace qwqdsp_filter
