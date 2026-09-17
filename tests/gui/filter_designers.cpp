// ------------------------------------------------------------
// 滤波器设计器幅度对比
//
// 同一组 GUI 参数(类型 / fc / Q / gain / sigma)同时喂给三个设计器:
//   - qwqdsp_filter::RBJ         双线性变换的 Audio EQ Cookbook 公式
//   - qwqdsp_filter::MatchBiquad Vicanek 的匹配拟合(脉冲响应不变式)
//   - qwqdsp_filter::Ivantsov    Ivantsov 的"理想双线性"(decramped)拟合
// 并把 qwqdsp_filter::AnalogResponce 的模拟原型幅度画成参考线。
//
// 参考线的频率轴约定是**直接轴**: 模拟频率 = 数字频率(ωa = 2πf), 不做预畸变。
// 各设计器相对参考的表现:
//   - MatchBiquad 与参考重合(shelf 实测相差 < 0.0001dB), 设计目标就是直接轴;
//   - Ivantsov 也瞄准直接轴, 同样贴得很近, 但它是**拟合**而不是精确映射:
//     等效模拟频率会随频率上飘(低通 fc=1k σ=2 时 20kHz 处约 1.105·f),
//     所以高频段仍会离开参考(低通 20kHz 处约 1.7dB), 且偏差随 sigma 变化;
//   - RBJ 的设计目标是**预畸变轴** ωa = 2·fs·tan(πf/fs), 所以必然在高频段
//     大幅离开参考(低通 20kHz 处约 18dB)。
// 这些都是频率轴的约定差异, 不是设计错误, 界面上固定文字说明了这一点。
//
// 覆盖范围与三个设计器的非对称性:
//   - 两极点共有族: lowpass / highpass / bandpass(norm) / notch /
//     peaking / lowshelf / highshelf / allpass;
//   - RBJ 没有: 两极点 tiltshelf; Ivantsov 没有: bandpass(峰值=Q)、tiltshelf;
//   - 一阶: Ivantsov 提供 lp / hp / ap / highshelf / lowshelf;
//     MatchBiquad 提供 highshelf / lowshelf / tiltshelf; RBJ 不提供一阶;
//   - 未收录: RBJ::BandpassKeep0Precise(需要两个频率参数)、RBJ::Dicimate(固定 Q)。
//
// ⚠ raylib 内置字体只有 ASCII 字形: 界面文字必须用 ASCII, 中文只出现在注释里。
// ------------------------------------------------------------
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <format>
#include <functional>
#include <numbers>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "raylib.h"
#include "slider.hpp"

#include "qwqdsp/filter/analog_responce.hpp"
#include "qwqdsp/filter/biquad_coeff.hpp"
#include "qwqdsp/filter/ivantsov.hpp"
#include "qwqdsp/filter/match_biquad.hpp"
#include "qwqdsp/filter/rbj.hpp"

// ------------------------------------------------------------
// 常量与布局
// ------------------------------------------------------------
static constexpr int kWidth = 960;
static constexpr int kHeight = 620;

/// 采样率, 决定数字角频率与频率轴上界
static constexpr double kFs = 48000.0;
static constexpr double kFreqMin = 20.0;
static constexpr double kFreqMax = 20000.0;

/// 幅度面板纵轴范围 (dB), 固定不随类型变化, 便于横向比较不同设计器
static constexpr double kMagTopDb = 24.0;
static constexpr double kMagBottomDb = -96.0;

/// 幅度低于这个值时不再参与"与参考的最大偏差"统计: 陷波零点与阻带深处本身
/// 就在 -100dB 量级, 那里的 dB 差值是浮点噪声放大出来的, 没有意义
static constexpr float kDevFloorDb = -60.0f;

/// 兜底下限, 避免 log10(0) 给出 -inf
static constexpr float kVeryLowDb = -300.0f;

static constexpr float kPlotLeft = 70.0f;
static constexpr float kPlotWidth = 870.0f;
static constexpr float kPlotTop = 248.0f;
static constexpr float kPlotHeight = 310.0f;
static constexpr float kPlotBottom = kPlotTop + kPlotHeight;

/// 绘图区像素列数: 一列一条竖线段
static constexpr int kPlotColumns = static_cast<int>(kPlotWidth);

/// 每个像素列内的子采样数, 防止窄陷波被漏画
static constexpr int kSubSamples = 8;

static constexpr int kTitleFontSize = 14;
static constexpr int kLabelFontSize = 12;
static constexpr int kButtonFontSize = 12;
static constexpr int kButtonRows = 4;
static constexpr int kButtonColumns = 4;

/// 网格线颜色; 0dB 线单独给亮一点的颜色
static constexpr Color kGridColor{44, 44, 44, 255};
static constexpr Color kZeroLineColor{92, 92, 92, 255};
static constexpr Color kMarkerColor{122, 122, 122, 255};
static constexpr Color kRefColor{255, 255, 255, 255};
static constexpr Color kRbjColor{0, 228, 48, 255};
static constexpr Color kMbColor{255, 161, 0, 255};
static constexpr Color kIvColor{80, 180, 255, 255};

/**
 * @brief 曲线的虚线样式
 * @note 两个设计器在通带等处会完全重合, 若都画实线则后画的会把先画的盖掉,
 *       看不出"两者一致"。让两者各占虚线周期的一半且**按列互补**: 平坦段上就
 *       表现为两色短划交替出现, 重合与否一眼可辨; 两者合起来仍覆盖每一列,
 *       所以陡峭段上的曲线也不会断开。
 * @note 无论当前类型只有一个还是两个设计器都用虚线: 实线会把参考白带完全盖住,
 *       反而看不出它是否贴住参考; 留出的间隔正是"两者重合"的直观证据
 */
struct Dash {
    int on;     ///< 短划占几列
    int period; ///< 周期占几列
    int phase;  ///< 相位偏移(列), 让两个设计器错开
};
/// 实线样式
static constexpr Dash kSolidDash{0, 0, 0};
/// 三个设计器的短划周期相同、相位依次错开, 保证重合处能交替显示
static constexpr int kDashPeriod = 18;
static constexpr Dash kFirstDash{6, kDashPeriod, 0};
static constexpr Dash kSecondDash{6, kDashPeriod, 6};
static constexpr Dash kThirdDash{6, kDashPeriod, 12};
/// 参考线画成上下各外扩这么多像素的"带", 让设计器曲线落在带内也仍然看得见。
/// ⚠ DrawLineEx 的 thickness 是沿垂线方向外扩的: 对竖直线段是**横向**变粗,
///   所以这里必须用手工外扩, 不能靠 thickness。
static constexpr float kRefBandPx = 1.5f;

// ------------------------------------------------------------
// 滤波器类型表
// ------------------------------------------------------------

/// 面板覆盖的类型
/// @note 顺序必须与 kKinds 表严格一一对应: 按钮下标直接当作 Kind 用
enum class Kind : int {
    Lowpass = 0,
    Highpass,
    Bandpass,
    BandpassNorm,
    Notch,
    Peaking,
    Lowshelf,
    Highshelf,
    Tiltshelf,
    Allpass,
    OnepoleLowpass,
    OnepoleHighpass,
    OnepoleAllpass,
    OnepoleHighshelf,
    OnepoleLowshelf,
    OnepoleTiltshelf,
};

/// 每个类型的元信息: 名字、极点数、各设计器是否提供、哪些旋钮参与设计
struct KindInfo {
    char const* name;      ///< 按钮上的名字
    std::string_view note; ///< 额外提示(可为空)
    int poles;             ///< 极点个数, 1 或 2
    bool has_rbj;          ///< RBJ 是否提供该类型
    bool has_ivantsov;     ///< Ivantsov 是否提供该类型
    bool uses_q;           ///< Q 是否参与设计
    bool uses_gain;        ///< gain 是否参与设计
};

static constexpr std::array<KindInfo, 16> kKinds{
    {
     {"lowpass", "", 2, true, true, true, false},
     {"highpass", "", 2, true, true, true, false},
     {"bandpass Q", "peak gain = Q", 2, false, true, false},
     {"bandpass norm", "peak gain = 1 (Ivantsov/RBJ keep0)", 2, true, true, true, false},
     {"notch", "", 2, true, true, true, false},
     {"peaking", "peak gain = gain", 2, true, true, true, true},
     {"lowshelf", "DC gain = gain", 2, true, true, true, true},
     {"highshelf", "Nyquist gain = gain", 2, true, true, true, true},
     {"tiltshelf", "DC -g/2, Nyquist +g/2", 2, false, false, true, true},
     {"allpass", "magnitude is 0 dB everywhere; only phase differs", 2, true, true, true, false},
     {"onepole lowpass", "", 1, false, true, false, false},
     {"onepole highpass", "", 1, false, true, false, false},
     {"onepole allpass", "magnitude is 0 dB everywhere; only phase differs", 1, false, true, false, false},
     {"onepole highshelf", "DC 0 dB, Nyquist gain = gain", 1, false, true, false, true},
     {"onepole lowshelf", "DC gain = gain, Nyquist 0 dB", 1, false, true, false, true},
     {"onepole tiltshelf", "DC -g/2, Nyquist +g/2", 1, false, false, false, true},
     }
};

static constexpr std::array<char const*, kKinds.size()> kKindNames = [] {
    std::array<char const*, kKinds.size()> names{};
    for (size_t i = 0; i < kKinds.size(); ++i) {
        names[i] = kKinds[i].name;
    }
    return names;
}();

// ------------------------------------------------------------
// 设计与曲线缓存
// ------------------------------------------------------------

/// GUI 参数
struct Params {
    Kind kind = Kind::Lowpass;
    float fc = 1000.0f;   ///< 设计频率 (Hz)
    float q = 0.707f;     ///< 品质因子
    float gain_db = 6.0f; ///< peaking / shelf 的增益 (dB)
    float sigma = 2.0f;   ///< Ivantsov 的形状参数
    bool dirty = true;    ///< 参数变了需要重算
};

/// 一次设计的产物
struct DesignSet {
    /// 模拟原型在直接轴上的幅度 (dB), 输入数字角频率 (rad/sample)
    std::function<float(float)> analog_db;
    bool has_analog = false;
    qwqdsp_filter::BiquadCoeff rbj{};
    bool has_rbj = false;
    qwqdsp_filter::BiquadCoeff mb{};
    bool has_mb = false;
    qwqdsp_filter::BiquadCoeff iv{};
    bool has_iv = false;
};

/// 采样后的曲线: 每列一条竖直幅度范围(列内 min/max)
struct CurveCache {
    std::array<float, kPlotColumns> analog_low{};
    std::array<float, kPlotColumns> analog_high{};
    std::array<float, kPlotColumns> rbj_low{};
    std::array<float, kPlotColumns> rbj_high{};
    std::array<float, kPlotColumns> mb_low{};
    std::array<float, kPlotColumns> mb_high{};
    std::array<float, kPlotColumns> iv_low{};
    std::array<float, kPlotColumns> iv_high{};
    /// 与参考的最大偏差, 以及出现偏差的频率; dev_hz 保持 0 表示偏差恒为 0
    float rbj_dev_db = 0.0f;
    float rbj_dev_hz = 0.0f;
    float mb_dev_db = 0.0f;
    float mb_dev_hz = 0.0f;
    float iv_dev_db = 0.0f;
    float iv_dev_hz = 0.0f;
};

static Params g_params{};
static DesignSet g_design{};
static CurveCache g_curve{};

// ------------------------------------------------------------
// 坐标映射
// ------------------------------------------------------------

/**
 * @brief 把频率映射到绘图区横坐标
 * @param freq 频率 (Hz)
 * @return 绘图区内的 x 坐标
 */
static float freqToX(double freq) noexcept {
    double const t = std::log(freq / kFreqMin) / std::log(kFreqMax / kFreqMin);
    return kPlotLeft + static_cast<float>(t) * kPlotWidth;
}

/**
 * @brief 把绘图区横坐标映射回频率
 * @param x 绘图区内的 x 坐标
 * @return 对应的频率 (Hz)
 */
static double xToFreq(double x) noexcept {
    double const t = (x - static_cast<double>(kPlotLeft)) / static_cast<double>(kPlotWidth);
    return kFreqMin * std::pow(kFreqMax / kFreqMin, t);
}

/**
 * @brief 把幅度映射到绘图区纵坐标
 * @param db 幅度 (dB)
 * @return 绘图区内的 y 坐标
 */
static float dbToY(double db) noexcept {
    double const clamped = std::clamp(db, kMagBottomDb, kMagTopDb);
    double const t = (kMagTopDb - clamped) / (kMagTopDb - kMagBottomDb);
    return kPlotTop + static_cast<float>(t) * kPlotHeight;
}

// ------------------------------------------------------------
// 设计
// ------------------------------------------------------------

/**
 * @brief 复数频响转 dB, 并对零响应兜底
 * @param h 复频响
 * @return 幅度 (dB), 不低于 kVeryLowDb
 */
static float toDb(std::complex<float> h) noexcept {
    float const magnitude = std::abs(h);
    if (!(magnitude > 0.0f)) {
        return kVeryLowDb;
    }
    return std::max(20.0f * std::log10(magnitude), kVeryLowDb);
}

/**
 * @brief 双二阶节在指定数字频率上的幅度
 * @param coeff 双二阶系数
 * @param w 数字角频率 (rad/sample)
 * @return 幅度 (dB)
 */
static float digitalDb(qwqdsp_filter::BiquadCoeff const& coeff, float w) noexcept {
    std::complex<float> const z{std::cos(w), std::sin(w)};
    return toDb(coeff.DigitalResonpoce(z));
}

/**
 * @brief 按当前 GUI 参数做一次设计
 * @param params 设计参数
 * @return 模拟原型参考与两个设计器的系数; has_* 表示该项是否可用
 * @note 参考取**直接轴**: 模拟原型的截止频率就是 fc 对应的数字角频率 wc,
 *       不做预畸变。这正是 MatchBiquad 的设计目标(RBJ 的目标是预畸变轴)。
 */
static DesignSet makeDesign(Params const& params) {
    using qwqdsp_filter::AnalogResponce;
    using qwqdsp_filter::Ivantsov;
    using qwqdsp_filter::MatchBiquad;
    using qwqdsp_filter::RBJ;

    float const wc = static_cast<float>(2.0 * std::numbers::pi * static_cast<double>(params.fc) / kFs);
    float const q = params.q;
    float const g = params.gain_db;
    // 模拟原型的幅度参数约定: 两极点 shelf 传 10^(db/80), peaking 与单极点 shelf 传 10^(db/40)
    float const sqrt_a80 = std::pow(10.0f, g / 80.0f);
    float const a40 = std::pow(10.0f, g / 40.0f);

    RBJ rbj;
    MatchBiquad mb;
    Ivantsov iv;
    float const sigma = params.sigma;
    DesignSet out;

    switch (params.kind) {
        case Kind::Lowpass:
            out.analog_db = [wc, q](float wa) {
                AnalogResponce a;
                return toDb(a.Lowpass(wa, wc, q));
            };
            rbj.Lowpass(wc, q);
            out.rbj = rbj.ToBiquadCoeff();
            out.has_rbj = true;
            out.mb = mb.Lowpass(wc, q);
            out.has_mb = true;
            out.iv = iv.Lowpass(wc, q, sigma);
            out.has_iv = true;
            break;
        case Kind::Highpass:
            out.analog_db = [wc, q](float wa) {
                AnalogResponce a;
                return toDb(a.Highpass(wa, wc, q));
            };
            rbj.Highpass(wc, q);
            out.rbj = rbj.ToBiquadCoeff();
            out.has_rbj = true;
            out.mb = mb.Highpass(wc, q);
            out.has_mb = true;
            out.iv = iv.Highpass(wc, q, sigma);
            out.has_iv = true;
            break;
        case Kind::Bandpass:
            out.analog_db = [wc, q](float wa) {
                AnalogResponce a;
                return toDb(a.Bandpass(wa, wc, q));
            };
            rbj.Bandpass(wc, q);
            out.rbj = rbj.ToBiquadCoeff();
            out.has_rbj = true;
            out.mb = mb.Bandpass(wc, q);
            out.has_mb = true;
            break;
        case Kind::BandpassNorm:
            out.analog_db = [wc, q](float wa) {
                AnalogResponce a;
                return toDb(a.NormBandpass(wa, wc, q));
            };
            rbj.BandpassKeep0(wc, q);
            out.rbj = rbj.ToBiquadCoeff();
            out.has_rbj = true;
            out.mb = mb.NormBandpass(wc, q);
            out.has_mb = true;
            out.iv = iv.Bandpass(wc, q, sigma);
            out.has_iv = true;
            break;
        case Kind::Notch:
            out.analog_db = [wc, q](float wa) {
                AnalogResponce a;
                return toDb(a.Notch(wa, wc, q));
            };
            rbj.Notch(wc, q);
            out.rbj = rbj.ToBiquadCoeff();
            out.has_rbj = true;
            out.mb = mb.Notch(wc, q);
            out.has_mb = true;
            out.iv = iv.Notch(wc, q, sigma);
            out.has_iv = true;
            break;
        case Kind::Peaking:
            out.analog_db = [wc, q, a40](float wa) {
                AnalogResponce a;
                return toDb(a.Peaking(wa, wc, q, a40));
            };
            rbj.Peak(wc, q, g);
            out.rbj = rbj.ToBiquadCoeff();
            out.has_rbj = true;
            out.mb = mb.Peaking(wc, q, g);
            out.has_mb = true;
            out.iv = iv.Peaking(wc, q, g, sigma);
            out.has_iv = true;
            break;
        case Kind::Lowshelf:
            out.analog_db = [wc, q, sqrt_a80](float wa) {
                AnalogResponce a;
                return toDb(a.Lowshelf(wa, wc, q, sqrt_a80));
            };
            rbj.Lowshelf(wc, q, g);
            out.rbj = rbj.ToBiquadCoeff();
            out.has_rbj = true;
            out.mb = mb.Lowshelf(wc, q, g);
            out.has_mb = true;
            out.iv = iv.Lowshelf(wc, q, g, sigma);
            out.has_iv = true;
            break;
        case Kind::Highshelf:
            out.analog_db = [wc, q, sqrt_a80](float wa) {
                AnalogResponce a;
                return toDb(a.Highshelf(wa, wc, q, sqrt_a80));
            };
            rbj.HighShelf(wc, q, g);
            out.rbj = rbj.ToBiquadCoeff();
            out.has_rbj = true;
            out.mb = mb.Highshelf(wc, q, g);
            out.has_mb = true;
            out.iv = iv.Highshelf(wc, q, g, sigma);
            out.has_iv = true;
            break;
        case Kind::Tiltshelf:
            // RBJ 没有倾斜 shelf, 这一格只有 MatchBiquad 与模拟原型
            out.analog_db = [wc, q, sqrt_a80](float wa) {
                AnalogResponce a;
                return toDb(a.Tiltshelf(wa, wc, q, sqrt_a80));
            };
            out.mb = mb.Tiltshelf(wc, q, g);
            out.has_mb = true;
            break;
        case Kind::Allpass:
            out.analog_db = [wc, q](float wa) {
                AnalogResponce a;
                return toDb(a.Allpass(wa, wc, q));
            };
            rbj.Allpass(wc, q);
            out.rbj = rbj.ToBiquadCoeff();
            out.has_rbj = true;
            out.mb = mb.Allpass(wc, q);
            out.has_mb = true;
            out.iv = iv.Allpass(wc, q, sigma);
            out.has_iv = true;
            break;
        case Kind::OnepoleHighshelf:
            out.analog_db = [wc, a40](float wa) {
                AnalogResponce a;
                return toDb(a.HighshelfOnepole(wa, wc, a40));
            };
            out.mb = mb.HighshelfOnepole(wc, g);
            out.has_mb = true;
            out.iv = iv.HighshelfOnepole(wc, g, sigma);
            out.has_iv = true;
            break;
        case Kind::OnepoleLowshelf:
            out.analog_db = [wc, a40](float wa) {
                AnalogResponce a;
                return toDb(a.LowshelfOnepole(wa, wc, a40));
            };
            out.mb = mb.LowshelfOnepole(wc, g);
            out.has_mb = true;
            out.iv = iv.LowshelfOnepole(wc, g, sigma);
            out.has_iv = true;
            break;
        case Kind::OnepoleTiltshelf:
            out.analog_db = [wc, a40](float wa) {
                AnalogResponce a;
                return toDb(a.TiltshelfOnepole(wa, wc, a40));
            };
            out.mb = mb.TiltshelfOnepole(wc, g);
            out.has_mb = true;
            break;
        case Kind::OnepoleLowpass:
            out.analog_db = [wc](float wa) {
                AnalogResponce a;
                return toDb(a.LowpassOnepole(wa, wc));
            };
            out.iv = iv.LowpassOnepole(wc, sigma);
            out.has_iv = true;
            break;
        case Kind::OnepoleHighpass:
            out.analog_db = [wc](float wa) {
                AnalogResponce a;
                return toDb(a.HighpassOnepole(wa, wc));
            };
            out.iv = iv.HighpassOnepole(wc, sigma);
            out.has_iv = true;
            break;
        case Kind::OnepoleAllpass:
            out.analog_db = [wc](float wa) {
                AnalogResponce a;
                return toDb(a.AllpassOnepole(wa, wc));
            };
            out.iv = iv.AllpassOnepole(wc, sigma);
            out.has_iv = true;
            break;
    }
    out.has_analog = true;
    return out;
}

/**
 * @brief 逐像素列采样三条曲线
 * @param design 设计结果
 * @param curve 输出的曲线缓存
 * @note 每列取 kSubSamples 个子采样; 幅度按列取 min/max 画成竖线段,
 *       窄陷波与陡沿才不会被漏掉。偏差统计只算参考高于 kDevFloorDb 的频点。
 */
static void sampleCurves(DesignSet const& design, CurveCache& curve) {
    curve.rbj_dev_db = 0.0f;
    curve.rbj_dev_hz = 0.0f;
    curve.mb_dev_db = 0.0f;
    curve.mb_dev_hz = 0.0f;
    curve.iv_dev_db = 0.0f;
    curve.iv_dev_hz = 0.0f;

    for (int column = 0; column < kPlotColumns; ++column) {
        float analog_low = 1.0e30f;
        float analog_high = -1.0e30f;
        float rbj_low = 1.0e30f;
        float rbj_high = -1.0e30f;
        float mb_low = 1.0e30f;
        float mb_high = -1.0e30f;
        float iv_low = 1.0e30f;
        float iv_high = -1.0e30f;

        for (int sub = 0; sub < kSubSamples; ++sub) {
            double const x = static_cast<double>(kPlotLeft) + static_cast<double>(column)
                           + (static_cast<double>(sub) + 0.5) / static_cast<double>(kSubSamples);
            double const freq = xToFreq(x);
            float const w = static_cast<float>(2.0 * std::numbers::pi * freq / kFs);

            float const analog_db = design.analog_db(w);
            float const rbj_db = design.has_rbj ? digitalDb(design.rbj, w) : 0.0f;
            float const mb_db = design.has_mb ? digitalDb(design.mb, w) : 0.0f;
            float const iv_db = design.has_iv ? digitalDb(design.iv, w) : 0.0f;

            analog_low = std::min(analog_low, analog_db);
            analog_high = std::max(analog_high, analog_db);
            if (design.has_rbj) {
                rbj_low = std::min(rbj_low, rbj_db);
                rbj_high = std::max(rbj_high, rbj_db);
            }
            if (design.has_mb) {
                mb_low = std::min(mb_low, mb_db);
                mb_high = std::max(mb_high, mb_db);
            }
            if (design.has_iv) {
                iv_low = std::min(iv_low, iv_db);
                iv_high = std::max(iv_high, iv_db);
            }

            if (analog_db > kDevFloorDb) {
                float const freq_f = static_cast<float>(freq);
                if (design.has_rbj) {
                    float const dev = std::abs(rbj_db - analog_db);
                    if (dev > curve.rbj_dev_db) {
                        curve.rbj_dev_db = dev;
                        curve.rbj_dev_hz = freq_f;
                    }
                }
                if (design.has_mb) {
                    float const dev = std::abs(mb_db - analog_db);
                    if (dev > curve.mb_dev_db) {
                        curve.mb_dev_db = dev;
                        curve.mb_dev_hz = freq_f;
                    }
                }
                if (design.has_iv) {
                    float const dev = std::abs(iv_db - analog_db);
                    if (dev > curve.iv_dev_db) {
                        curve.iv_dev_db = dev;
                        curve.iv_dev_hz = freq_f;
                    }
                }
            }
        }

        curve.analog_low[column] = analog_low;
        curve.analog_high[column] = analog_high;
        curve.rbj_low[column] = rbj_low;
        curve.rbj_high[column] = rbj_high;
        curve.mb_low[column] = mb_low;
        curve.mb_high[column] = mb_high;
        curve.iv_low[column] = iv_low;
        curve.iv_high[column] = iv_high;
    }
}

// ------------------------------------------------------------
// 绘制辅助
// ------------------------------------------------------------

/**
 * @brief 画按网格排布的一排可点击选择按钮
 * @param area 整片按钮占据的矩形
 * @param names 每项名称
 * @param columns 每行几个
 * @param selected 当前选中项下标
 * @return 点击后应选中的项, 未点击则原样返回
 */
static size_t drawSelectorGrid(Rectangle area, std::span<char const* const> names, size_t columns, size_t selected) {
    auto const mouse = GetMousePosition();
    size_t result = selected;
    size_t const rows = (names.size() + columns - 1) / columns;
    float const cell_w = area.width / static_cast<float>(columns);
    float const cell_h = area.height / static_cast<float>(rows);

    for (size_t i = 0; i < names.size(); ++i) {
        Rectangle const item{
            area.x + static_cast<float>(i % columns) * cell_w,
            area.y + static_cast<float>(i / columns) * cell_h,
            cell_w - 3.0f,
            cell_h - 2.0f,
        };
        bool const active = (i == selected);
        bool const hover = CheckCollisionPointRec(mouse, item);
        if (hover && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            result = i;
        }
        if (active) {
            DrawRectangleRec(item, RAYWHITE);
        }
        else {
            DrawRectangleLinesEx(item, 1.0f, hover ? WHITE : GRAY);
        }
        DrawText(names[i], static_cast<int>(item.x) + 4, static_cast<int>(item.y) + 5, kButtonFontSize,
                 active ? BLACK : (hover ? WHITE : GRAY));
    }
    return result;
}

/**
 * @brief 画一条虚线竖线
 * @param x 横坐标
 * @param y0 起点纵坐标
 * @param y1 终点纵坐标
 * @param color 颜色
 */
static void drawDashedVertical(float x, float y0, float y1, Color color) {
    constexpr float kDash = 6.0f;
    for (float y = y0; y < y1; y += 2.0f * kDash) {
        DrawLineV({x, y}, {x, std::min(y + kDash, y1)}, color);
    }
}

/**
 * @brief 画一条虚线横线
 * @param y 纵坐标
 * @param x0 起点横坐标
 * @param x1 终点横坐标
 * @param color 颜色
 */
static void drawDashedHorizontal(float y, float x0, float x1, Color color) {
    constexpr float kDash = 6.0f;
    for (float x = x0; x < x1; x += 2.0f * kDash) {
        DrawLineV({x, y}, {std::min(x + kDash, x1), y}, color);
    }
}

/**
 * @brief 画频率轴与幅度轴的网格和刻度
 * @note 频率刻度画在面板下方, 幅度刻度画在面板左侧
 */
static void drawGrid() {
    static constexpr std::array<double, 10> kFreqTicks{20.0,   50.0,   100.0,  200.0,   500.0,
                                                       1000.0, 2000.0, 5000.0, 10000.0, 20000.0};
    static constexpr std::array<double, 6> kDbTicks{20.0, 0.0, -20.0, -40.0, -60.0, -80.0};

    DrawRectangleLinesEx({kPlotLeft, kPlotTop, kPlotWidth, kPlotHeight}, 1.0f, GRAY);

    for (double db : kDbTicks) {
        float const y = dbToY(db);
        bool const is_zero = (db == 0.0);
        DrawLineV({kPlotLeft, y}, {kPlotLeft + kPlotWidth, y}, is_zero ? kZeroLineColor : kGridColor);
        auto const text = std::format("{:+.0f}", db);
        int const text_width = MeasureText(text.c_str(), kLabelFontSize);
        DrawText(text.c_str(), static_cast<int>(kPlotLeft) - 8 - text_width, static_cast<int>(y) - 6, kLabelFontSize,
                 is_zero ? LIGHTGRAY : GRAY);
    }

    for (double freq : kFreqTicks) {
        float const x = freqToX(freq);
        DrawLineV({x, kPlotTop}, {x, kPlotTop + kPlotHeight}, kGridColor);
        auto const text = freq >= 1000.0 ? std::format("{:.0f}k", freq / 1000.0) : std::format("{:.0f}", freq);
        int const text_width = MeasureText(text.c_str(), kLabelFontSize);
        DrawText(text.c_str(), static_cast<int>(x) - text_width / 2, static_cast<int>(kPlotBottom) + 4, kLabelFontSize,
                 GRAY);
    }

    DrawText("Hz", static_cast<int>(kPlotLeft + kPlotWidth) + 4, static_cast<int>(kPlotBottom) + 4, kLabelFontSize,
             GRAY);
    DrawText("dB", static_cast<int>(kPlotLeft) - 32, static_cast<int>(kPlotTop) - 15, kLabelFontSize, GRAY);
}

/**
 * @brief 画设计点标记: fc 竖线与目标增益横线
 * @param params 当前参数
 * @param info 当前类型的元信息
 */
static void drawMarkers(Params const& params, KindInfo const& info) {
    float const x = freqToX(params.fc);
    drawDashedVertical(x, kPlotTop, kPlotBottom, kMarkerColor);
    DrawText("fc", static_cast<int>(x) + 3, static_cast<int>(kPlotTop) + 3, kLabelFontSize, LIGHTGRAY);

    if (info.uses_gain) {
        float const y = dbToY(params.gain_db);
        drawDashedHorizontal(y, kPlotLeft, kPlotLeft + kPlotWidth, kMarkerColor);
        auto const text = std::format("gain {:+.1f} dB", params.gain_db);
        int const text_width = MeasureText(text.c_str(), kLabelFontSize);
        DrawText(text.c_str(), static_cast<int>(kPlotLeft + kPlotWidth) - text_width - 4, static_cast<int>(y) - 15,
                 kLabelFontSize, LIGHTGRAY);
    }
}

/**
 * @brief 把一条曲线画成一串竖线段
 * @param low 每列的幅度下界 (dB)
 * @param high 每列的幅度上界 (dB)
 * @param color 颜色
 * @param band_px 每列在上下各外扩多少像素, 0 表示不扩
 * @param dash 虚线样式, 传 kSolidDash 画实线
 * @note 列内 min/max 撑开成竖线段: 相邻列在频率上重叠, 画出来是连续的,
 *       而同列内的剧烈变化(陡沿、陷波)也不会被漏掉
 */
static void drawCurve(std::array<float, kPlotColumns> const& low, std::array<float, kPlotColumns> const& high,
                      Color color, float band_px, Dash dash = kSolidDash) {
    for (int column = 0; column < kPlotColumns; ++column) {
        if (dash.period > 0 && (column + dash.phase) % dash.period >= dash.on) {
            continue;
        }
        float const x = kPlotLeft + static_cast<float>(column) + 0.5f;
        float const y0 = dbToY(low[column]);
        float const y1 = dbToY(high[column]);
        DrawLineEx({x, std::min(y0, y1) - 0.5f - band_px}, {x, std::max(y0, y1) + 0.5f + band_px}, 1.0f, color);
    }
}

/**
 * @brief 画图例
 * @param bound 图例占据的矩形
 * @param design 设计结果, 用来标注某设计器是否提供该类型
 * @note 图例画在绘图区**外面**(旋钮右侧的信息区): 放在面板里会盖住曲线,
 *       半透明底也会把穿过的曲线压暗, 反而看不清要对比的东西
 */
static void drawLegend(Rectangle bound, DesignSet const& design) {
    constexpr float kSwatch = 26.0f;
    constexpr float kLineHeight = 17.0f;
    constexpr int kRows = 4;

    struct Row {
        char const* name;
        Color color;
        bool available;
        Dash dash;
    };
    std::array<Row, kRows> const rows{
        {
         {"analog prototype", kRefColor, true, kSolidDash},
         {"RBJ", kRbjColor, design.has_rbj, kFirstDash},
         {"MatchBiquad", kMbColor, design.has_mb, kSecondDash},
         {"Ivantsov", kIvColor, design.has_iv, kThirdDash},
         }
    };

    for (size_t i = 0; i < rows.size(); ++i) {
        float const y = bound.y + static_cast<float>(i) * kLineHeight + 8.0f;
        bool const available = rows[i].available;
        Color const color = available ? rows[i].color : DARKGRAY;
        float const x0 = bound.x;
        float const x1 = x0 + kSwatch;
        if (rows[i].dash.period <= 0) {
            DrawLineEx({x0, y}, {x1, y}, 3.0f, color);
        }
        else {
            for (float x = x0; x < x1; x += 1.0f) {
                if (static_cast<int>(x - x0) % rows[i].dash.period < rows[i].dash.on) {
                    DrawLineEx({x, y}, {x + 1.0f, y}, 3.0f, color);
                }
            }
        }
        auto const text = available ? std::string{rows[i].name} : std::format("{} (n/a)", rows[i].name);
        DrawText(text.c_str(), static_cast<int>(x1 + 6.0f), static_cast<int>(y) - 6, kLabelFontSize,
                 available ? RAYWHITE : DARKGRAY);
    }
}

/**
 * @brief 拼出当前类型的一句话说明(哪些旋钮参与、RBJ 是否有该类型)
 * @param info 类型元信息
 * @return 说明文字
 * @note 返回的是 ASCII 文本: raylib 内置字体没有中文字形
 */
static std::string kindHint(KindInfo const& info) {
    std::vector<std::string> parts;
    if (!info.has_rbj) {
        parts.emplace_back("RBJ does not provide it");
    }
    if (!info.has_ivantsov) {
        parts.emplace_back("Ivantsov does not provide it");
    }
    if (!info.uses_q) {
        parts.emplace_back("Q unused");
    }
    if (!info.uses_gain) {
        parts.emplace_back("gain unused");
    }

    std::string hint;
    for (size_t i = 0; i < parts.size(); ++i) {
        if (i > 0) {
            hint += ", ";
        }
        hint += parts[i];
    }
    if (!info.note.empty()) {
        if (!hint.empty()) {
            hint += ". ";
        }
        hint += info.note;
    }
    if (hint.empty()) {
        hint = "both designers provide it";
    }
    return hint;
}

/**
 * @brief 拼出某个设计器与参考的偏差说明
 * @param available 该设计器是否提供当前类型
 * @param dev_db 最大偏差 (dB)
 * @param dev_hz 偏差出现的频率 (Hz); 保持 0 表示偏差恒为 0
 * @return 说明文字
 */
static std::string deviationText(bool available, float dev_db, float dev_hz) {
    if (!available) {
        return "n/a for this kind";
    }
    if (dev_hz <= 0.0f || dev_db < 0.0005f) {
        return "coincides with the ref (< 0.001 dB)";
    }
    return std::format("max dev {:.2f} dB @ {:.0f} Hz", dev_db, dev_hz);
}

/**
 * @brief 按像素宽度把文字折成若干行, 尽量在空格处断开
 * @param text 原始文字
 * @param font_size 字号
 * @param max_width 每行最大像素宽度
 * @param max_lines 最多几行, 超出部分直接丢弃
 * @return 折行后的各行(至少一行)
 */
static std::vector<std::string> wrapText(std::string const& text, int font_size, float max_width, size_t max_lines) {
    std::vector<std::string> lines;
    std::string current;
    size_t i = 0;
    while (i < text.size()) {
        size_t const space = text.find(' ', i);
        std::string const word = text.substr(i, space == std::string::npos ? std::string::npos : space - i);
        std::string candidate = current.empty() ? word : current + " " + word;
        if (MeasureText(candidate.c_str(), font_size) <= static_cast<int>(max_width) || current.empty()) {
            current = candidate;
        }
        else {
            lines.push_back(current);
            if (lines.size() >= max_lines) {
                return lines;
            }
            current = word;
        }
        if (space == std::string::npos) {
            break;
        }
        i = space + 1;
    }
    if (!current.empty() && lines.size() < max_lines) {
        lines.push_back(current);
    }
    if (lines.empty()) {
        lines.emplace_back();
    }
    return lines;
}

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
int main() {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWidth, kHeight, "filter designers: magnitude vs analog prototype");
    SetTargetFPS(60);

    // ----- 类型选择按钮 -----
    Rectangle const kind_area{
        10.0f,
        22.0f,
        940.0f,
        static_cast<float>(kButtonRows) * 22.0f,
    };

    // ----- 参数旋钮 -----
    float const knob_width = 104.0f;
    float const knob_height = 74.0f;
    float const knob_y = 116.0f;
    auto knobBound = [&](size_t index) {
        return Rectangle{
            10.0f + static_cast<float>(index) * knob_width,
            knob_y,
            knob_width - 4.0f,
            knob_height,
        };
    };
    auto setupKnob = [](Knob& knob, Rectangle bound, char const* title) {
        knob.set_bound(bound);
        knob.set_bg_color(BLACK);
        knob.set_fore_color(RAYWHITE);
        knob.set_title(title);
        knob.set_name_font_size(12);
        knob.set_number_font_size(12);
    };

    // 频率旋钮走 log2 域: 量程 4.32 个倍频程, 步长 0.005 -> 满量程约 864 像素
    float const fc_log_max = static_cast<float>(std::log2(kFreqMax / kFreqMin));
    Knob fc_knob;
    setupKnob(fc_knob, knobBound(0), "fc");
    fc_knob.set_range(0.0f, fc_log_max, 0.005f, static_cast<float>(std::log2(g_params.fc / kFreqMin)));
    fc_knob.value_to_text_function = [](float value) {
        return std::format("{:.1f} Hz", kFreqMin * std::exp2(static_cast<double>(value)));
    };
    fc_knob.on_value_change = [](float value) {
        g_params.fc = static_cast<float>(kFreqMin * std::exp2(static_cast<double>(value)));
        g_params.dirty = true;
    };

    Knob q_knob;
    setupKnob(q_knob, knobBound(1), "Q");
    q_knob.set_range(0.1f, 10.0f, 0.01f, g_params.q);
    q_knob.value_to_text_function = [](float value) { return std::format("{:.3f}", value); };
    q_knob.on_value_change = [](float value) {
        g_params.q = value;
        g_params.dirty = true;
    };

    Knob sigma_knob;
    setupKnob(sigma_knob, knobBound(3), "sigma");
    sigma_knob.set_range(1.0f, 4.0f, 0.005f, g_params.sigma);
    sigma_knob.value_to_text_function = [](float value) { return std::format("{:.3f}", value); };
    sigma_knob.on_value_change = [](float value) {
        g_params.sigma = value;
        g_params.dirty = true;
    };

    Knob gain_knob;
    setupKnob(gain_knob, knobBound(2), "gain");
    gain_knob.set_range(-18.0f, 18.0f, 0.1f, g_params.gain_db);
    gain_knob.value_to_text_function = [](float value) { return std::format("{:+.1f} dB", value); };
    gain_knob.on_value_change = [](float value) {
        g_params.gain_db = value;
        g_params.dirty = true;
    };

    while (!WindowShouldClose()) {
        BeginDrawing();
        {
            ClearBackground(BLACK);

            DrawText("filter designers: magnitude response vs analog prototype", 10, 4, kTitleFontSize, RAYWHITE);

            // ----- 选类型 -----
            size_t const kind_index =
                drawSelectorGrid(kind_area, std::span<char const* const>{kKindNames.data(), kKindNames.size()},
                                 static_cast<size_t>(kButtonColumns), static_cast<size_t>(g_params.kind));
            if (kind_index != static_cast<size_t>(g_params.kind)) {
                g_params.kind = static_cast<Kind>(kind_index);
                g_params.dirty = true;
            }
            KindInfo const& info = kKinds[static_cast<size_t>(g_params.kind)];

            // ----- 旋钮: 不参与当前设计的旋钮直接不画(不画就不接受输入) -----
            // 同时把**可见的**旋钮依次排到最左边, 不留空洞: 固定位置的话,
            // 被隐藏的旋钮会空出一块, 后面的 sigma 就会被推到信息文字上
            q_knob.SetEnable(info.uses_q);
            gain_knob.SetEnable(info.uses_gain);
            sigma_knob.SetEnable(info.has_ivantsov);
            {
                std::array<Knob*, 4> const knobs{&fc_knob, &q_knob, &gain_knob, &sigma_knob};
                std::array<bool, 4> const visible{true, info.uses_q, info.uses_gain, info.has_ivantsov};
                size_t slot = 0;
                for (size_t i = 0; i < knobs.size(); ++i) {
                    if (visible[i]) {
                        knobs[i]->set_bound(knobBound(slot));
                        ++slot;
                    }
                }
            }
            fc_knob.display();
            q_knob.display();
            gain_knob.display();
            sigma_knob.display();

            // ----- 参数变了才重算 -----
            if (g_params.dirty) {
                g_design = makeDesign(g_params);
                sampleCurves(g_design, g_curve);
                g_params.dirty = false;
            }

            // ----- 面板 -----
            drawGrid();
            drawMarkers(g_params, info);
            // 参考线画成一条"带"(实线, 上下各外扩 1.5px), 设计器的虚线正好画在带
            // **中心**: 全部重合时表现为带中央三色短划交替, 一眼能看出"都贴着参考";
            // 偏离时短划就跑到带外面去了。
            // 三个设计器各占虚线周期的 1/3 且相位错开, 这样任意两者重合时都能交替
            // 显示; 它们合起来仍覆盖每一列, 所以陡峭段上的曲线也不会断开。
            drawCurve(g_curve.analog_low, g_curve.analog_high, kRefColor, kRefBandPx);
            if (g_design.has_rbj) {
                drawCurve(g_curve.rbj_low, g_curve.rbj_high, kRbjColor, 0.0f, kFirstDash);
            }
            if (g_design.has_mb) {
                drawCurve(g_curve.mb_low, g_curve.mb_high, kMbColor, 0.0f, kSecondDash);
            }
            if (g_design.has_iv) {
                drawCurve(g_curve.iv_low, g_curve.iv_high, kIvColor, 0.0f, kThirdDash);
            }

            // ----- 文字信息(旋钮右侧) 与 图例(最右侧) -----
            {
                float const info_x = 440.0f;
                float const info_width = 296.0f;
                float const line_height = 14.0f;
                int y = static_cast<int>(knob_y) + 1;
                auto drawLine = [&](std::string const& text, Color color) {
                    DrawText(text.c_str(), static_cast<int>(info_x), y, kLabelFontSize, color);
                    y += static_cast<int>(line_height);
                };
                // 折行后逐行画, 避免把说明截断成半句话
                auto drawWrapped = [&](std::string const& text, Color color, size_t max_lines) {
                    for (auto const& line : wrapText(text, kLabelFontSize, info_width, max_lines)) {
                        drawLine(line, color);
                    }
                };

                drawWrapped(std::format("fs {} Hz  |  y {:+.0f} .. {:+.0f} dB  |  x log {:.0f} .. {:.0f} Hz",
                                        static_cast<int>(kFs), kMagTopDb, kMagBottomDb, kFreqMin, kFreqMax),
                            LIGHTGRAY, 2);
                drawLine(std::format("kind: {} ({}-pole)", info.name, info.poles), RAYWHITE);
                // 每种类型的说明都实测在 info_width 内(最长的两条约 340px), 所以限一行
                drawWrapped(kindHint(info), RAYWHITE, 1);
                drawLine(std::format("RBJ         : {}",
                                     deviationText(g_design.has_rbj, g_curve.rbj_dev_db, g_curve.rbj_dev_hz)),
                         g_design.has_rbj ? kRbjColor : DARKGRAY);
                drawLine(std::format("MatchBiquad : {}",
                                     deviationText(g_design.has_mb, g_curve.mb_dev_db, g_curve.mb_dev_hz)),
                         g_design.has_mb ? kMbColor : DARKGRAY);
                drawLine(std::format("Ivantsov    : {}",
                                     deviationText(g_design.has_iv, g_curve.iv_dev_db, g_curve.iv_dev_hz)),
                         g_design.has_iv ? kIvColor : DARKGRAY);

                drawLegend(Rectangle{752.0f, static_cast<float>(knob_y), 200.0f, 4.0f * 17.0f + 8.0f}, g_design);
            }

            // ----- 脚注: 频率轴约定 + 偏差定义 + 操作提示 -----
            DrawText(
                "ref = analog prototype on the DIRECT axis (wa = 2*pi*f, no prewarp). MatchBiquad and Ivantsov "
                "target this axis; RBJ targets the PREWARPED axis (wa = 2*fs*tan(pi*f/fs)) and so must leave the "
                "ref at high frequencies.",
                10, 574, kLabelFontSize, GRAY);
            DrawText(
                "Ivantsov's match is a fit, not an exact map: its effective analog frequency reaches about 1.1*f "
                "near 20 kHz (sigma shapes that warping). That is a frequency-axis convention issue, not a design "
                "error.",
                10, 590, kLabelFontSize, GRAY);
            DrawText(
                "dev = max |designer - ref| over 20 Hz..20 kHz, counted only where ref > -60 dB.    "
                "drag knob to change / right click knob: reset / click a button to switch kind",
                10, 606, kLabelFontSize, GRAY);
        }
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
