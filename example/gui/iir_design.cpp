// ------------------------------------------------------------
// IIR 设计曲线查看器
//
// 用 qwqdsp_filter::IIRDesign / IIRDesignExtra / IIRDesignShelf 设计 IIR 滤波器,
// 画出幅度曲线与相位曲线, 并用竖线标出设计频率、横线标出设计电平
// (-3.01dB / 通带纹波 / 阻带电平 / shelf 的三个电平)。
//
// 所有曲线与标记都来自库的真实输出: 原型 -> 频率映射 -> 双线性变换 ->
// 双二阶系数 -> BiquadCoeff::DigitalResonpoce 连乘。
//
// 两个数字参数在不同原型下的含义(与库里的语义一致):
//
//   原型                      ripple 旋钮            atten 旋钮
//   ------------------------  ---------------------  -----------------------
//   IIRDesign::Butterworth    -                      -
//   IIRDesign::Chebyshev1     通带纹波 (0.01~6 dB)   -
//   IIRDesign::Chebyshev2     -                      阻带衰减 (10~120 dB)
//   IIRDesign::Elliptic       通带纹波 (0.01~6 dB)   阻带衰减 (10~120 dB)
//   Extra::ButterworthAttenDb -                      截止处幅度 (0.5~60 dB)
//   Extra::Chebyshev1         通带纹波 (0.01~6 dB)   截止处幅度 (0.5~60 dB)
//   Extra::Chebyshev2         阻带涟漪 (10~120 dB)   截止处幅度 (0.5~60 dB)
//
// Extra 的 Chebyshev2 与 Chebyshev1 参数含义相反: ripple 是阻带等波纹深度,
// atten 是 (1)rad/sec 处的幅度; atten=3.01 时与 IIRDesign::Chebyshev2 形状相同。
//
// ResponseKind::HighShelf 走另一条路(IIRDesignShelf): 中点频率由 fc 决定,
// 提升量由 boost 旋钮给(可正可负), 低频端恒为 0dB、高频端恒为 boost dB、
// 中点恒为 boost/2 dB。此时原型只有 Butterworth / Chebyshev1 / Elliptic 有效
// (另有对应关系: bw+atten -> bw, cheby1 extra -> cheby1), 其余原型会提示无解。
// ------------------------------------------------------------
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <format>
#include <numbers>
#include <optional>
#include <span>
#include <string>

#include "raylib.h"
#include "slider.hpp"

#include "qwqdsp/filter/biquad_coeff.hpp"
#include "qwqdsp/filter/iir_design.hpp"
#include "qwqdsp/filter/iir_design_extra.hpp"
#include "qwqdsp/filter/iir_design_shelf.hpp"

// ------------------------------------------------------------
// 常量与布局
// ------------------------------------------------------------
static constexpr int kWidth = 960;
static constexpr int kHeight = 620;

/// 采样率, 决定频率轴上限与双线性变换
static constexpr double kFs = 48000.0;
static constexpr double kFreqMin = 20.0;
static constexpr double kFreqMax = 20000.0;

/// 幅度面板纵轴范围 (dB), 下界是最低标记电平再往下留出余量
static constexpr double kMagTopDb = 6.0;
static constexpr double kMagMinFloorDb = -120.0;
static constexpr double kMagFloorMarginDb = 12.0;

/// 最多 8 个极点对; 带通/带阻映射后节数翻倍
static constexpr size_t kMaxPairs = 8;
static constexpr size_t kMaxSections = 2 * kMaxPairs;

static constexpr float kPlotLeft = 70.0f;
static constexpr float kPlotWidth = 860.0f;
static constexpr float kPlotRight = kPlotLeft + kPlotWidth;
static constexpr int kPlotColumns = static_cast<int>(kPlotWidth);

/// 每个像素列内的子采样数, 防止窄陷波被漏画
static constexpr int kSubSamples = 8;

static constexpr Rectangle kMagPanel{kPlotLeft, 148.0f, kPlotWidth, 240.0f};
static constexpr Rectangle kPhasePanel{kPlotLeft, 418.0f, kPlotWidth, 150.0f};

static constexpr int kSelectorFontSize = 12;
static constexpr int kLabelFontSize = 12;

// ------------------------------------------------------------
// 设计参数与结果
// ------------------------------------------------------------
enum class Prototype : int {
    Butterworth = 0,
    Chebyshev1,
    Chebyshev2,
    Elliptic,
    ButterworthAtten,
    ExtraChebyshev1,
    ExtraChebyshev2,
};
static constexpr std::array<char const*, 7> kPrototypeNames{
    "butterworth",
    "cheby1",
    "cheby2",
    "elliptic",
    "bw + atten",
    "cheby1 extra",
    "cheby2 extra",
};

enum class ResponseKind : int {
    Lowpass = 0,
    Highpass,
    Bandpass,
    Bandstop,
    HighShelf,
};
static constexpr std::array<char const*, 5> kKindNames{"lp", "hp", "bp", "bs", "h shelf"};
static constexpr std::array<char const*, 2> kEvenNames{"even off", "even on"};

struct DesignParams {
    Prototype prototype = Prototype::Butterworth;
    ResponseKind kind = ResponseKind::Lowpass;
    int num_pairs = 3;
    double fc = 1000.0;
    double q = 1.0;
    bool even_modify = false;
    // 按"含义"分别存值: 切换原型时同名旋钮会接到对应的值上
    double passband_ripple_db = 1.0;   // 通带纹波 (cheby1 / elliptic / cheby1 extra)
    double stopband_ripple_db = 40.0;  // 阻带涟漪 (cheby2 extra)
    double stopband_atten_db = 40.0;   // 阻带衰减 (cheby2 / elliptic)
    double cutoff_level_db = 3.0;      // 截止频率处的幅度 (extra 的三个原型)
    double boost_db = 12.0;            // shelf 的提升量(负值即衰减)
};

/// 设计出来的标记线: 竖线是设计频率, 横线是设计电平
struct MarkerLines {
    std::array<double, 3> vertical_hz{};
    size_t num_vertical = 0;
    std::array<double, 3> level_db{};
    size_t num_level = 0;
    size_t main_level_index = 0;
};

struct FilterDesign {
    std::array<qwqdsp_filter::BiquadCoeff, kMaxSections> biquad{};
    size_t num_sections = 0;
    MarkerLines markers{};
    bool valid = false;
    double mag_floor_db = kMagMinFloorDb;
    double mag_top_db = kMagTopDb;
};

/// 缓存下来的曲线: 每列一条竖直幅度范围 + 一个相位值
struct CurveCache {
    std::array<float, kPlotColumns> db_low{};
    std::array<float, kPlotColumns> db_high{};
    std::array<float, kPlotColumns> phase_deg{};
};

static DesignParams g_params{};
static FilterDesign g_design{};
static CurveCache g_curve{};
static bool g_dirty = true;

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
 * @brief 把幅度映射到幅度面板纵坐标
 * @param db 幅度 (dB)
 * @param floor_db 纵轴下界 (dB)
 * @param top_db 纵轴上界 (dB)
 * @return 面板内的 y 坐标
 */
static float dbToY(double db, double floor_db, double top_db) noexcept {
    double const clamped = std::clamp(db, floor_db, top_db);
    double const t = (top_db - clamped) / (top_db - floor_db);
    return kMagPanel.y + static_cast<float>(t) * kMagPanel.height;
}

/**
 * @brief 把相位映射到相位面板纵坐标
 * @param phase_deg 相位 (度), 卷绕在 [-180, 180]
 * @return 面板内的 y 坐标
 */
static float phaseToY(double phase_deg) noexcept {
    double const t = (180.0 - std::clamp(phase_deg, -180.0, 180.0)) / 360.0;
    return kPhasePanel.y + static_cast<float>(t) * kPhasePanel.height;
}

// ------------------------------------------------------------
// 设计与频响
// ------------------------------------------------------------

/// ripple 旋钮在不同原型下的含义
enum class RippleRole {
    Unused,
    PassbandRipple,  ///< 通带纹波, 通带边沿电平
    StopbandRipple,  ///< 阻带涟漪, 阻带等波纹深度
};

/// atten 旋钮在不同原型下的含义
enum class LevelRole {
    Unused,
    StopbandAtten,  ///< 阻带衰减指标
    CutoffLevel,    ///< 截止频率(1rad/sec)处的幅度
};

/**
 * @brief 取 ripple 旋钮在当前原型下的含义
 * @param prototype 原型类型
 * @return 对应的角色
 */
static RippleRole rippleRoleOf(Prototype prototype) noexcept {
    switch (prototype) {
    case Prototype::Chebyshev1:
    case Prototype::Elliptic:
    case Prototype::ExtraChebyshev1:
        return RippleRole::PassbandRipple;
    case Prototype::ExtraChebyshev2:
        return RippleRole::StopbandRipple;
    default:
        return RippleRole::Unused;
    }
}

/**
 * @brief 取 atten 旋钮在当前原型下的含义
 * @param prototype 原型类型
 * @return 对应的角色
 */
static LevelRole levelRoleOf(Prototype prototype) noexcept {
    switch (prototype) {
    case Prototype::Chebyshev2:
    case Prototype::Elliptic:
        return LevelRole::StopbandAtten;
    case Prototype::ButterworthAtten:
    case Prototype::ExtraChebyshev1:
    case Prototype::ExtraChebyshev2:
        return LevelRole::CutoffLevel;
    default:
        return LevelRole::Unused;
    }
}

/**
 * @brief 判断响应类型是不是 shelf
 * @param kind 响应类型
 * @return 是 shelf 返回 true
 * @note shelf 只在低频/高频两端定增益, 没有高通/带通那种频率映射
 */
static bool isShelfKind(ResponseKind kind) noexcept {
    return kind == ResponseKind::HighShelf;
}

/**
 * @brief 判断原型在 shelf 模式下是否有对应的构造
 *
 * shelf 走 IIRDesignShelf, 已实现巴特沃斯/切比雪夫I/切比雪夫II/椭圆四种。
 * Extra 的两个原型与基础原型同属一族, 这里一并接受(参数含义按基础原型取),
 * 这样在原型的各档之间切换时不必先改响应类型。
 *
 * @param prototype 原型类型
 * @return 有 shelf 构造返回 true
 */
static bool isShelfPrototype(Prototype prototype) noexcept {
    switch (prototype) {
    case Prototype::Butterworth:
    case Prototype::ButterworthAtten:
    case Prototype::Chebyshev1:
    case Prototype::ExtraChebyshev1:
    case Prototype::Chebyshev2:
    case Prototype::ExtraChebyshev2:
    case Prototype::Elliptic:
        return true;
    default:
        return false;
    }
}

/**
 * @brief 判断 shelf 模式下这个原型的纹波旋钮是不是"棚带等波纹"
 *
 * Chebyshev2 shelf 用 ripple 旋钮控制棚带的等波纹深度(单位 dB, 会被夹到
 * (0, |boost|) 内); Chebyshev1 shelf 用它控制通带纹波; 其余原型不用这个旋钮。
 *
 * @param prototype 原型类型
 * @return 是 Chebyshev2 一族返回 true
 */
static bool usesShelfBandRipple(Prototype prototype) noexcept {
    return prototype == Prototype::Chebyshev2 || prototype == Prototype::ExtraChebyshev2;
}

/**
 * @brief 判断原型是否有偶数极点修正开关
 * @param prototype 原型类型
 * @return 有则返回 true
 * @note shelf 用的是不带修正的原型, 所以 shelf 模式下这个开关不显示
 */
static bool usesEvenModify(Prototype prototype) noexcept {
    return prototype == Prototype::Chebyshev1
        || prototype == Prototype::Chebyshev2
        || prototype == Prototype::ExtraChebyshev1
        || prototype == Prototype::ExtraChebyshev2
        || prototype == Prototype::Elliptic;
}

/**
 * @brief 取竖线所在频率处的设计电平
 *
 * 该电平就是主横线的位置, 竖线与横线的交点即设计意图所在。
 * Butterworth / Chebyshev2 是 -3.01dB, Chebyshev1 / Elliptic 是通带边沿电平
 * -passband_ripple dB, Extra 的三个原型是截止频率处的幅度 -cutoff_level dB,
 * shelf 是中点的 boost/2 dB。
 *
 * @param params 设计参数
 * @return 设计频率处的目标电平 (dB)
 */
static double edgeLevelDb(DesignParams const& params) noexcept {
    if (isShelfKind(params.kind)) {
        return 0.5 * params.boost_db;
    }
    switch (params.prototype) {
    case Prototype::Butterworth:
    case Prototype::Chebyshev2:
        return -10.0 * std::log10(2.0);
    case Prototype::Chebyshev1:
    case Prototype::Elliptic:
        return -params.passband_ripple_db;
    default:
        return -params.cutoff_level_db;
    }
}

/**
 * @brief 取需要额外画的那条电平横线(通带边沿或阻带电平)
 * @param params 设计参数
 * @return 电平(dB); 与主横线重合或该原型不需要时返回空
 * @note Chebyshev2 / Elliptic 画阻带电平; Extra::Chebyshev1 画通带边沿;
 *       Extra::Chebyshev2 画阻带涟漪(注意它是阻带电平而不是通带电平)
 */
static std::optional<double> extraLevelDb(DesignParams const& params) noexcept {
    if (isShelfKind(params.kind)) {
        // shelf 的两端电平由 main_level 之外的循环补齐, 不走这里
        return std::nullopt;
    }
    switch (params.prototype) {
    case Prototype::Chebyshev2:
    case Prototype::Elliptic:
        return -params.stopband_atten_db;
    case Prototype::ExtraChebyshev1:
        return -params.passband_ripple_db;
    case Prototype::ExtraChebyshev2:
        return -params.stopband_ripple_db;
    default:
        return std::nullopt;
    }
}

/**
 * @brief 按参数设计滤波器, 并算好标记线
 * @param params 设计参数
 * @return 双二阶系数与标记信息, valid 为 false 表示参数无解
 */
static FilterDesign makeDesign(DesignParams const& params) {
    using qwqdsp_filter::BiquadCoeff;
    using qwqdsp_filter::IIRDesign;
    using qwqdsp_filter::IIRDesignExtra;
    using qwqdsp_filter::IIRDesignShelf;
    using ZPK = IIRDesign::ZPK;

    FilterDesign out;
    size_t const num_pairs = static_cast<size_t>(std::clamp(params.num_pairs, 1, static_cast<int>(kMaxPairs)));
    bool const is_shelf = isShelfKind(params.kind);
    bool const is_band = params.kind == ResponseKind::Bandpass || params.kind == ResponseKind::Bandstop;
    size_t const num_sections = (is_band ? 2 : 1) * num_pairs;
    out.num_sections = num_sections;

    // 值初始化: 零点默认在无穷远处 (原型赋值只覆盖用到的节)
    std::array<ZPK, kMaxSections> zpk{};
    std::span<ZPK> proto{zpk.data(), num_sections};

    // ----- 原型滤波器 -----
    bool prototype_ok = true;
    if (is_shelf) {
        // shelf 原型: 中点在 (1)rad/sec, 低频/中点/高频分别是 0 / boost/2 / boost dB
        if (!isShelfPrototype(params.prototype)) {
            prototype_ok = false;
        }
        else {
            switch (params.prototype) {
            case Prototype::Butterworth:
            case Prototype::ButterworthAtten:
                IIRDesignShelf::Butterworth(proto, num_pairs, params.boost_db);
                break;
            case Prototype::Chebyshev1:
            case Prototype::ExtraChebyshev1:
                IIRDesignShelf::Chebyshev1(proto, num_pairs, params.boost_db, params.passband_ripple_db);
                break;
            case Prototype::Chebyshev2:
            case Prototype::ExtraChebyshev2:
                prototype_ok = IIRDesignShelf::Chebyshev2(
                    proto, num_pairs, params.boost_db, params.passband_ripple_db);
                break;
            case Prototype::Elliptic:
                prototype_ok = IIRDesignShelf::Elliptic(
                    proto, num_pairs, params.boost_db, params.passband_ripple_db, params.stopband_atten_db);
                break;
            default:
                prototype_ok = false;
                break;
            }
        }
    }
    else {
        switch (params.prototype) {
        case Prototype::Butterworth:
            IIRDesign::Butterworth(proto, num_pairs);
            break;
        case Prototype::Chebyshev1:
            IIRDesign::Chebyshev1(proto, num_pairs, params.passband_ripple_db, params.even_modify);
            break;
        case Prototype::Chebyshev2:
            // IIRDesign::Chebyshev2 的 ripple 是负的阻带衰减
            IIRDesign::Chebyshev2(proto, num_pairs, -params.stopband_atten_db, params.even_modify);
            break;
        case Prototype::Elliptic:
            // 规格过陡时(阻带边沿与通带边沿重合)会返回 false, 此时不该继续算曲线
            prototype_ok = IIRDesign::Elliptic(
                proto, num_pairs, params.passband_ripple_db, params.stopband_atten_db, params.even_modify);
            break;
        case Prototype::ButterworthAtten:
            // 这个原型的 atten 是截止频率处的幅度
            IIRDesignExtra::ButterworthAttenDb(proto, num_pairs, params.cutoff_level_db);
            break;
        case Prototype::ExtraChebyshev1:
            // 通带纹波 + 截止处幅度; 库里有 atten >= ripple 的断言, 这里先夹住
            IIRDesignExtra::Chebyshev1(
                proto,
                num_pairs,
                params.passband_ripple_db,
                std::max(params.cutoff_level_db, params.passband_ripple_db),
                params.even_modify);
            break;
        case Prototype::ExtraChebyshev2:
            // 注意: 这个原型的 ripple 是阻带涟漪, atten 是截止处的幅度(与 Chebyshev1 相反)
            IIRDesignExtra::Chebyshev2(
                proto,
                num_pairs,
                params.stopband_ripple_db,
                params.cutoff_level_db,
                params.even_modify);
            break;
        }
    }

    // ----- 频率映射 -----
    if (!prototype_ok) {
        out.valid = false;
        return out;
    }
    double const wo = IIRDesign::Digital2AnalogW(params.fc, kFs);
    if (is_shelf) {
        // shelf 原型的"中点"就是 (1)rad/sec, 搬到 fc 即可
        IIRDesign::ProtyleToLowpass(proto, num_pairs, wo);
    }
    else {
        switch (params.kind) {
        case ResponseKind::Lowpass:
            IIRDesign::ProtyleToLowpass(proto, num_pairs, wo);
            break;
        case ResponseKind::Highpass:
            IIRDesign::ProtyleToHighpass(proto, num_pairs, wo);
            break;
        case ResponseKind::Bandpass:
            IIRDesign::ProtyleToBandpass(proto, num_pairs, wo, params.q);
            break;
        case ResponseKind::Bandstop:
            IIRDesign::ProtyleToBandstop(proto, num_pairs, wo, params.q);
            break;
        default:
            break;
        }
    }

    // ----- 离散化 -----
    IIRDesign::Bilinear(proto, kFs);

    std::span<BiquadCoeff> biquad{out.biquad.data(), num_sections};
    IIRDesign::TfToBiquad(proto, biquad);

    out.valid = true;
    for (size_t i = 0; i < num_pairs; ++i) {
        if (!std::isfinite(zpk[i].p.real()) || !std::isfinite(zpk[i].p.imag())) {
            out.valid = false;
        }
    }
    for (auto const& coeff : biquad) {
        if (!std::isfinite(coeff.b0) || !std::isfinite(coeff.b1) || !std::isfinite(coeff.b2)
            || !std::isfinite(coeff.a1) || !std::isfinite(coeff.a2)) {
            out.valid = false;
        }
    }

    // ----- 标记线 -----
    // 竖线: 原型通带边沿(u = 1rad/s)经映射后的数字频率。
    // LP/HP 就是设计频率本身, BP/BS 是 s^2 -/+ bw*s + wo^2 = 0 的两个根。
    if (is_band) {
        double const bw = wo / params.q;
        double const root = std::sqrt(bw * bw + 4.0 * wo * wo);
        double const w1 = (root - bw) * 0.5;
        double const w2 = (root + bw) * 0.5;
        // 逆双线性: omega = 2fs*tan(pi*f/fs) -> f = fs/pi*atan(omega/(2fs))
        out.markers.vertical_hz[0] = kFs / std::numbers::pi * std::atan(w1 / (2.0 * kFs));
        out.markers.vertical_hz[1] = kFs / std::numbers::pi * std::atan(w2 / (2.0 * kFs));
        out.markers.num_vertical = 2;
    }
    else {
        out.markers.vertical_hz[0] = params.fc;
        out.markers.num_vertical = 1;
    }

    // 横线: 设计频率处的电平永远画; 额外的通带/阻带电平与它不同时再补一条;
    // shelf 则画全三个参考电平(0 / boost/2 / boost)
    double const main_level = edgeLevelDb(params);
    out.markers.level_db[0] = main_level;
    out.markers.num_level = 1;
    out.markers.main_level_index = 0;
    if (is_shelf) {
        for (double const level : {0.0, params.boost_db}) {
            bool const duplicated = [&] {
                for (size_t i = 0; i < out.markers.num_level; ++i) {
                    if (std::abs(out.markers.level_db[i] - level) <= 0.01) {
                        return true;
                    }
                }
                return false;
            }();
            if (!duplicated && out.markers.num_level < out.markers.level_db.size()) {
                out.markers.level_db[out.markers.num_level++] = level;
            }
        }
    }
    else if (auto const extra = extraLevelDb(params); extra && std::abs(*extra - main_level) > 0.01) {
        out.markers.level_db[out.markers.num_level++] = *extra;
    }

    // 纵轴范围: 保证最低与最高的标记横线都看得见
    double lowest = main_level;
    double highest = main_level;
    for (size_t i = 0; i < out.markers.num_level; ++i) {
        lowest = std::min(lowest, out.markers.level_db[i]);
        highest = std::max(highest, out.markers.level_db[i]);
    }
    if (is_shelf) {
        // shelf 两端都是平的, 没有深阻带, 纵轴按台阶范围收紧, 否则台阶会被压成一条线
        out.mag_floor_db = lowest - kMagFloorMarginDb;
        out.mag_top_db = highest + kMagFloorMarginDb;
    }
    else {
        out.mag_floor_db = std::min(kMagMinFloorDb, lowest - kMagFloorMarginDb);
        out.mag_top_db = std::max(kMagTopDb, highest + kMagFloorMarginDb);
    }
    return out;
}

/**
 * @brief 计算双二阶串联在指定数字频率上的复频响
 * @param biquad 双二阶系数
 * @param w 数字角频率 (rad/sample)
 * @return 复频响 H(e^jw)
 */
static std::complex<double> frequencyResponse(
    std::span<qwqdsp_filter::BiquadCoeff const> biquad,
    double w) {
    std::complex<float> const z{static_cast<float>(std::cos(w)), static_cast<float>(std::sin(w))};
    std::complex<double> h{1.0, 0.0};
    for (auto const& coeff : biquad) {
        h *= static_cast<std::complex<double>>(coeff.DigitalResonpoce(z));
    }
    return h;
}

/**
 * @brief 逐像素列采样频响, 填入曲线缓存
 * @param design 设计结果
 * @param curve 输出的曲线缓存
 * @note 每列取多个子采样, 幅度用最小/最大值画成竖线, 相位取最后一个子采样
 */
static void sampleCurve(FilterDesign const& design, CurveCache& curve) {
    std::span<qwqdsp_filter::BiquadCoeff const> const biquad{
        design.biquad.data(), design.num_sections};
    constexpr double kVeryLowDb = -300.0;

    for (int column = 0; column < kPlotColumns; ++column) {
        double db_low = 1.0e30;
        double db_high = -1.0e30;
        double phase_deg = 0.0;
        for (int sub = 0; sub < kSubSamples; ++sub) {
            double const x = kPlotLeft + static_cast<double>(column)
                           + (static_cast<double>(sub) + 0.5) / kSubSamples;
            double const freq = xToFreq(x);
            double const w = 2.0 * std::numbers::pi * freq / kFs;
            auto const h = frequencyResponse(biquad, w);
            double const magnitude = std::abs(h);
            double db = magnitude > 0.0 ? 20.0 * std::log10(magnitude) : kVeryLowDb;
            db = std::max(db, kVeryLowDb);
            db_low = std::min(db_low, db);
            db_high = std::max(db_high, db);
            phase_deg = std::arg(h) * 180.0 / std::numbers::pi;
        }
        curve.db_low[column] = static_cast<float>(db_low);
        curve.db_high[column] = static_cast<float>(db_high);
        curve.phase_deg[column] = static_cast<float>(phase_deg);
    }
}

// ------------------------------------------------------------
// 绘制辅助
// ------------------------------------------------------------

/**
 * @brief 画一排可点击的选择按钮
 * @param bound 整排按钮占据的矩形
 * @param names 每项名称
 * @param selected 当前选中项下标
 * @param enabled 为 false 时画成灰色且不响应点击
 * @return 点击后应选中的项, 未点击则原样返回
 */
static size_t drawSelector(
    Rectangle bound,
    std::span<char const* const> names,
    size_t selected,
    bool enabled = true) {
    auto const mouse_pos = GetMousePosition();
    size_t result = selected;
    float const item_width = bound.width / static_cast<float>(names.size());

    for (size_t i = 0; i < names.size(); ++i) {
        Rectangle const item{
            bound.x + static_cast<float>(i) * item_width,
            bound.y,
            item_width - 2.0f,
            bound.height,
        };
        bool const active = enabled && (i == selected);
        bool const hover = enabled && CheckCollisionPointRec(mouse_pos, item);
        if (hover && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            result = i;
        }

        Color const fore = active ? BLACK : (hover ? WHITE : (enabled ? GRAY : DARKGRAY));
        if (active) {
            DrawRectangleRec(item, RAYWHITE);
        }
        else {
            DrawRectangleLinesEx(item, 1.0f, fore);
        }
        int const text_width = MeasureText(names[i], kSelectorFontSize);
        DrawText(
            names[i],
            static_cast<int>(item.x + (item.width - static_cast<float>(text_width)) * 0.5f),
            static_cast<int>(item.y + (item.height - static_cast<float>(kSelectorFontSize)) * 0.5f),
            kSelectorFontSize,
            fore);
    }
    return result;
}

/// 频率轴的刻度 (标数字的那些)
static constexpr std::array<double, 10> kFreqTicks{20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0, 20000.0};
/// 频率轴的细网格
static constexpr std::array<double, 16> kFreqGrid{
    20.0, 30.0, 50.0, 70.0, 100.0, 200.0, 300.0, 500.0, 700.0, 1000.0, 2000.0, 3000.0, 5000.0, 7000.0, 10000.0, 20000.0};

/**
 * @brief 格式化频率刻度文字
 * @param freq 频率 (Hz)
 * @return 例如 "1k" / "200"
 */
static std::string formatTick(double freq) {
    if (freq >= 1000.0) {
        return std::format("{:g}k", freq / 1000.0);
    }
    return std::format("{:g}", freq);
}

/**
 * @brief 画两个面板共用的频率网格与刻度文字
 */
static void drawFrequencyGrid() {
    for (double const freq : kFreqGrid) {
        float const x = freqToX(freq);
        DrawLineV({x, kMagPanel.y}, {x, kMagPanel.y + kMagPanel.height}, Color{40, 40, 40, 255});
        DrawLineV({x, kPhasePanel.y}, {x, kPhasePanel.y + kPhasePanel.height}, Color{40, 40, 40, 255});
    }
    for (double const freq : kFreqTicks) {
        float const x = freqToX(freq);
        DrawLineV({x, kMagPanel.y}, {x, kMagPanel.y + kMagPanel.height}, Color{70, 70, 70, 255});
        DrawLineV({x, kPhasePanel.y}, {x, kPhasePanel.y + kPhasePanel.height}, Color{70, 70, 70, 255});
        auto const label = formatTick(freq);
        int const label_width = MeasureText(label.c_str(), kLabelFontSize);
        DrawText(
            label.c_str(),
            static_cast<int>(x) - label_width / 2,
            static_cast<int>(kPhasePanel.y + kPhasePanel.height) + 4,
            kLabelFontSize,
            LIGHTGRAY);
    }
    DrawRectangleLinesEx(kMagPanel, 1.0f, GRAY);
    DrawRectangleLinesEx(kPhasePanel, 1.0f, GRAY);
}

/**
 * @brief 画幅度面板: 坐标标签、曲线、横线标记
 * @param design 设计结果
 * @param curve 曲线缓存
 */
static void drawMagnitudePanel(FilterDesign const& design, CurveCache const& curve) {
    double const floor_db = design.mag_floor_db;
    double const top_db = design.mag_top_db;
    DrawText("dB", 4, static_cast<int>(kMagPanel.y) - 2, kLabelFontSize, LIGHTGRAY);

    // 每 12dB 一条水平网格; 起点对齐到 12dB 的整数倍, 保证 0dB 线一定被画到
    double const grid_step = 12.0;
    for (double db = std::floor(top_db / grid_step) * grid_step; db >= floor_db; db -= grid_step) {
        float const y = dbToY(db, floor_db, top_db);
        DrawLineV({kPlotLeft, y}, {kPlotRight, y}, Color{45, 45, 45, 255});
        auto const label = std::format("{:g}", db);
        DrawText(label.c_str(), static_cast<int>(kMagPanel.x) - 42, static_cast<int>(y) - 7, kLabelFontSize, LIGHTGRAY);
    }

    // 曲线: 每列一条竖直线段, 相邻列自然连成曲线
    for (int column = 0; column < kPlotColumns; ++column) {
        float const x = kPlotLeft + static_cast<float>(column);
        float const y_high = dbToY(curve.db_high[column], floor_db, top_db) - 0.5f;
        float const y_low = dbToY(curve.db_low[column], floor_db, top_db) + 0.5f;
        DrawLineV({x, y_high}, {x, y_low}, GREEN);
    }
}

/**
 * @brief 画相位面板: 坐标标签与曲线
 * @param curve 曲线缓存
 */
static void drawPhasePanel(CurveCache const& curve) {
    DrawText("deg", 2, static_cast<int>(kPhasePanel.y) - 2, kLabelFontSize, LIGHTGRAY);
    for (int deg = 180; deg >= -180; deg -= 45) {
        float const y = phaseToY(static_cast<double>(deg));
        Color const color = (deg == 0) ? Color{90, 90, 90, 255} : Color{45, 45, 45, 255};
        DrawLineV({kPlotLeft, y}, {kPlotRight, y}, color);
        auto const label = std::format("{}", deg);
        DrawText(label.c_str(), static_cast<int>(kPhasePanel.x) - 42, static_cast<int>(y) - 7, kLabelFontSize, LIGHTGRAY);
    }

    // 相位曲线: 相邻列一律连线, 不做断线处理。
    // 陷波附近的相位一列内就能变化上百度和卷绕, 断开会让曲线出现大片空白;
    // 连起来以后卷绕处自然画成一条竖线, 与常见工具的相位图一致。
    float prev_x = kPlotLeft;
    float prev_y = phaseToY(curve.phase_deg[0]);
    DrawLineV({prev_x, prev_y - 0.5f}, {prev_x, prev_y + 0.5f}, ORANGE);
    for (int column = 1; column < kPlotColumns; ++column) {
        float const x = kPlotLeft + static_cast<float>(column);
        float const y = phaseToY(curve.phase_deg[column]);
        DrawLineV({prev_x, prev_y}, {x, y}, ORANGE);
        prev_x = x;
        prev_y = y;
    }
}

/**
 * @brief 画截止频率竖线与设计电平横线(只有线, 不含文字)
 * @param design 设计结果
 * @note 必须在曲线之前调用: shelf 的通带/棚带是平的, 曲线正好压在参考线上,
 *       参考线后画就会把曲线整段盖住
 */
static void drawMarkerLines(FilterDesign const& design) {
    double const floor_db = design.mag_floor_db;
    double const top_db = design.mag_top_db;

    for (size_t i = 0; i < design.markers.num_vertical; ++i) {
        float const x = freqToX(design.markers.vertical_hz[i]);
        DrawLineV({x, kMagPanel.y}, {x, kMagPanel.y + kMagPanel.height}, RED);
        DrawLineV({x, kPhasePanel.y}, {x, kPhasePanel.y + kPhasePanel.height}, RED);
    }
    for (size_t i = 0; i < design.markers.num_level; ++i) {
        float const y = dbToY(design.markers.level_db[i], floor_db, top_db);
        Color const color = (i == design.markers.main_level_index) ? RED : BLUE;
        DrawLineV({kPlotLeft, y}, {kPlotRight, y}, color);
    }
}

/**
 * @brief 画标记线的文字与设计点圆圈
 * @param design 设计结果
 * @note 在曲线之后调用, 保证文字不被曲线压住;
 *       竖线与主横线的交点画一个小圈, 表示设计意图所在的位置
 */
static void drawMarkerLabels(FilterDesign const& design) {
    double const floor_db = design.mag_floor_db;
    double const top_db = design.mag_top_db;

    for (size_t i = 0; i < design.markers.num_vertical; ++i) {
        float const x = freqToX(design.markers.vertical_hz[i]);
        auto const label = std::format("fc {:.1f} Hz", design.markers.vertical_hz[i]);
        int const label_width = MeasureText(label.c_str(), kLabelFontSize);
        DrawText(
            label.c_str(),
            static_cast<int>(x) - label_width / 2,
            static_cast<int>(kMagPanel.y + kMagPanel.height) + 4,
            kLabelFontSize,
            RED);
    }

    // 电平标签按行号从高到低排序, 相邻太近就把下面的往下推, 避免几条标签叠在一起。
    // 排序的是索引而不是坐标, 否则推完以后标签会贴到别的线上。
    std::array<size_t, 3> order{};
    for (size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    std::sort(order.begin(), order.begin() + design.markers.num_level, [&](size_t a, size_t b) {
        return design.markers.level_db[a] > design.markers.level_db[b];
    });

    float const label_gap = static_cast<float>(kLabelFontSize + 2);
    std::array<float, 3> label_y{};
    for (size_t i = 0; i < design.markers.num_level; ++i) {
        label_y[order[i]] = dbToY(design.markers.level_db[order[i]], floor_db, top_db) - kLabelFontSize - 2;
    }
    for (size_t i = 1; i < design.markers.num_level; ++i) {
        label_y[order[i]] = std::max(label_y[order[i]], label_y[order[i - 1]] + label_gap);
    }
    // 推下去以后可能超出面板, 从底部再往上挤一次
    float const label_bottom = kMagPanel.y + kMagPanel.height - static_cast<float>(kLabelFontSize);
    for (size_t i = design.markers.num_level; i-- > 1;) {
        label_y[order[i]] = std::min(label_y[order[i]], label_bottom);
        label_y[order[i - 1]] = std::min(label_y[order[i - 1]], label_y[order[i]] - label_gap);
    }

    for (size_t i = 0; i < design.markers.num_level; ++i) {
        double const level = design.markers.level_db[i];
        bool const is_main = (i == design.markers.main_level_index);
        Color const color = is_main ? RED : BLUE;

        auto const label = std::format("{:.2f} dB", level);
        int const label_width = MeasureText(label.c_str(), kLabelFontSize);
        DrawText(
            label.c_str(),
            static_cast<int>(kPlotRight) - label_width - 4,
            static_cast<int>(label_y[i]),
            kLabelFontSize,
            color);

        if (is_main) {
            // 设计频率处的实际曲线位置: 设计恰好达标时圆圈落在横线交点上,
            // 例如 IIRDesignExtra::Chebyshev2 的阻带优于指标时圆圈会低于横线
            std::span<qwqdsp_filter::BiquadCoeff const> const biquad{design.biquad.data(), design.num_sections};
            for (size_t v = 0; v < design.markers.num_vertical; ++v) {
                double const freq = design.markers.vertical_hz[v];
                auto const response = frequencyResponse(biquad, 2.0 * std::numbers::pi * freq / kFs);
                double const magnitude = std::abs(response);
                double const curve_db = magnitude > 0.0 ? 20.0 * std::log10(magnitude) : floor_db;
                float const x = freqToX(freq);
                DrawCircleLines(static_cast<int>(x), static_cast<int>(dbToY(curve_db, floor_db, top_db)), 5.0f, RED);
            }
        }
    }
}

// ------------------------------------------------------------
// 旋钮与参数角色的绑定
// ------------------------------------------------------------

/**
 * @brief 把 ripple / atten 旋钮接到当前原型对应的参数角色上
 *
 * 切换原型时调用: 重设旋钮标题、量程与当前值。量程按角色取:
 * 通带纹波 0.01~6dB, 阻带涟漪与阻带衰减 10~120dB, 截止处幅度 0.5~60dB。
 *
 * @param params 设计参数, 值取自对应角色的字段
 * @param ripple_knob 纹波旋钮
 * @param level_knob 幅度/衰减旋钮
 */
static void bindParamKnobs(DesignParams& params, Knob& ripple_knob, Knob& level_knob) {
    // shelf 模式: ripple 旋钮的含义按 shelf 构造重新指定
    if (isShelfKind(params.kind)) {
        if (usesShelfBandRipple(params.prototype)) {
            ripple_knob.set_title("band ripple");
            ripple_knob.set_range(0.01f, 6.0f, 0.01f, static_cast<float>(params.passband_ripple_db));
            ripple_knob.SetEnable(true);
        }
        else if (isShelfPrototype(params.prototype)) {
            // Chebyshev1 shelf: 通带纹波; 巴特沃斯/椭圆不用这个旋钮(椭圆用下面的阻带电平)
            bool const uses_ripple = params.prototype == Prototype::Chebyshev1
                                  || params.prototype == Prototype::ExtraChebyshev1;
            if (uses_ripple) {
                ripple_knob.set_title("pb ripple");
                ripple_knob.set_range(0.01f, 6.0f, 0.01f, static_cast<float>(params.passband_ripple_db));
            }
            ripple_knob.SetEnable(uses_ripple);
        }
        else {
            ripple_knob.SetEnable(false);
        }
        // shelf 模式下 level 旋钮只有椭圆用(当阻带衰减)
        bool const ell = params.prototype == Prototype::Elliptic;
        if (ell) {
            level_knob.set_title("sb atten");
            level_knob.set_range(10.0f, 120.0f, 0.5f, static_cast<float>(params.stopband_atten_db));
        }
        level_knob.SetEnable(ell);
        return;
    }

    switch (rippleRoleOf(params.prototype)) {
    case RippleRole::PassbandRipple:
        ripple_knob.set_title("pb ripple");
        ripple_knob.set_range(0.01f, 6.0f, 0.01f, static_cast<float>(params.passband_ripple_db));
        ripple_knob.SetEnable(true);
        break;
    case RippleRole::StopbandRipple:
        ripple_knob.set_title("sb ripple");
        ripple_knob.set_range(10.0f, 120.0f, 0.5f, static_cast<float>(params.stopband_ripple_db));
        ripple_knob.SetEnable(true);
        break;
    case RippleRole::Unused:
        ripple_knob.SetEnable(false);
        break;
    }

    switch (levelRoleOf(params.prototype)) {
    case LevelRole::StopbandAtten:
        level_knob.set_title("sb atten");
        level_knob.set_range(10.0f, 120.0f, 0.5f, static_cast<float>(params.stopband_atten_db));
        level_knob.SetEnable(true);
        break;
    case LevelRole::CutoffLevel:
        level_knob.set_title("cut level");
        level_knob.set_range(0.5f, 60.0f, 0.1f, static_cast<float>(params.cutoff_level_db));
        level_knob.SetEnable(true);
        break;
    case LevelRole::Unused:
        level_knob.SetEnable(false);
        break;
    }
}

/**
 * @brief 把 boost 旋钮接到 shelf 的提升量上
 *
 * 只有 shelf 模式用它; 其他响应类型下旋钮置灰。
 *
 * @param params 设计参数
 * @param boost_knob 提升量旋钮
 */
static void bindBoostKnob(DesignParams& params, Knob& boost_knob) {
    bool const enabled = isShelfKind(params.kind) && isShelfPrototype(params.prototype);
    boost_knob.SetEnable(enabled);
    if (enabled) {
        boost_knob.set_title("boost");
        boost_knob.set_range(-40.0f, 40.0f, 0.1f, static_cast<float>(params.boost_db));
    }
}

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
int main() {
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(kWidth, kHeight, "IIR design viewer");
    SetTargetFPS(60);

    // ----- 参数名称行 -----
    Rectangle const proto_bound{10.0f, 6.0f, 940.0f, 26.0f};
    Rectangle const kind_bound{10.0f, 36.0f, 530.0f, 26.0f};
    Rectangle const even_bound{548.0f, 36.0f, 220.0f, 26.0f};

    // ----- 数字参数旋钮 -----
    float const knob_width = 104.0f;
    float const knob_height = 74.0f;
    float const knob_y = 68.0f;
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
    };

    Knob pairs_knob;
    setupKnob(pairs_knob, knobBound(0), "pairs");
    pairs_knob.set_range(1.0f, static_cast<float>(kMaxPairs), 1.0f, static_cast<float>(g_params.num_pairs));
    pairs_knob.value_to_text_function = [](float value) { return std::format("order {:.0f}", 2.0f * value); };
    pairs_knob.on_value_change = [](float value) {
        g_params.num_pairs = static_cast<int>(value);
        g_dirty = true;
    };

    Knob fc_knob;
    setupKnob(fc_knob, knobBound(1), "fc");
    float const fc_log_max = static_cast<float>(std::log2(kFreqMax / kFreqMin));
    fc_knob.set_range(0.0f, fc_log_max, 0.001f, static_cast<float>(std::log2(g_params.fc / kFreqMin)));
    fc_knob.value_to_text_function = [](float value) {
        return std::format("{:.1f} Hz", kFreqMin * std::exp2(static_cast<double>(value)));
    };
    fc_knob.on_value_change = [](float value) {
        g_params.fc = kFreqMin * std::exp2(static_cast<double>(value));
        g_dirty = true;
    };

    Knob ripple_knob;
    setupKnob(ripple_knob, knobBound(2), "pb ripple");
    ripple_knob.value_to_text_function = [](float value) { return std::format("{:.2f} dB", value); };
    ripple_knob.on_value_change = [](float value) {
        // shelf 的 Chebyshev1/Chebyshev2 都用 passband_ripple_db 这个字段存纹波
        if (isShelfKind(g_params.kind)
            && (rippleRoleOf(g_params.prototype) == RippleRole::PassbandRipple
                || usesShelfBandRipple(g_params.prototype))) {
            g_params.passband_ripple_db = value;
            g_dirty = true;
            return;
        }
        switch (rippleRoleOf(g_params.prototype)) {
        case RippleRole::PassbandRipple:
            g_params.passband_ripple_db = value;
            break;
        case RippleRole::StopbandRipple:
            g_params.stopband_ripple_db = value;
            break;
        case RippleRole::Unused:
            break;
        }
        g_dirty = true;
    };

    Knob level_knob;
    setupKnob(level_knob, knobBound(3), "cut level");
    level_knob.value_to_text_function = [](float value) { return std::format("{:.1f} dB", value); };
    level_knob.on_value_change = [](float value) {
        switch (levelRoleOf(g_params.prototype)) {
        case LevelRole::StopbandAtten:
            g_params.stopband_atten_db = value;
            break;
        case LevelRole::CutoffLevel:
            g_params.cutoff_level_db = value;
            break;
        case LevelRole::Unused:
            break;
        }
        g_dirty = true;
    };

    Knob q_knob;
    setupKnob(q_knob, knobBound(4), "q");
    q_knob.set_range(0.3f, 10.0f, 0.01f, static_cast<float>(g_params.q));
    q_knob.value_to_text_function = [](float value) { return std::format("{:.2f}", value); };
    q_knob.on_value_change = [](float value) {
        g_params.q = value;
        g_dirty = true;
    };

    Knob boost_knob;
    setupKnob(boost_knob, knobBound(5), "boost");
    boost_knob.value_to_text_function = [](float value) { return std::format("{:+.1f} dB", value); };
    boost_knob.on_value_change = [](float value) {
        g_params.boost_db = value;
        g_dirty = true;
    };

    // 初始把两个参数旋钮接到当前原型的角色上
    bindParamKnobs(g_params, ripple_knob, level_knob);
    bindBoostKnob(g_params, boost_knob);

    while (!WindowShouldClose()) {
        BeginDrawing();
        {
            ClearBackground(BLACK);

            // ----- 选择行 -----
            size_t const prototype_index = drawSelector(proto_bound, kPrototypeNames, static_cast<size_t>(g_params.prototype));
            if (prototype_index != static_cast<size_t>(g_params.prototype)) {
                g_params.prototype = static_cast<Prototype>(prototype_index);
                // 两个数字参数的含义随原型变化, 重新把旋钮接到对应角色上
                bindParamKnobs(g_params, ripple_knob, level_knob);
                bindBoostKnob(g_params, boost_knob);
                g_dirty = true;
            }

            size_t const kind_index = drawSelector(kind_bound, kKindNames, static_cast<size_t>(g_params.kind));
            if (kind_index != static_cast<size_t>(g_params.kind)) {
                g_params.kind = static_cast<ResponseKind>(kind_index);
                // 进入/离开 shelf 时 boost 旋钮的可用性会变
                bindBoostKnob(g_params, boost_knob);
                g_dirty = true;
            }

            bool const even_enabled = usesEvenModify(g_params.prototype) && !isShelfKind(g_params.kind);
            size_t const even_index = drawSelector(
                even_bound,
                kEvenNames,
                g_params.even_modify ? 1 : 0,
                even_enabled);
            bool const even_modify = (even_index == 1);
            if (even_modify != g_params.even_modify) {
                g_params.even_modify = even_modify;
                g_dirty = true;
            }

            // ----- 旋钮 -----
            bool const is_band = g_params.kind == ResponseKind::Bandpass || g_params.kind == ResponseKind::Bandstop;
            q_knob.SetEnable(is_band);
            pairs_knob.display();
            fc_knob.display();
            ripple_knob.display();
            level_knob.display();
            q_knob.display();
            boost_knob.display();

            // ----- 参数变了就重新设计 -----
            if (g_dirty || g_params.fc < kFreqMin || g_params.fc > kFreqMax || g_params.q <= 0.0) {
                g_params.fc = std::clamp(g_params.fc, kFreqMin, kFreqMax);
                g_params.q = std::max(g_params.q, 0.05);
                g_design = makeDesign(g_params);
                if (g_design.valid) {
                    sampleCurve(g_design, g_curve);
                }
                g_dirty = false;
            }

            // ----- 曲线 -----
            drawFrequencyGrid();
            if (g_design.valid) {
                // 参考线先画, 曲线后画: shelf 的平坦段正好落在参考线上, 否则会被盖住
                drawMarkerLines(g_design);
                drawMagnitudePanel(g_design, g_curve);
                drawPhasePanel(g_curve);
                drawMarkerLabels(g_design);
            }
            else if (isShelfKind(g_params.kind) && !isShelfPrototype(g_params.prototype)) {
                DrawText(
                    "this prototype has no shelf construction; try butterworth / cheby1 / cheby2 / elliptic",
                    static_cast<int>(kPlotLeft),
                    static_cast<int>(kMagPanel.y) + 20,
                    20,
                    RED);
            }
            else {
                DrawText("parameter set has no finite solution", static_cast<int>(kPlotLeft), static_cast<int>(kMagPanel.y) + 20, 20, RED);
            }

            // ----- 文字信息 -----
            size_t const order = 2 * g_design.num_sections;
            auto const info = std::format(
                "order {}  |  {} biquad  |  fs {} Hz  |  markers: {} (vertical) + design level (horizontal)",
                order,
                g_design.num_sections,
                static_cast<int>(kFs),
                isShelfKind(g_params.kind) ? "midpoint fc" : "cutoff fc");
            DrawText(info.c_str(), 10, 588, kLabelFontSize, LIGHTGRAY);

            auto const used = std::format(
                "used: pairs, fc{}{}{}{}{}{}{}",
                isShelfKind(g_params.kind)
                    ? ((g_params.prototype == Prototype::Chebyshev1
                        || g_params.prototype == Prototype::ExtraChebyshev1
                        || usesShelfBandRipple(g_params.prototype)) ? ", band ripple" : "")
                    : (rippleRoleOf(g_params.prototype) == RippleRole::PassbandRipple ? ", passband ripple" : ""),
                isShelfKind(g_params.kind)
                    ? (g_params.prototype == Prototype::Elliptic ? ", stopband atten" : "")
                    : (rippleRoleOf(g_params.prototype) == RippleRole::StopbandRipple ? ", stopband ripple" : ""),
                levelRoleOf(g_params.prototype) == LevelRole::StopbandAtten && !isShelfKind(g_params.kind)
                    ? ", stopband atten" : "",
                levelRoleOf(g_params.prototype) == LevelRole::CutoffLevel && !isShelfKind(g_params.kind)
                    ? ", cutoff level" : "",
                is_band ? ", q" : "",
                usesEvenModify(g_params.prototype) && !isShelfKind(g_params.kind) ? ", even modify" : "",
                isShelfKind(g_params.kind) && isShelfPrototype(g_params.prototype) ? ", boost" : "");
            DrawText(used.c_str(), 10, 606, kLabelFontSize, DARKGRAY);
            DrawText("drag knob / right click: reset", 640, 606, kLabelFontSize, DARKGRAY);
        }
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
