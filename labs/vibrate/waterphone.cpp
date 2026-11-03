#include <raylib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <numbers>
#include <string>

#include "wav_export.hpp"

// 水琴：刚性金属圆盘边缘固定 24 根长短不一的钢制悬臂杆。
// 每根杆取 4 阶弯曲模态；点击是敲击冲量，按住拖动是弓弦摩擦（stick-slip）。
namespace qwqdsp_vibrate {

constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kTwoPi = 2.0f * kPi;

constexpr int kRodCount = 24;
constexpr int kModesPerRod = 4;
constexpr float kDiskRadius = 0.55f;
constexpr float kDiskThickness = 0.035f;
constexpr float kRodRadius = 0.004f;  // m
constexpr float kSteelDensity = 7850.0f;
constexpr float kYoungsModulus = 200.0e9f;
constexpr float kRodArea = kPi * kRodRadius * kRodRadius;
constexpr float kRodInertia =
    0.25f * kPi * kRodRadius * kRodRadius * kRodRadius * kRodRadius;
const float kRodWave =
    std::sqrt(kYoungsModulus * kRodInertia / (kSteelDensity * kRodArea));

constexpr std::array<float, kRodCount> kRodLengths = {
    0.46f, 0.31f, 0.19f, 0.37f, 0.24f, 0.52f, 0.28f, 0.21f,
    0.41f, 0.18f, 0.33f, 0.49f, 0.26f, 0.23f, 0.44f, 0.30f,
    0.20f, 0.38f, 0.54f, 0.25f, 0.35f, 0.22f, 0.47f, 0.29f};

constexpr std::array<float, kModesPerRod> kBeta = {
    1.87510407f, 4.69409113f, 7.85475744f, 10.99554073f};
// 悬臂模态形状中的 sigma，形状按 int(phi^2) dx = 1 归一化。
constexpr std::array<float, kModesPerRod> kSigma = {
    0.73409551f, 1.01846732f, 0.99922450f, 1.00003355f};

constexpr int kStepsPerSecond = 16000;
constexpr float kFixedDt = 1.0f / static_cast<float>(kStepsPerSecond);
constexpr float kMaxFrameTime = 0.05f;
constexpr int kMaxStepsPerFrame = 900;
constexpr float kSpeedMin = 0.10f;
constexpr float kSpeedMax = 1.00f;
constexpr float kSpeedDefault = 0.25f;

constexpr float kStrikeVelocity = 3.0f;  // m/s
constexpr float kStrikeWidth = 0.025f;   // m
constexpr float kStrikeContact = 0.72f;
constexpr float kInitialContact = 0.78f;
constexpr int kInitialRod = 0;
constexpr float kDefaultDamping = 0.08f;

constexpr float kBowNormal = 14.0f;     // N
constexpr float kBowMu = 0.75f;
constexpr float kBowSpeedMax = 2.0f;    // m/s
constexpr float kBowContactMin = 0.10f;
constexpr float kBowContactMax = 0.92f;

constexpr float kBendVisualGain = 16.0f;
constexpr int kRenderSegments = 12;

constexpr int kHistorySeconds = 4;
constexpr int kHistoryRate = 500;
constexpr int kHistorySize = 2048;
constexpr int kHistoryWindowSamples = kHistorySeconds * kHistoryRate;
constexpr int kStepsPerSample = kStepsPerSecond / kHistoryRate;

float ModeShape(int mode, float x) {
    const float b = kBeta[mode];
    return std::cosh(b * x) - std::cos(b * x) -
           kSigma[mode] * (std::sinh(b * x) - std::sin(b * x));
}

float RodFundamental(int rod) {
    const float length = kRodLengths[rod];
    return kBeta[0] * kBeta[0] * kRodWave / (kTwoPi * length * length);
}

struct ModeCoefficient {
    float omega = 0.0f;
    float omega_square = 0.0f;
    float gamma = 0.0f;
    float omega_d = 1.0f;
    float decay = 1.0f;
    float cos_wd = 1.0f;
    float sin_wd = 0.0f;
    float inv_omega_square = 0.0f;
    float mass = 1.0f;
    float tip_shape = 0.0f;
};

class Waterphone {
public:
    Waterphone() {
        q_.fill(0.0f);
        qdot_.fill(0.0f);
        RebuildCoefficients();
    }

    void Reset() {
        q_.fill(0.0f);
        qdot_.fill(0.0f);
        time_ = 0.0f;
        ClearBow();
    }

    void SetDamping(float damping) {
        damping_ = std::max(0.0f, damping);
        RebuildCoefficients();
    }

    // 在 fraction 处给杆一个沿切向的初始速度冲量。
    void Strike(int rod, float fraction, float velocity) {
        if (rod < 0 || rod >= kRodCount || velocity == 0.0f) {
            return;
        }
        const float x = std::clamp(fraction, 0.02f, 0.98f);
        const float length = kRodLengths[rod];
        for (int n = 0; n < kModesPerRod; ++n) {
            const int index = rod * kModesPerRod + n;
            qdot_[index] += velocity * ModeShape(n, x) * kStrikeWidth / length;
        }
    }

    // 弓在 fraction 处以速度 speed 摩擦杆；normal 为法向压力。
    void SetBow(int rod, float fraction, float speed, float normal) {
        if (rod < 0 || rod >= kRodCount || normal <= 0.0f) {
            ClearBow();
            return;
        }
        bow_rod_ = rod;
        bow_fraction_ = std::clamp(fraction, kBowContactMin, kBowContactMax);
        bow_speed_ = std::clamp(speed, -kBowSpeedMax, kBowSpeedMax);
        bow_normal_ = normal;
    }

    void ClearBow() {
        bow_rod_ = -1;
        bow_speed_ = 0.0f;
        bow_normal_ = 0.0f;
    }

    void Step() {
        float bow_force = 0.0f;
        float bow_contact_fraction = 0.0f;
        float bow_phi[kModesPerRod]{};
        float bow_accel[kModesPerRod]{};
        if (bow_rod_ >= 0 && bow_normal_ > 0.0f) {
            bow_contact_fraction = bow_fraction_;
            float contact_velocity = 0.0f;
            float inverse_effective_mass = 0.0f;
            for (int n = 0; n < kModesPerRod; ++n) {
                const int index = bow_rod_ * kModesPerRod + n;
                const float phi = ModeShape(n, bow_contact_fraction);
                bow_phi[n] = phi;
                contact_velocity += qdot_[index] * phi;
                inverse_effective_mass +=
                    phi * phi / modes_[index].mass;
            }
            if (inverse_effective_mass > 1e-12f) {
                const float required = (bow_speed_ - contact_velocity) /
                                       (inverse_effective_mass * kFixedDt);
                const float limit = kBowMu * bow_normal_;
                bow_force = std::clamp(required, -limit, limit);
            }
            for (int n = 0; n < kModesPerRod; ++n) {
                const int index = bow_rod_ * kModesPerRod + n;
                bow_accel[n] =
                    bow_force * bow_phi[n] / modes_[index].mass;
            }
        }

        for (int rod = 0; rod < kRodCount; ++rod) {
            const bool bowed = rod == bow_rod_;
            for (int n = 0; n < kModesPerRod; ++n) {
                const int index = rod * kModesPerRod + n;
                Advance(index, bowed ? bow_accel[n] : 0.0f);
            }
        }
        time_ += kFixedDt;
    }

    float RodLateral(int rod, float fraction) const {
        if (rod < 0 || rod >= kRodCount) {
            return 0.0f;
        }
        const float x = std::clamp(fraction, 0.0f, 1.0f);
        float value = 0.0f;
        for (int n = 0; n < kModesPerRod; ++n) {
            value += q_[rod * kModesPerRod + n] * ModeShape(n, x);
        }
        return value;
    }

    float TipDisplacement(int rod) const {
        if (rod < 0 || rod >= kRodCount) {
            return 0.0f;
        }
        float value = 0.0f;
        for (int n = 0; n < kModesPerRod; ++n) {
            value += q_[rod * kModesPerRod + n] *
                     modes_[rod * kModesPerRod + n].tip_shape;
        }
        return value;
    }

    float TotalEnergy() const {
        float energy = 0.0f;
        for (int index = 0; index < kRodCount * kModesPerRod; ++index) {
            const ModeCoefficient& mode = modes_[index];
            energy += 0.5f * mode.mass *
                      (qdot_[index] * qdot_[index] +
                       mode.omega_square * q_[index] * q_[index]);
        }
        return energy;
    }

    float PeakTipDisplacement() const {
        float peak = 0.0f;
        for (int rod = 0; rod < kRodCount; ++rod) {
            peak = std::max(peak, std::abs(TipDisplacement(rod)));
        }
        return peak;
    }

    float Time() const { return time_; }
    int BowingRod() const { return bow_rod_; }
    float BowingFraction() const { return bow_fraction_; }
    bool IsBowing() const { return bow_rod_ >= 0 && bow_normal_ > 0.0f; }

    bool IsFinite() const {
        for (int index = 0; index < kRodCount * kModesPerRod; ++index) {
            if (!std::isfinite(q_[index]) || !std::isfinite(qdot_[index])) {
                return false;
            }
        }
        return true;
    }

private:
    void RebuildCoefficients() {
        for (int rod = 0; rod < kRodCount; ++rod) {
            const float length = kRodLengths[rod];
            const float length_square = length * length;
            for (int n = 0; n < kModesPerRod; ++n) {
                const int index = rod * kModesPerRod + n;
                ModeCoefficient& mode = modes_[index];
                mode.omega = kBeta[n] * kBeta[n] * kRodWave / length_square;
                mode.omega_square = mode.omega * mode.omega;
                mode.gamma = 0.5f * damping_;
                const float under = mode.omega_square -
                                    mode.gamma * mode.gamma;
                mode.omega_d = std::sqrt(std::max(under, 1e-6f));
                mode.decay = std::exp(-mode.gamma * kFixedDt);
                mode.cos_wd = std::cos(mode.omega_d * kFixedDt);
                mode.sin_wd = std::sin(mode.omega_d * kFixedDt);
                mode.inv_omega_square = 1.0f / mode.omega_square;
                mode.mass = kSteelDensity * kRodArea * length;
                mode.tip_shape = ModeShape(n, 1.0f);
            }
        }
    }

    // q'' + 2 gamma q' + omega^2 q = force，force 在步内视为常数。
    void Advance(int index, float force) {
        const ModeCoefficient& mode = modes_[index];
        const float equilibrium = force * mode.inv_omega_square;
        const float y0 = q_[index] - equilibrium;
        const float v0 = qdot_[index];
        const float b = (v0 + mode.gamma * y0) / mode.omega_d;
        const float y = mode.decay *
                        (y0 * mode.cos_wd + b * mode.sin_wd);
        const float v = mode.decay *
                        (v0 * mode.cos_wd -
                         ((mode.omega_square * y0 + mode.gamma * v0) /
                          mode.omega_d) * mode.sin_wd);
        q_[index] = equilibrium + y;
        qdot_[index] = v;
    }

    std::array<ModeCoefficient, kRodCount * kModesPerRod> modes_{};
    std::array<float, kRodCount * kModesPerRod> q_{};
    std::array<float, kRodCount * kModesPerRod> qdot_{};
    float damping_ = kDefaultDamping;
    float time_ = 0.0f;
    int bow_rod_ = -1;
    float bow_fraction_ = 0.0f;
    float bow_speed_ = 0.0f;
    float bow_normal_ = 0.0f;
};

class ProbeHistory {
public:
    void Clear() {
        head_ = 0;
        filled_ = 0;
        std::fill(samples_.begin(), samples_.end(), 0.0f);
    }

    void Push(float value) {
        samples_[head_] = value;
        head_ = (head_ + 1) % kHistorySize;
        filled_ = std::min(filled_ + 1, kHistorySize);
    }

    float At(int age) const {
        return samples_[(head_ - 1 - age + kHistorySize * 2) % kHistorySize];
    }

    int Size() const { return filled_; }

private:
    std::array<float, kHistorySize> samples_{};
    int head_ = 0;
    int filled_ = 0;
};

struct Layout {
    Rectangle scene{};
    Rectangle scope{};
    Vector2 center{};
    float axis_scale = 300.0f;
    float depth_x = 165.0f;
    float depth_y = -90.0f;
    float z_scale = 650.0f;
};

Layout MakeLayout(int width, int height) {
    Layout layout;
    const float w = static_cast<float>(width);
    const float h = static_cast<float>(height);
    const float scene_height = h * 0.63f;
    layout.scene = Rectangle{0.0f, 0.0f, w, scene_height};
    layout.center = Vector2{w * 0.5f, scene_height * 0.74f};
    layout.axis_scale = std::min(w * 0.26f, scene_height * 0.45f);
    layout.depth_x = layout.axis_scale * 0.55f;
    layout.depth_y = -layout.axis_scale * 0.30f;
    layout.z_scale = scene_height * 0.85f;
    const float scope_top = scene_height + 22.0f;
    layout.scope = Rectangle{w * 0.07f, scope_top, w * 0.86f,
                             h - 24.0f - scope_top};
    return layout;
}

Vector2 Project(const Layout& layout, float wx, float wy, float wz) {
    return Vector2{
        layout.center.x + wx * layout.axis_scale + wy * layout.depth_x,
        layout.center.y + wy * layout.depth_y - wz * layout.z_scale,
    };
}

Vector2 RodDisplayPoint(const Waterphone& model, const Layout& layout,
                        int rod, float fraction) {
    const float x = std::clamp(fraction, 0.0f, 1.0f);
    const float angle = kTwoPi * static_cast<float>(rod) /
                        static_cast<float>(kRodCount);
    const float cosine = std::cos(angle);
    const float sine = std::sin(angle);
    const float lateral = model.RodLateral(rod, x) * kBendVisualGain;
    return Project(layout,
                   kDiskRadius * cosine - sine * lateral,
                   kDiskRadius * sine + cosine * lateral,
                   kRodLengths[rod] * x);
}

Rectangle SpeedSliderTrack(int width) {
    return Rectangle{static_cast<float>(width) - 250.0f, 40.0f, 200.0f, 4.0f};
}

Rectangle SpeedSliderHit(int width) {
    return Rectangle{static_cast<float>(width) - 260.0f, 26.0f, 220.0f, 30.0f};
}

void DrawSpeedSlider(int width, float speed) {
    const Rectangle track = SpeedSliderTrack(width);
    const float fraction = (speed - kSpeedMin) / (kSpeedMax - kSpeedMin);
    const float knob_x = track.x + fraction * track.width;
    const std::string label = std::format("speed {:.2f}x", speed);
    DrawText(label.c_str(), static_cast<int>(track.x), 14, 13,
             Color{120, 140, 165, 230});
    DrawLineEx(Vector2{track.x, track.y},
               Vector2{track.x + track.width, track.y}, 2.0f,
               Color{42, 56, 74, 255});
    DrawLineEx(Vector2{track.x, track.y}, Vector2{knob_x, track.y}, 2.0f,
               Color{110, 180, 245, 255});
    DrawCircleV(Vector2{knob_x, track.y}, 5.0f, Color{190, 225, 255, 255});
}

Rectangle PadBox(int width) {
    return Rectangle{static_cast<float>(width) - 168.0f, 70.0f, 148.0f, 148.0f};
}

Rectangle PadInner(const Rectangle& box) {
    return Rectangle{box.x + 10.0f, box.y + 10.0f, box.width - 20.0f,
                     box.height - 20.0f};
}

int PadRod(const Rectangle& inner, Vector2 mouse) {
    const Vector2 center{inner.x + inner.width * 0.5f,
                         inner.y + inner.height * 0.5f};
    const float angle =
        std::atan2(center.y - mouse.y, mouse.x - center.x);
    const float wrapped = angle < 0.0f ? angle + kTwoPi : angle;
    int best = 0;
    float best_error = kTwoPi;
    for (int rod = 0; rod < kRodCount; ++rod) {
        const float rod_angle = kTwoPi * static_cast<float>(rod) /
                                static_cast<float>(kRodCount);
        float error = std::abs(wrapped - rod_angle);
        error = std::min(error, kTwoPi - error);
        if (error < best_error) {
            best_error = error;
            best = rod;
        }
    }
    return best;
}

void DrawPad(const Waterphone& model, const Rectangle& box, int probe_rod) {
    const Rectangle inner = PadInner(box);
    const Vector2 center{inner.x + inner.width * 0.5f,
                         inner.y + inner.height * 0.5f};
    const float rim = inner.width * 0.36f;
    DrawRectangleRec(box, Color{7, 9, 14, 255});
    DrawRectangleLinesEx(box, 1.0f, Color{44, 60, 82, 230});
    DrawText("rim", static_cast<int>(box.x) + 8,
             static_cast<int>(box.y) + 6, 12, Color{100, 120, 145, 220});

    DrawCircleLinesV(center, rim, Color{255, 225, 175, 160});
    DrawCircleLinesV(center, rim * 0.55f, Color{80, 100, 125, 100});
    for (int rod = 0; rod < kRodCount; ++rod) {
        const float angle = kTwoPi * static_cast<float>(rod) /
                            static_cast<float>(kRodCount);
        const Vector2 rim_point{center.x + std::cos(angle) * rim,
                                center.y - std::sin(angle) * rim};
        const float tick = 3.0f + 12.0f * kRodLengths[rod] / 0.54f;
        const Vector2 end{center.x + std::cos(angle) * (rim + tick),
                          center.y - std::sin(angle) * (rim + tick)};
        const bool selected = rod == probe_rod;
        const Color color = selected ? Color{255, 200, 130, 255}
                                     : Color{150, 175, 205, 150};
        DrawLineEx(rim_point, end, selected ? 2.0f : 1.0f, color);
        DrawCircleV(rim_point, selected ? 3.0f : 1.8f, color);
    }
    if (model.IsBowing()) {
        DrawCircleLinesV(center, rim + 4.0f, Color{255, 190, 120, 150});
    }
}

void DrawSceneFloor(const Layout& layout) {
    constexpr int kSegments = 72;
    for (int i = 0; i < kSegments; ++i) {
        const float a0 = kTwoPi * static_cast<float>(i) / kSegments;
        const float a1 = kTwoPi * static_cast<float>(i + 1) / kSegments;
        DrawLineV(Project(layout, kDiskRadius * 1.35f * std::cos(a0),
                          kDiskRadius * 1.35f * std::sin(a0),
                          -kDiskThickness),
                  Project(layout, kDiskRadius * 1.35f * std::cos(a1),
                          kDiskRadius * 1.35f * std::sin(a1),
                          -kDiskThickness),
                  Color{24, 32, 45, 110});
    }
}

void DrawDisk(const Layout& layout) {
    constexpr int kSegments = 96;
    std::array<Vector2, kSegments> top{};
    std::array<Vector2, kSegments> bottom{};
    for (int i = 0; i < kSegments; ++i) {
        const float angle = kTwoPi * static_cast<float>(i) / kSegments;
        top[i] = Project(layout, kDiskRadius * std::cos(angle),
                         kDiskRadius * std::sin(angle), 0.0f);
        bottom[i] = Project(layout, kDiskRadius * std::cos(angle),
                            kDiskRadius * std::sin(angle), -kDiskThickness);
    }
    for (int i = 0; i < kSegments; ++i) {
        const int next = (i + 1) % kSegments;
        const Color side{42, 50, 61, 255};
        DrawTriangle(top[i], top[next], bottom[next], side);
        DrawTriangle(top[i], bottom[next], bottom[i], side);
    }
    const Vector2 center = Project(layout, 0.0f, 0.0f, 0.0f);
    for (int i = 0; i < kSegments; ++i) {
        const int next = (i + 1) % kSegments;
        const float angle = kTwoPi * (static_cast<float>(i) + 0.5f) /
                            kSegments;
        const float sheen = 0.5f + 0.5f * std::cos(angle - 2.2f);
        DrawTriangle(center, top[i], top[next],
                     ColorLerp(Color{55, 64, 76, 255},
                               Color{102, 116, 134, 255}, sheen * 0.65f));
    }
    for (float radius : {0.32f, 0.58f, 0.82f}) {
        for (int i = 0; i < kSegments; ++i) {
            const int next = (i + 1) % kSegments;
            const float a0 = kTwoPi * static_cast<float>(i) / kSegments;
            const float a1 = kTwoPi * static_cast<float>(next) / kSegments;
            DrawLineV(Project(layout, kDiskRadius * radius * std::cos(a0),
                              kDiskRadius * radius * std::sin(a0), 0.0f),
                      Project(layout, kDiskRadius * radius * std::cos(a1),
                              kDiskRadius * radius * std::sin(a1), 0.0f),
                      Color{125, 145, 168, 45});
        }
    }
    for (int rod = 0; rod < kRodCount; ++rod) {
        const float angle = kTwoPi * static_cast<float>(rod) /
                            static_cast<float>(kRodCount);
        DrawLineEx(center,
                   Project(layout, kDiskRadius * std::cos(angle),
                           kDiskRadius * std::sin(angle), 0.0f),
                   1.0f, Color{90, 108, 128, 35});
    }
    for (int i = 0; i < kSegments; ++i) {
        const int next = (i + 1) % kSegments;
        DrawLineV(top[i], top[next], Color{145, 165, 188, 180});
    }
    DrawCircleV(center, 4.0f, Color{130, 150, 172, 180});
}

void DrawRod(const Waterphone& model, const Layout& layout, int rod,
             bool selected) {
    const Color main_color = selected ? Color{255, 205, 140, 255}
                                       : Color{190, 205, 225, 210};
    const Color glow_color = selected ? Color{255, 190, 120, 38}
                                       : Color{150, 175, 210, 20};
    Vector2 previous = RodDisplayPoint(model, layout, rod, 0.0f);
    const Vector2 base = previous;
    for (int pass = 0; pass < 2; ++pass) {
        previous = base;
        for (int i = 1; i <= kRenderSegments; ++i) {
            const float fraction = static_cast<float>(i) /
                                   static_cast<float>(kRenderSegments);
            const Vector2 point =
                RodDisplayPoint(model, layout, rod, fraction);
            DrawLineEx(previous, point, pass == 0 ? 5.0f : 1.7f,
                       pass == 0 ? glow_color : main_color);
            previous = point;
        }
    }
    const Vector2 tip = RodDisplayPoint(model, layout, rod, 1.0f);
    const Vector2 mid = RodDisplayPoint(model, layout, rod, 0.55f);
    DrawCircleV(base, selected ? 3.0f : 2.2f,
                selected ? Color{255, 205, 140, 255}
                         : Color{145, 165, 190, 190});
    DrawLineEx(mid, tip, 1.0f,
               selected ? Color{255, 235, 210, 120}
                        : Color{220, 232, 245, 40});
    DrawCircleV(tip, selected ? 3.5f : 2.5f, main_color);
}

struct RodHit {
    int rod = -1;
    float fraction = 0.0f;
    float distance = 1e9f;
};

RodHit HitTestRod(const Waterphone& model, const Layout& layout,
                  Vector2 mouse) {
    RodHit hit;
    for (int rod = 0; rod < kRodCount; ++rod) {
        Vector2 previous = RodDisplayPoint(model, layout, rod, 0.0f);
        for (int i = 1; i <= kRenderSegments; ++i) {
            const float fraction = static_cast<float>(i) /
                                   static_cast<float>(kRenderSegments);
            const Vector2 point =
                RodDisplayPoint(model, layout, rod, fraction);
            const float dx = point.x - previous.x;
            const float dy = point.y - previous.y;
            const float length_square = dx * dx + dy * dy;
            const float t = length_square <= 1e-12f
                                ? 0.0f
                                : std::clamp(
                                      ((mouse.x - previous.x) * dx +
                                       (mouse.y - previous.y) * dy) /
                                          length_square,
                                      0.0f, 1.0f);
            const float distance =
                std::hypot(mouse.x - (previous.x + dx * t),
                           mouse.y - (previous.y + dy * t));
            if (distance < hit.distance) {
                hit.distance = distance;
                hit.rod = rod;
                hit.fraction = (static_cast<float>(i - 1) + t) /
                               static_cast<float>(kRenderSegments);
            }
            previous = point;
        }
    }
    if (hit.distance > 11.0f) {
        hit.rod = -1;
    }
    return hit;
}

float BowSpeed(const Layout& layout, int rod, Vector2 mouse_delta,
               float frame_time) {
    const float angle = kTwoPi * static_cast<float>(rod) /
                        static_cast<float>(kRodCount);
    const float tx = -std::sin(angle);
    const float ty = std::cos(angle);
    Vector2 tangent{
        tx * layout.axis_scale + ty * layout.depth_x,
        ty * layout.depth_y,
    };
    const float length = std::hypot(tangent.x, tangent.y);
    if (length < 1e-4f) {
        return 0.0f;
    }
    tangent.x /= length;
    tangent.y /= length;
    const float along = mouse_delta.x * tangent.x + mouse_delta.y * tangent.y;
    const float dt = std::max(frame_time, 1.0f / 240.0f);
    return std::clamp(along / layout.axis_scale / dt, -kBowSpeedMax,
                      kBowSpeedMax);
}

void DrawSceneHud(const Waterphone& model) {
    const std::string line = std::format(
        "t {:.2f} s      click rod: strike      drag rod: bow      "
        "R: reset      E: export wav",
        model.Time());
    DrawText(line.c_str(), 24, 20, 15, Color{120, 140, 165, 230});
}

std::string ProbeTitle(int rod) {
    return std::format("PROBE   ROD {:02d}   L {:.2f} m   f1 {:.1f} Hz", rod,
                       kRodLengths[rod], RodFundamental(rod));
}

void DrawScope(const ProbeHistory& history, const Waterphone& model,
               int probe_rod, float& scale) {
    const Layout layout = MakeLayout(GetScreenWidth(), GetScreenHeight());
    const Rectangle panel = layout.scope;
    const Rectangle plot{panel.x + 18.0f, panel.y + 44.0f, panel.width - 36.0f,
                         panel.height - 76.0f};
    const float mid = plot.y + plot.height * 0.5f;
    const int count = std::min(history.Size(), kHistoryWindowSamples);
    const float value = model.TipDisplacement(probe_rod);

    float peak = 0.0f;
    for (int age = 0; age < count; ++age) {
        peak = std::max(peak, std::abs(history.At(age)));
    }
    if (count > 0) {
        const float target = std::max(peak * 1.15f, 0.0005f);
        if (target > scale) {
            scale = target;
        } else {
            scale += (target - scale) * 0.04f;
        }
    }

    DrawRectangleRec(panel, Color{7, 9, 14, 255});
    DrawRectangleLinesEx(panel, 1.0f, Color{42, 58, 80, 220});

    const std::string title = ProbeTitle(probe_rod);
    DrawText(title.c_str(), static_cast<int>(panel.x) + 16,
             static_cast<int>(panel.y) + 12, 15, Color{200, 215, 235, 255});

    const std::string live = std::format("u {:+.4f} m      scale {:.4f} m",
                                         value, scale);
    DrawText(live.c_str(),
             static_cast<int>(panel.x + panel.width) - 16 -
                 MeasureText(live.c_str(), 14),
             static_cast<int>(panel.y) + 13, 14, Color{130, 150, 175, 255});

    DrawRectangleRec(plot, Color{5, 8, 14, 220});
    DrawRectangleLinesEx(plot, 1.0f, Color{38, 56, 82, 200});

    for (float x = plot.x; x < plot.x + plot.width; x += 12.0f) {
        DrawLineV(Vector2{x, mid},
                  Vector2{std::min(x + 6.0f, plot.x + plot.width), mid},
                  Color{60, 80, 110, 120});
    }
    for (int k = 1; k <= 3; ++k) {
        const float x = plot.x + plot.width * static_cast<float>(k) / 4.0f;
        DrawLineV(Vector2{x, plot.y + 2.0f},
                  Vector2{x, plot.y + plot.height - 2.0f},
                  Color{36, 52, 74, 130});
    }

    const int label_y = static_cast<int>(plot.y + plot.height + 6.0f);
    for (int k = 1; k <= 3; ++k) {
        const std::string label = std::format("-{}s", kHistorySeconds - k);
        const int label_width = MeasureText(label.c_str(), 13);
        const float x = plot.x + plot.width * static_cast<float>(k) / 4.0f;
        DrawText(label.c_str(), static_cast<int>(x) - label_width / 2, label_y,
                 13, Color{100, 120, 145, 220});
    }
    DrawText("now",
             static_cast<int>(plot.x + plot.width - MeasureText("now", 13)),
             label_y, 13, Color{100, 120, 145, 220});

    if (count >= 2) {
        const Color trace_color{190, 220, 255, 235};
        Vector2 previous{};
        bool has_previous = false;
        for (int age = count - 1; age >= 0; --age) {
            const float x =
                plot.x + plot.width *
                             (1.0f - static_cast<float>(age) /
                                         static_cast<float>(
                                             kHistoryWindowSamples - 1));
            const float y = std::clamp(
                mid - history.At(age) / scale * (plot.height * 0.5f - 4.0f),
                plot.y + 1.0f, plot.y + plot.height - 1.0f);
            const Vector2 current{x, y};
            if (has_previous) {
                DrawLineEx(previous, current, 1.7f, trace_color);
            }
            previous = current;
            has_previous = true;
        }
    } else {
        const char* waiting = "click a rod to strike it...";
        DrawText(waiting,
                 static_cast<int>(plot.x + plot.width * 0.5f) -
                     MeasureText(waiting, 14) / 2,
                 static_cast<int>(mid) - 7, 14, Color{90, 110, 135, 255});
    }
}

int RunSelfTest() {
    int failures = 0;
    const auto check = [&failures](bool ok, const char* name,
                                   const std::string& detail) {
        std::cout << std::format("[{}] {:<28} {}\n", ok ? "PASS" : "FAIL",
                                 name, detail);
        if (!ok) {
            ++failures;
        }
    };

    {
        Waterphone model;
        model.SetDamping(0.0f);
        model.Reset();
        model.Strike(0, kInitialContact, kStrikeVelocity);
        const float energy0 = model.TotalEnergy();
        float peak = 0.0f;
        for (int step = 0; step < kStepsPerSecond; ++step) {
            model.Step();
            peak = std::max(peak, model.PeakTipDisplacement());
        }
        const float energy1 = model.TotalEnergy();
        const float drift = std::abs(energy1 - energy0) / energy0;
        check(model.IsFinite(), "finite state", "1 s undamped run");
        check(peak > 1e-5f && peak < 0.05f, "rod impulse response",
              std::format("peak {:.5f} m", peak));
        check(drift < 1e-3f, "undamped energy drift",
              std::format("relative drift {:.2e}", drift));
    }

    {
        Waterphone model;
        model.SetDamping(0.12f);
        model.Reset();
        model.Strike(0, kInitialContact, kStrikeVelocity);
        for (int step = 0; step < kStepsPerSecond / 10; ++step) {
            model.Step();
        }
        const float energy0 = model.TotalEnergy();
        for (int step = 0; step < kStepsPerSecond * 4; ++step) {
            model.Step();
        }
        const float energy1 = model.TotalEnergy();
        check(energy1 > 0.0f && energy1 < 0.8f * energy0,
              "damped energy decay",
              std::format("E {:.4f} -> {:.4f} mJ", energy0 * 1000.0f,
                          energy1 * 1000.0f));
        check(model.IsFinite(), "damped finite state", "4 s run");
    }

    {
        Waterphone model;
        model.SetDamping(0.12f);
        model.Reset();
        model.SetBow(1, 0.82f, 0.55f, kBowNormal);
        for (int step = 0; step < kStepsPerSecond; ++step) {
            model.Step();
        }
        const float energy = model.TotalEnergy();
        const float peak = model.PeakTipDisplacement();
        check(energy > 1e-8f && peak > 1e-6f, "bow excitation",
              std::format("E {:.4f} mJ peak {:.5f} m", energy * 1000.0f,
                          peak));
        check(model.IsFinite(), "bow finite state", "1 s bowed run");
    }

    std::cout << (failures == 0 ? "selftest ok\n" : "selftest failed\n");
    return failures == 0 ? 0 : 1;
}

}  // namespace qwqdsp_vibrate

int main(int argc, char** argv) {
    if (argc > 1 && std::string(argv[1]) == "--selftest") {
        return qwqdsp_vibrate::RunSelfTest();
    }

    using namespace qwqdsp_vibrate;

    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_MSAA_4X_HINT | FLAG_VSYNC_HINT);
    InitWindow(1280, 900, "Waterphone");
    SetWindowMinSize(960, 720);
    SetTargetFPS(60);

    Waterphone model;
    model.Strike(kInitialRod, kInitialContact, kStrikeVelocity);
    ProbeHistory history;
    int probe_rod = kInitialRod;
    float wave_scale = 0.01f;
    int sample_countdown = kStepsPerSample;
    float accumulator = 0.0f;
    float speed = kSpeedDefault;
    bool speed_dragging = false;
    bool pad_dragging = false;
    int press_rod = -1;
    float press_fraction = 0.0f;
    Vector2 press_position{};
    bool press_moved = false;
    float strike_flash = 0.0f;
    int strike_rod = -1;
    float strike_fraction = 0.0f;

    const auto select_rod = [&](int rod) {
        if (rod != probe_rod) {
            probe_rod = rod;
            history.Clear();
            wave_scale = 0.01f;
        }
    };

    while (!WindowShouldClose()) {
        const int width = GetScreenWidth();
        const int height = GetScreenHeight();
        const Layout layout = MakeLayout(width, height);
        const Vector2 mouse{static_cast<float>(GetMouseX()),
                            static_cast<float>(GetMouseY())};
        const Rectangle speed_hit = SpeedSliderHit(width);
        const Rectangle speed_track = SpeedSliderTrack(width);
        const Rectangle pad_box = PadBox(width);
        const Rectangle pad_inner = PadInner(pad_box);

        if (IsKeyPressed(KEY_R)) {
            model.Reset();
            model.Strike(kInitialRod, kInitialContact, kStrikeVelocity);
            probe_rod = kInitialRod;
            history.Clear();
            wave_scale = 0.01f;
            sample_countdown = kStepsPerSample;
            accumulator = 0.0f;
            press_rod = -1;
            press_moved = false;
            pad_dragging = false;
            model.ClearBow();
        }
        if (IsKeyPressed(KEY_E)) {
            std::cout << std::format("wav export: rendering {:.0f} s...\n",
                                     kExportSeconds)
                      << std::flush;
            Waterphone export_model;
            export_model.SetDamping(kExportDamping);
            export_model.Reset();
            export_model.Strike(probe_rod, kStrikeContact, kStrikeVelocity);
            const bool ok = ExportProbeWav(
                "vibrate_waterphone.wav", export_model,
                kExportStepsPerSecond1D, kExportSeconds,
                [](Waterphone& m, float) { m.Step(); },
                [&](const Waterphone& m) {
                    return m.TipDisplacement(probe_rod);
                });
            std::cout << std::format("wav export: {}\n",
                                     ok ? "vibrate_waterphone.wav" : "failed");
        }

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            if (CheckCollisionPointRec(mouse, speed_hit)) {
                speed_dragging = true;
            } else if (CheckCollisionPointRec(mouse, pad_box)) {
                pad_dragging = true;
                select_rod(PadRod(pad_inner, mouse));
            } else {
                const RodHit hit = HitTestRod(model, layout, mouse);
                if (hit.rod >= 0) {
                    select_rod(hit.rod);
                    press_rod = hit.rod;
                    press_fraction = hit.fraction;
                    press_position = mouse;
                    press_moved = false;
                    model.ClearBow();
                }
            }
        }

        if (speed_dragging) {
            const bool near_slider =
                mouse.y >= 20.0f && mouse.y <= 62.0f &&
                mouse.x >= speed_track.x - 20.0f &&
                mouse.x <= speed_track.x + speed_track.width + 20.0f;
            if (!near_slider || !IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                speed_dragging = false;
            } else {
                const float fraction = std::clamp(
                    (mouse.x - speed_track.x) / speed_track.width, 0.0f,
                    1.0f);
                speed = kSpeedMin + fraction * (kSpeedMax - kSpeedMin);
            }
        }

        if (pad_dragging) {
            if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                pad_dragging = false;
            } else if (CheckCollisionPointRec(mouse, pad_box)) {
                select_rod(PadRod(pad_inner, mouse));
            }
        }

        if (press_rod >= 0) {
            if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                if (!press_moved) {
                    model.Strike(press_rod, press_fraction, kStrikeVelocity);
                    history.Clear();
                    wave_scale = 0.01f;
                    strike_flash = 0.18f;
                    strike_rod = press_rod;
                    strike_fraction = press_fraction;
                }
                model.ClearBow();
                press_rod = -1;
                press_moved = false;
            } else {
                const float moved =
                    std::hypot(mouse.x - press_position.x,
                               mouse.y - press_position.y);
                if (!press_moved && moved > 4.0f) {
                    press_moved = true;
                    history.Clear();
                    wave_scale = 0.01f;
                }
                if (press_moved) {
                    const float bow_speed =
                        BowSpeed(layout, press_rod, GetMouseDelta(),
                                 GetFrameTime());
                    model.SetBow(press_rod, press_fraction, bow_speed,
                                 kBowNormal);
                }
            }
        } else {
            model.ClearBow();
        }

        strike_flash = std::max(0.0f, strike_flash - GetFrameTime());

        const float dt = std::min(GetFrameTime(), kMaxFrameTime) * speed;
        accumulator += dt;
        int steps = static_cast<int>(accumulator / kFixedDt);
        steps = std::min(steps, kMaxStepsPerFrame);
        accumulator -= static_cast<float>(steps) * kFixedDt;
        for (int i = 0; i < steps; ++i) {
            model.Step();
            if (--sample_countdown <= 0) {
                sample_countdown = kStepsPerSample;
                history.Push(model.TipDisplacement(probe_rod));
            }
        }

        BeginDrawing();
        ClearBackground(Color{5, 7, 11, 255});
        DrawSceneFloor(layout);
        DrawDisk(layout);
        for (int rod = 0; rod < kRodCount; ++rod) {
            if (rod != probe_rod) {
                DrawRod(model, layout, rod, false);
            }
        }
        DrawRod(model, layout, probe_rod, true);

        if (model.IsBowing()) {
            const int rod = model.BowingRod();
            const Vector2 point = RodDisplayPoint(
                model, layout, rod, model.BowingFraction());
            DrawCircleV(point, 3.5f, Color{255, 235, 190, 255});
            DrawCircleLinesV(point, 7.0f, Color{255, 200, 140, 190});
            DrawLineEx(Vector2{point.x - 10.0f, point.y - 4.0f},
                       Vector2{point.x + 10.0f, point.y + 4.0f}, 1.5f,
                       Color{255, 225, 175, 170});
        }
        if (strike_flash > 0.0f && strike_rod >= 0) {
            const Vector2 point =
                RodDisplayPoint(model, layout, strike_rod, strike_fraction);
            const float t = strike_flash / 0.18f;
            DrawCircleLinesV(point, 5.0f + 18.0f * (1.0f - t),
                             Color{255, 210, 160,
                                   static_cast<unsigned char>(180.0f * t)});
        }

        DrawPad(model, pad_box, probe_rod);
        DrawSceneHud(model);
        DrawSpeedSlider(width, speed);
        DrawScope(history, model, probe_rod, wave_scale);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
