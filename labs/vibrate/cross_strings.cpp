#include <raylib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

#include "wav_export.hpp"

// 两根长度相同、互相垂直的一维弦，均沿 z 方向振动。
// 两根弦的中心共用一个无质量连接点，因此中心位移与速度始终相同。
namespace qwqdsp_vibrate {

constexpr int kNodeCount = 321;
constexpr int kCenter = (kNodeCount - 1) / 2;
constexpr float kLength = 3.0f;
constexpr float kSpacing = kLength / static_cast<float>(kNodeCount - 1);
constexpr float kTension = 80.0f;
constexpr float kLinearDensity = 0.25f;
constexpr float kNodeMass = kLinearDensity * kSpacing;
constexpr float kDefaultDamping = 0.12f;

constexpr int kStepsPerSecond = 16000;
constexpr float kFixedDt = 1.0f / static_cast<float>(kStepsPerSecond);
constexpr float kMaxFrameTime = 0.05f;
constexpr int kMaxStepsPerFrame = 900;

// 短时升余弦力脉冲，总冲量 1.2 N·s；8 ms 对 3 m 长的弦近似 delta。
constexpr float kImpulse = 1.2f;
constexpr float kImpulseDuration = 0.008f;
constexpr float kImpulsePeak = kImpulse / (0.5f * kImpulseDuration);

constexpr float kHistorySeconds = 4.0f;
constexpr int kHistoryRate = 500;
constexpr int kHistorySize = 2048;
constexpr int kHistoryWindowSamples = static_cast<int>(kHistorySeconds * kHistoryRate);
constexpr int kStepsPerSample = kStepsPerSecond / kHistoryRate;

constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kTwoPi = 2.0f * kPi;

class CrossStrings {
public:
    CrossStrings()
        : x_(kNodeCount, 0.0f),
          y_(kNodeCount, 0.0f),
          vx_(kNodeCount, 0.0f),
          vy_(kNodeCount, 0.0f),
          ax_(kNodeCount, 0.0f),
          ay_(kNodeCount, 0.0f) {
        ResetImpulse();
    }

    void ResetImpulse() {
        std::fill(x_.begin(), x_.end(), 0.0f);
        std::fill(y_.begin(), y_.end(), 0.0f);
        std::fill(vx_.begin(), vx_.end(), 0.0f);
        std::fill(vy_.begin(), vy_.end(), 0.0f);
        std::fill(ax_.begin(), ax_.end(), 0.0f);
        std::fill(ay_.begin(), ay_.end(), 0.0f);
        time_ = 0.0f;
    }

    void Step(float dt) {
        float center_accel = 0.0f;
        ComputeAccelerations(time_, center_accel);

        for (int i = 1; i < kNodeCount - 1; ++i) {
            if (i == kCenter) {
                continue;
            }
            vx_[i] += 0.5f * ax_[i] * dt;
            x_[i] += vx_[i] * dt;
            vy_[i] += 0.5f * ay_[i] * dt;
            y_[i] += vy_[i] * dt;
        }
        const float center_velocity_half = vx_[kCenter] + 0.5f * center_accel * dt;
        const float center_position = x_[kCenter] + center_velocity_half * dt;
        vx_[kCenter] = center_velocity_half;
        vy_[kCenter] = center_velocity_half;
        x_[kCenter] = center_position;
        y_[kCenter] = center_position;

        ComputeAccelerations(time_ + dt, center_accel);
        for (int i = 1; i < kNodeCount - 1; ++i) {
            if (i == kCenter) {
                continue;
            }
            vx_[i] += 0.5f * ax_[i] * dt;
            vy_[i] += 0.5f * ay_[i] * dt;
        }
        const float center_velocity = vx_[kCenter] + 0.5f * center_accel * dt;
        vx_[kCenter] = center_velocity;
        vy_[kCenter] = center_velocity;

        time_ += dt;
    }

    float TotalEnergy() const {
        float kinetic = 0.0f;
        for (int i = 0; i < kNodeCount; ++i) {
            kinetic += 0.5f * kNodeMass *
                       (vx_[i] * vx_[i] + vy_[i] * vy_[i]);
        }

        // 与线性弦力对应的二次势能：U = T / (2 dx) * sum(dy^2)。
        float potential = 0.0f;
        for (int i = 0; i < kNodeCount - 1; ++i) {
            const float dx = x_[i + 1] - x_[i];
            const float dy = y_[i + 1] - y_[i];
            potential += 0.5f * kTension / kSpacing * (dx * dx + dy * dy);
        }
        return kinetic + potential;
    }

    float PeakDisplacement() const {
        float peak = 0.0f;
        for (int i = 0; i < kNodeCount; ++i) {
            peak = std::max({peak, std::abs(x_[i]), std::abs(y_[i])});
        }
        return peak;
    }

    float MaxAsymmetry() const {
        float asymmetry = 0.0f;
        for (int i = 0; i < kNodeCount; ++i) {
            asymmetry = std::max(asymmetry, std::abs(x_[i] - y_[i]));
        }
        return asymmetry;
    }

    bool IsFinite() const {
        for (int i = 0; i < kNodeCount; ++i) {
            if (!std::isfinite(x_[i]) || !std::isfinite(y_[i]) ||
                !std::isfinite(vx_[i]) || !std::isfinite(vy_[i])) {
                return false;
            }
        }
        return true;
    }

    void SetDamping(float damping) { damping_ = damping; }
    float Time() const { return time_; }
    float CenterDisplacement() const { return x_[kCenter]; }
    float CenterVelocity() const { return vx_[kCenter]; }
    const std::vector<float>& X() const { return x_; }
    const std::vector<float>& Y() const { return y_; }

private:
    void ComputeAccelerations(float time, float& center_accel) {
        const float wave_speed_sq = kTension / kLinearDensity;
        const float inv_spacing_sq = 1.0f / (kSpacing * kSpacing);
        for (int i = 1; i < kNodeCount - 1; ++i) {
            if (i == kCenter) {
                continue;
            }
            ax_[i] = wave_speed_sq * inv_spacing_sq *
                         (x_[i + 1] - 2.0f * x_[i] + x_[i - 1]) -
                     damping_ * vx_[i];
            ay_[i] = wave_speed_sq * inv_spacing_sq *
                         (y_[i + 1] - 2.0f * y_[i] + y_[i - 1]) -
                     damping_ * vy_[i];
        }

        const float curvature_x =
            x_[kCenter + 1] - 2.0f * x_[kCenter] + x_[kCenter - 1];
        const float curvature_y =
            y_[kCenter + 1] - 2.0f * y_[kCenter] + y_[kCenter - 1];
        center_accel = 0.5f * wave_speed_sq * inv_spacing_sq *
                           (curvature_x + curvature_y) -
                       damping_ * vx_[kCenter];

        if (time < kImpulseDuration) {
            const float envelope =
                0.5f * (1.0f - std::cos(kTwoPi * time / kImpulseDuration));
            center_accel += (kImpulsePeak * envelope) / (2.0f * kNodeMass);
        }
    }

    std::vector<float> x_;
    std::vector<float> y_;
    std::vector<float> vx_;
    std::vector<float> vy_;
    std::vector<float> ax_;
    std::vector<float> ay_;
    float time_ = 0.0f;
    float damping_ = kDefaultDamping;
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
    float z_scale = 420.0f;
};

Layout MakeLayout(int width, int height) {
    Layout layout;
    const float w = static_cast<float>(width);
    const float h = static_cast<float>(height);
    const float scene_height = h * 0.63f;
    layout.scene = Rectangle{0.0f, 0.0f, w, scene_height};
    layout.center = Vector2{w * 0.5f, scene_height * 0.55f};
    layout.axis_scale = std::min(w * 0.26f, scene_height * 0.60f);
    layout.depth_x = layout.axis_scale * 0.55f;
    layout.depth_y = -layout.axis_scale * 0.30f;
    layout.z_scale = scene_height * 0.75f;

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

int ScreenWidth() {
    return GetScreenWidth();
}

int ScreenHeight() {
    return GetScreenHeight();
}

Color ProbeColor(int axis, int node) {
    if (node == kCenter) {
        return Color{255, 240, 205, 235};
    }
    return axis == 0 ? Color{140, 205, 255, 235} : Color{255, 195, 120, 235};
}

void DrawPlane(const Layout& layout) {
    const Color grid{24, 32, 45, 110};
    const Color border{38, 51, 70, 180};
    for (int i = 0; i <= 4; ++i) {
        const float w = -kLength * 0.5f + kLength * static_cast<float>(i) / 4.0f;
        const Color color = (i == 0 || i == 4) ? border : grid;
        DrawLineV(Project(layout, -kLength * 0.5f, w, 0.0f),
                  Project(layout, kLength * 0.5f, w, 0.0f), color);
        DrawLineV(Project(layout, w, -kLength * 0.5f, 0.0f),
                  Project(layout, w, kLength * 0.5f, 0.0f), color);
    }
    DrawLineV(Project(layout, -kLength * 0.5f, 0.0f, 0.0f),
              Project(layout, kLength * 0.5f, 0.0f, 0.0f),
              Color{52, 70, 96, 170});
    DrawLineV(Project(layout, 0.0f, -kLength * 0.5f, 0.0f),
              Project(layout, 0.0f, kLength * 0.5f, 0.0f),
              Color{52, 70, 96, 170});
}

void DrawString(const std::vector<float>& z, bool along_x, const Layout& layout,
                Color color) {
    for (int i = 0; i < kNodeCount - 1; ++i) {
        const float w0 = static_cast<float>(i - kCenter) * kSpacing;
        const float w1 = static_cast<float>(i + 1 - kCenter) * kSpacing;
        const Vector2 p0 = along_x ? Project(layout, w0, 0.0f, z[i])
                                   : Project(layout, 0.0f, w0, z[i]);
        const Vector2 p1 = along_x ? Project(layout, w1, 0.0f, z[i + 1])
                                   : Project(layout, 0.0f, w1, z[i + 1]);
        DrawLineEx(p0, p1, 7.0f, Color{color.r, color.g, color.b, 30});
    }
    for (int i = 0; i < kNodeCount - 1; ++i) {
        const float w0 = static_cast<float>(i - kCenter) * kSpacing;
        const float w1 = static_cast<float>(i + 1 - kCenter) * kSpacing;
        const Vector2 p0 = along_x ? Project(layout, w0, 0.0f, z[i])
                                   : Project(layout, 0.0f, w0, z[i]);
        const Vector2 p1 = along_x ? Project(layout, w1, 0.0f, z[i + 1])
                                   : Project(layout, 0.0f, w1, z[i + 1]);
        DrawLineEx(p0, p1, 2.0f, color);
    }
}

void DrawEndpoints(const Layout& layout) {
    const std::array<Vector2, 4> points{
        Project(layout, -kLength * 0.5f, 0.0f, 0.0f),
        Project(layout, kLength * 0.5f, 0.0f, 0.0f),
        Project(layout, 0.0f, -kLength * 0.5f, 0.0f),
        Project(layout, 0.0f, kLength * 0.5f, 0.0f),
    };
    for (const Vector2 point : points) {
        DrawCircleV(point, 3.5f, Color{150, 172, 200, 255});
    }
}

void DrawJoint(const CrossStrings& model, const Layout& layout) {
    const Vector2 base = Project(layout, 0.0f, 0.0f, 0.0f);
    const Vector2 joint = Project(layout, 0.0f, 0.0f, model.CenterDisplacement());
    const float top = std::min(base.y, joint.y);
    const float bottom = std::max(base.y, joint.y);
    for (float y = top + 4.0f; y + 4.0f < bottom; y += 10.0f) {
        DrawLineV(Vector2{joint.x, y}, Vector2{joint.x, y + 4.0f},
                  Color{255, 245, 220, 90});
    }
    DrawCircleV(joint, 4.5f, Color{255, 248, 225, 255});
}

void DrawProbeGuide(const CrossStrings& model, const Layout& layout, int axis, int node,
                    Color color) {
    if (node == kCenter) {
        const Vector2 joint = Project(layout, 0.0f, 0.0f, model.CenterDisplacement());
        DrawCircleLinesV(joint, 8.0f, color);
        return;
    }
    const float w = static_cast<float>(node - kCenter) * kSpacing;
    const float value = axis == 0 ? model.X()[node] : model.Y()[node];
    const Vector2 base = axis == 0 ? Project(layout, w, 0.0f, 0.0f)
                                   : Project(layout, 0.0f, w, 0.0f);
    const Vector2 point = axis == 0 ? Project(layout, w, 0.0f, value)
                                    : Project(layout, 0.0f, w, value);
    const float top = std::min(base.y, point.y);
    const float bottom = std::max(base.y, point.y);
    for (float y = top + 4.0f; y + 4.0f < bottom; y += 10.0f) {
        DrawLineV(Vector2{point.x, y}, Vector2{point.x, y + 4.0f},
                  Color{color.r, color.g, color.b, 90});
    }
    DrawCircleV(point, 3.0f, color);
    DrawCircleLinesV(point, 7.0f, color);
}

Rectangle ProbeBox(int width) {
    return Rectangle{static_cast<float>(width) - 168.0f, 24.0f, 144.0f, 144.0f};
}

Rectangle ProbeBoxInner(const Rectangle& box) {
    return Rectangle{box.x + 10.0f, box.y + 10.0f, box.width - 20.0f,
                     box.height - 20.0f};
}

void DrawProbeBox(const Rectangle& box, int axis, int node) {
    const Rectangle inner = ProbeBoxInner(box);
    const float mid_x = inner.x + inner.width * 0.5f;
    const float mid_y = inner.y + inner.height * 0.5f;

    DrawRectangleRec(box, Color{6, 8, 13, 225});
    DrawRectangleLinesEx(box, 1.0f, Color{44, 60, 82, 230});
    DrawText("xy", static_cast<int>(box.x) + 8, static_cast<int>(box.y) + 6, 12,
             Color{100, 120, 145, 220});
    DrawLineV(Vector2{inner.x, mid_y},
              Vector2{inner.x + inner.width, mid_y}, Color{70, 140, 200, 150});
    DrawLineV(Vector2{mid_x, inner.y},
              Vector2{mid_x, inner.y + inner.height}, Color{220, 150, 80, 150});

    const float wx = node == kCenter || axis != 0
                         ? 0.0f
                         : static_cast<float>(node - kCenter) * kSpacing;
    const float wy = node == kCenter || axis != 1
                         ? 0.0f
                         : static_cast<float>(node - kCenter) * kSpacing;
    const float px = inner.x + (wx + kLength * 0.5f) / kLength * inner.width;
    const float py = inner.y + (kLength * 0.5f - wy) / kLength * inner.height;
    const Color color = ProbeColor(axis, node);
    DrawCircleV(Vector2{px, py}, 4.0f, color);
    DrawCircleLinesV(Vector2{px, py}, 7.0f, Color{255, 255, 255, 170});
}

void DrawSceneHud(const CrossStrings& model) {
    const std::string line = std::format(
        "t {:.2f} s      R: reset impulse      drag xy box: probe      E: export wav",
        model.Time());
    DrawText(line.c_str(), 24, 20, 15, Color{120, 140, 165, 230});
}

std::string ProbeTitle(int axis, int node) {
    if (node == kCenter) {
        return std::format("PROBE  JOINT  node {}", node);
    }
    const float offset = static_cast<float>(node - kCenter) * kSpacing;
    if (axis == 0) {
        return std::format("PROBE  X  node {}   x {:+.2f} m", node, offset);
    }
    return std::format("PROBE  Y  node {}   y {:+.2f} m", node, offset);
}

void DrawScope(const ProbeHistory& history, const CrossStrings& model, int axis,
               int node, float& scale) {
    const Layout layout = MakeLayout(ScreenWidth(), ScreenHeight());
    const Rectangle panel = layout.scope;
    const Rectangle plot{panel.x + 18.0f, panel.y + 44.0f, panel.width - 36.0f,
                         panel.height - 76.0f};
    const float mid = plot.y + plot.height * 0.5f;
    const int count = std::min(history.Size(), kHistoryWindowSamples);
    const float value = axis == 0 ? model.X()[node] : model.Y()[node];

    float peak = 0.0f;
    for (int age = 0; age < count; ++age) {
        peak = std::max(peak, std::abs(history.At(age)));
    }
    if (count > 0) {
        const float target = std::max(peak * 1.15f, 0.004f);
        if (target > scale) {
            scale = target;
        } else {
            scale += (target - scale) * 0.04f;
        }
    }

    DrawRectangleRec(panel, Color{7, 9, 14, 255});
    DrawRectangleLinesEx(panel, 1.0f, Color{42, 58, 80, 220});

    const std::string title = ProbeTitle(axis, node);
    DrawText(title.c_str(), static_cast<int>(panel.x) + 16,
             static_cast<int>(panel.y) + 12, 15, Color{200, 215, 235, 255});

    const std::string live = std::format("z {:+.3f} m      scale {:.3f} m",
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
    for (int i = 1; i <= 3; ++i) {
        const float x = plot.x + plot.width * static_cast<float>(i) / 4.0f;
        DrawLineV(Vector2{x, plot.y + 2.0f},
                  Vector2{x, plot.y + plot.height - 2.0f},
                  Color{36, 52, 74, 130});
    }

    const int label_y = static_cast<int>(plot.y + plot.height + 6.0f);
    for (int i = 1; i <= 3; ++i) {
        const std::string label = std::format("-{}s", 4 - i);
        const int label_width = MeasureText(label.c_str(), 13);
        const float x = plot.x + plot.width * static_cast<float>(i) / 4.0f;
        DrawText(label.c_str(), static_cast<int>(x) - label_width / 2, label_y, 13,
                 Color{100, 120, 145, 220});
    }
    DrawText("now",
             static_cast<int>(plot.x + plot.width - MeasureText("now", 13)),
             label_y, 13, Color{100, 120, 145, 220});

    const Color trace_color = ProbeColor(axis, node);

    if (count >= 2) {
        Vector2 previous{};
        bool has_previous = false;
        for (int age = count - 1; age >= 0; --age) {
            const float x = plot.x + plot.width *
                                        (1.0f - static_cast<float>(age) /
                                                    static_cast<float>(kHistoryWindowSamples - 1));
            const float y = std::clamp(mid - history.At(age) / scale *
                                                  (plot.height * 0.5f - 4.0f),
                                       plot.y + 1.0f, plot.y + plot.height - 1.0f);
            const Vector2 current{x, y};
            if (has_previous) {
                DrawLineEx(previous, current, 1.7f, trace_color);
            }
            previous = current;
            has_previous = true;
        }
    } else {
        const char* waiting = "waiting for samples...";
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
        std::cout << std::format("[{}] {:<28} {}\n", ok ? "PASS" : "FAIL", name,
                                 detail);
        if (!ok) {
            ++failures;
        }
    };

    {
        CrossStrings model;
        model.SetDamping(0.0f);
        model.ResetImpulse();
        for (int i = 0; i < kStepsPerSecond / 50; ++i) {
            model.Step(kFixedDt);
        }
        const float energy0 = model.TotalEnergy();
        float peak = 0.0f;
        for (int i = 0; i < kStepsPerSecond; ++i) {
            model.Step(kFixedDt);
            peak = std::max(peak, model.PeakDisplacement());
        }
        const float energy1 = model.TotalEnergy();
        const float drift = std::abs(energy1 - energy0) / energy0;
        check(model.IsFinite(), "finite state", "1 s undamped run");
        check(peak > 0.02f && peak < 0.4f, "impulse response",
              std::format("peak {:.4f} m", peak));
        check(model.MaxAsymmetry() < 1e-6f, "shared center symmetry",
              std::format("max |x-y| {:.2e} m", model.MaxAsymmetry()));
        check(drift < 0.02f, "undamped energy drift",
              std::format("relative drift {:.4f}", drift));
    }

    {
        CrossStrings model;
        model.SetDamping(0.15f);
        model.ResetImpulse();
        for (int i = 0; i < kStepsPerSecond / 50; ++i) {
            model.Step(kFixedDt);
        }
        const float energy0 = model.TotalEnergy();
        for (int i = 0; i < kStepsPerSecond * 4; ++i) {
            model.Step(kFixedDt);
        }
        const float energy1 = model.TotalEnergy();
        check(energy1 > 0.0f && energy1 < 0.8f * energy0, "damped energy decay",
              std::format("E {:.2f} -> {:.2f} mJ", energy0 * 1000.0f,
                          energy1 * 1000.0f));
        check(model.IsFinite(), "damped finite state", "4 s run");
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
    InitWindow(1280, 900, "Two Coupled Strings");
    SetWindowMinSize(960, 720);
    SetTargetFPS(60);

    CrossStrings model;
    ProbeHistory history;
    int probe_axis = 0;
    int probe_node = kCenter;
    float wave_scale = 0.05f;
    int sample_countdown = kStepsPerSample;
    float accumulator = 0.0f;
    bool probe_dragging = false;

    while (!WindowShouldClose()) {
        const int width = ScreenWidth();
        const int height = ScreenHeight();
        const Layout layout = MakeLayout(width, height);
        const Vector2 mouse{static_cast<float>(GetMouseX()),
                            static_cast<float>(GetMouseY())};

        if (IsKeyPressed(KEY_R)) {
            model.ResetImpulse();
            history.Clear();
            wave_scale = 0.05f;
            sample_countdown = kStepsPerSample;
            accumulator = 0.0f;
        }
        if (IsKeyPressed(KEY_E)) {
            std::cout << std::format("wav export: rendering {:.0f} s...\n",
                                     kExportSeconds)
                      << std::flush;
            CrossStrings export_model;
            export_model.SetDamping(kExportDamping);
            export_model.ResetImpulse();
            const bool ok = ExportProbeWav(
                "vibrate_cross.wav", export_model, kExportStepsPerSecond1D,
                kExportSeconds,
                [](CrossStrings& m, float dt) { m.Step(dt); },
                [&](const CrossStrings& m) {
                    return probe_axis == 0 ? m.X()[probe_node]
                                           : m.Y()[probe_node];
                });
            std::cout << std::format("wav export: {}\n",
                                     ok ? "vibrate_cross.wav" : "failed");
        }

        const Rectangle box = ProbeBox(width);
        const Rectangle inner = ProbeBoxInner(box);
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) &&
            CheckCollisionPointRec(mouse, box)) {
            probe_dragging = true;
        }
        if (probe_dragging && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            const float mx = std::clamp(mouse.x, inner.x, inner.x + inner.width);
            const float my = std::clamp(mouse.y, inner.y, inner.y + inner.height);
            const float wx =
                (mx - inner.x) / inner.width * kLength - kLength * 0.5f;
            const float wy =
                kLength * 0.5f - (my - inner.y) / inner.height * kLength;
            const bool x_closer = std::abs(wy) <= std::abs(wx);
            const int axis = x_closer ? 0 : 1;
            const float coordinate = x_closer ? wx : wy;
            const int node = std::clamp(
                static_cast<int>(std::lround(
                    (coordinate + kLength * 0.5f) / kLength *
                    static_cast<float>(kNodeCount - 1))),
                1, kNodeCount - 2);
            if (axis != probe_axis || node != probe_node) {
                probe_axis = axis;
                probe_node = node;
                history.Clear();
                wave_scale = 0.05f;
                sample_countdown = kStepsPerSample;
            }
        }
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) {
            probe_dragging = false;
        }

        accumulator += std::min(GetFrameTime(), kMaxFrameTime);
        int steps = 0;
        while (accumulator >= kFixedDt && steps < kMaxStepsPerFrame) {
            model.Step(kFixedDt);
            if (--sample_countdown <= 0) {
                const float value = probe_axis == 0 ? model.X()[probe_node]
                                                    : model.Y()[probe_node];
                history.Push(value);
                sample_countdown = kStepsPerSample;
            }
            accumulator -= kFixedDt;
            ++steps;
        }

        BeginDrawing();
        ClearBackground(Color{8, 10, 16, 255});

        BeginScissorMode(0, 0, width, static_cast<int>(layout.scene.height));
        DrawPlane(layout);
        DrawString(model.Y(), false, layout, Color{255, 190, 110, 255});
        DrawString(model.X(), true, layout, Color{120, 200, 255, 255});
        DrawEndpoints(layout);
        DrawJoint(model, layout);
        DrawProbeGuide(model, layout, probe_axis, probe_node,
                       ProbeColor(probe_axis, probe_node));
        EndScissorMode();

        DrawProbeBox(ProbeBox(width), probe_axis, probe_node);
        DrawSceneHud(model);
        DrawScope(history, model, probe_axis, probe_node, wave_scale);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
