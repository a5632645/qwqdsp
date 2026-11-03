#include <raylib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <string>
#include <vector>

#include "wav_export.hpp"

// 4 根等长弦（4 条直径）在中心共用一个连接点；
// 8 个外端连接在外圈圆环上，圆环自身也是可振动的 1D 环形弦。
namespace qwqdsp_vibrate {

constexpr int kStringCount = 4;
constexpr int kStringNodes = 241;
constexpr int kStringCenter = (kStringNodes - 1) / 2;
constexpr float kOuterRadius = 1.5f;
constexpr float kStringLength = 2.0f * kOuterRadius;
constexpr float kStringSpacing = kStringLength / static_cast<float>(kStringNodes - 1);
// 四根弦张力略有偏差，产生轻微 detune；整体张力水平决定音高。
constexpr std::array<float, kStringCount> kStringTensions = {
    79.0f, 80.6f, 78.6f, 81.4f};
constexpr float kStringDensity = 0.25f;
constexpr float kStringNodeMass = kStringDensity * kStringSpacing;
constexpr std::array<float, kStringCount> kStringWaveSq = {
    kStringTensions[0] / kStringDensity,
    kStringTensions[1] / kStringDensity,
    kStringTensions[2] / kStringDensity,
    kStringTensions[3] / kStringDensity,
};

constexpr int kRingNodes = 256;
constexpr int kRingTerminalStep = kRingNodes / 8;
// 4 个固定点位于相邻弦端之间（22.5° + 90°），抑制环的刚体上漂，
// 同时保留环的 n=4 模态与 8 个弦端的耦合。
constexpr int kRingFixedStride = kRingNodes / 4;
constexpr int kRingFixedOffset = kRingTerminalStep / 2;
bool IsRingFixed(int k) { return (k % kRingFixedStride) == kRingFixedOffset; }
constexpr float kRingCircumference = 6.28318530718f * kOuterRadius;
constexpr float kRingSpacing = kRingCircumference / static_cast<float>(kRingNodes);
constexpr float kRingTotalMass = 0.05f;
constexpr float kRingNodeMass = kRingTotalMass / static_cast<float>(kRingNodes);
constexpr float kRingTension = 125.0f;
// constexpr float kRingTension = 5.0f;
constexpr float kRingDensity = kRingTotalMass / kRingCircumference;
constexpr float kRingWaveSq = kRingTension / kRingDensity;
constexpr float kDefaultDamping = 0.12f;

constexpr int kStepsPerSecond = 16000;
constexpr float kFixedDt = 1.0f / static_cast<float>(kStepsPerSecond);
constexpr float kMaxFrameTime = 0.05f;
constexpr int kMaxStepsPerFrame = 900;
constexpr float kSpeedMin = 0.10f;
constexpr float kSpeedMax = 1.00f;
constexpr float kSpeedDefault = 0.25f;
constexpr float kImpulseVelocity = 20.0f;

constexpr float kHistorySeconds = 4.0f;
constexpr int kHistoryRate = 500;
constexpr int kHistorySize = 2048;
constexpr int kHistoryWindowSamples = static_cast<int>(kHistorySeconds * kHistoryRate);
constexpr int kStepsPerSample = kStepsPerSecond / kHistoryRate;

constexpr float kTwoPi = 6.28318530718f;

class RingDiameters {
public:
    RingDiameters() {
        for (int s = 0; s < kStringCount; ++s) {
            u_[s].assign(kStringNodes, 0.0f);
            v_[s].assign(kStringNodes, 0.0f);
            a_[s].assign(kStringNodes, 0.0f);
        }
        ring_u_.assign(kRingNodes, 0.0f);
        ring_v_.assign(kRingNodes, 0.0f);
        ring_a_.assign(kRingNodes, 0.0f);
        ResetImpulse();
    }

    void ResetImpulse() {
        for (int s = 0; s < kStringCount; ++s) {
            std::fill(u_[s].begin(), u_[s].end(), 0.0f);
            std::fill(v_[s].begin(), v_[s].end(), 0.0f);
            std::fill(a_[s].begin(), a_[s].end(), 0.0f);
        }
        std::fill(ring_u_.begin(), ring_u_.end(), 0.0f);
        std::fill(ring_v_.begin(), ring_v_.end(), 0.0f);
        std::fill(ring_a_.begin(), ring_a_.end(), 0.0f);
        center_z_ = 0.0f;
        center_v_ = kImpulseVelocity;
        center_a_ = 0.0f;
        time_ = 0.0f;
        SetSharedBoundaries();
    }

    void Step(float dt) {
        ComputeAccelerations();

        for (int s = 0; s < kStringCount; ++s) {
            for (int i = 1; i < kStringNodes - 1; ++i) {
                if (i == kStringCenter) {
                    continue;
                }
                v_[s][i] += 0.5f * a_[s][i] * dt;
                u_[s][i] += v_[s][i] * dt;
            }
        }
        center_v_ += 0.5f * center_a_ * dt;
        center_z_ += center_v_ * dt;
        for (int k = 0; k < kRingNodes; ++k) {
            if (IsRingFixed(k)) {
                ring_u_[k] = 0.0f;
                ring_v_[k] = 0.0f;
                continue;
            }
            ring_v_[k] += 0.5f * ring_a_[k] * dt;
            ring_u_[k] += ring_v_[k] * dt;
        }
        SetSharedBoundaries();

        ComputeAccelerations();
        for (int s = 0; s < kStringCount; ++s) {
            for (int i = 1; i < kStringNodes - 1; ++i) {
                if (i == kStringCenter) {
                    continue;
                }
                v_[s][i] += 0.5f * a_[s][i] * dt;
            }
        }
        center_v_ += 0.5f * center_a_ * dt;
        for (int k = 0; k < kRingNodes; ++k) {
            if (IsRingFixed(k)) {
                ring_u_[k] = 0.0f;
                ring_v_[k] = 0.0f;
                continue;
            }
            ring_v_[k] += 0.5f * ring_a_[k] * dt;
        }
        SetSharedBoundaries();
        time_ += dt;
    }

    float TotalEnergy() const {
        float kinetic = 0.0f;
        for (int s = 0; s < kStringCount; ++s) {
            for (int i = 1; i < kStringNodes - 1; ++i) {
                kinetic += 0.5f * kStringNodeMass * v_[s][i] * v_[s][i];
            }
        }
        for (int k = 0; k < kRingNodes; ++k) {
            kinetic += 0.5f * kRingNodeMass * ring_v_[k] * ring_v_[k];
        }

        float potential = 0.0f;
        for (int s = 0; s < kStringCount; ++s) {
            for (int i = 0; i < kStringNodes - 1; ++i) {
                const float difference = u_[s][i + 1] - u_[s][i];
                potential += 0.5f * kStringTensions[s] / kStringSpacing * difference *
                             difference;
            }
        }
        for (int k = 0; k < kRingNodes; ++k) {
            const float next = ring_u_[(k + 1) % kRingNodes];
            const float difference = next - ring_u_[k];
            potential += 0.5f * kRingTension / kRingSpacing * difference * difference;
        }
        return kinetic + potential;
    }

    float PeakDisplacement() const {
        float peak = std::abs(center_z_);
        for (int s = 0; s < kStringCount; ++s) {
            for (float value : u_[s]) {
                peak = std::max(peak, std::abs(value));
            }
        }
        for (float value : ring_u_) {
            peak = std::max(peak, std::abs(value));
        }
        return peak;
    }

    float MaxStringSpread() const {
        float spread = 0.0f;
        for (int i = 0; i < kStringNodes; ++i) {
            float low = u_[0][i];
            float high = u_[0][i];
            for (int s = 1; s < kStringCount; ++s) {
                low = std::min(low, u_[s][i]);
                high = std::max(high, u_[s][i]);
            }
            spread = std::max(spread, high - low);
        }
        return spread;
    }

    float MaxRingFixedDrift() const {
        float drift = 0.0f;
        for (int fixed = 0; fixed < 4; ++fixed) {
            const int k = fixed * kRingFixedStride + kRingFixedOffset;
            drift = std::max({drift, std::abs(ring_u_[k]), std::abs(ring_v_[k])});
        }
        return drift;
    }

    float MaxRingDisplacement() const {
        float peak = 0.0f;
        for (float value : ring_u_) {
            peak = std::max(peak, std::abs(value));
        }
        return peak;
    }

    bool IsFinite() const {
        if (!std::isfinite(center_z_) || !std::isfinite(center_v_)) {
            return false;
        }
        for (int s = 0; s < kStringCount; ++s) {
            for (int i = 1; i < kStringNodes - 1; ++i) {
                if (!std::isfinite(u_[s][i]) || !std::isfinite(v_[s][i])) {
                    return false;
                }
            }
        }
        for (int k = 0; k < kRingNodes; ++k) {
            if (!std::isfinite(ring_u_[k]) || !std::isfinite(ring_v_[k])) {
                return false;
            }
        }
        return true;
    }

    void SetDamping(float damping) { damping_ = damping; }
    float Time() const { return time_; }
    float CenterZ() const { return center_z_; }
    const std::vector<float>& String(int s) const { return u_[s]; }
    const std::vector<float>& Ring() const { return ring_u_; }

private:
    void SetSharedBoundaries() {
        for (int s = 0; s < kStringCount; ++s) {
            const int near = s * kRingTerminalStep;
            const int far = (near + kRingNodes / 2) % kRingNodes;
            u_[s][0] = ring_u_[near];
            u_[s][kStringNodes - 1] = ring_u_[far];
            u_[s][kStringCenter] = center_z_;
            v_[s][kStringCenter] = center_v_;
        }
    }

    void ComputeAccelerations() {
        const float inv_string_spacing_sq = 1.0f / (kStringSpacing * kStringSpacing);
        for (int s = 0; s < kStringCount; ++s) {
            for (int i = 1; i < kStringNodes - 1; ++i) {
                if (i == kStringCenter) {
                    continue;
                }
                a_[s][i] = kStringWaveSq[s] * inv_string_spacing_sq *
                               (u_[s][i + 1] - 2.0f * u_[s][i] + u_[s][i - 1]) -
                           damping_ * v_[s][i];
            }
        }

        float curvature_sum = 0.0f;
        for (int s = 0; s < kStringCount; ++s) {
            curvature_sum +=
                kStringWaveSq[s] * (u_[s][kStringCenter + 1] -
                                    2.0f * center_z_ +
                                    u_[s][kStringCenter - 1]);
        }
        center_a_ = inv_string_spacing_sq * 0.25f * curvature_sum -
                    damping_ * center_v_;

        const float inv_ring_spacing_sq = 1.0f / (kRingSpacing * kRingSpacing);
        for (int k = 0; k < kRingNodes; ++k) {
            if (IsRingFixed(k)) {
                ring_a_[k] = 0.0f;
                continue;
            }
            const float previous = ring_u_[(k - 1 + kRingNodes) % kRingNodes];
            const float next = ring_u_[(k + 1) % kRingNodes];
            ring_a_[k] = kRingWaveSq * inv_ring_spacing_sq *
                             (previous - 2.0f * ring_u_[k] + next) -
                         damping_ * ring_v_[k];
        }

        for (int s = 0; s < kStringCount; ++s) {
            const int near = s * kRingTerminalStep;
            const int far = (near + kRingNodes / 2) % kRingNodes;
            if (!IsRingFixed(near)) {
                ring_a_[near] +=
                    kStringTensions[s] / kStringSpacing *
                    (u_[s][1] - ring_u_[near]) / kRingNodeMass;
            }
            if (!IsRingFixed(far)) {
                ring_a_[far] +=
                    kStringTensions[s] / kStringSpacing *
                    (u_[s][kStringNodes - 2] - ring_u_[far]) / kRingNodeMass;
            }
        }
    }

    std::array<std::vector<float>, kStringCount> u_;
    std::array<std::vector<float>, kStringCount> v_;
    std::array<std::vector<float>, kStringCount> a_;
    std::vector<float> ring_u_;
    std::vector<float> ring_v_;
    std::vector<float> ring_a_;
    float center_z_ = 0.0f;
    float center_v_ = 0.0f;
    float center_a_ = 0.0f;
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
    float z_scale = 650.0f;
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
    layout.z_scale = scene_height * 3.0f;
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

Color StringColor(int string_index) {
    return (string_index == 0 || string_index == 2)
               ? Color{140, 205, 255, 255}
               : Color{255, 195, 120, 255};
}

Color ProbeColor(int string_index) {
    if (string_index < 0) {
        return Color{255, 240, 205, 235};
    }
    return StringColor(string_index);
}

void DrawPlane(const Layout& layout) {
    const Color grid{24, 32, 45, 110};
    const Color border{38, 51, 70, 180};
    for (int i = 0; i <= 4; ++i) {
        const float w = -kOuterRadius + 2.0f * kOuterRadius *
                                            static_cast<float>(i) / 4.0f;
        const Color color = (i == 0 || i == 4) ? border : grid;
        DrawLineV(Project(layout, -kOuterRadius, w, 0.0f),
                  Project(layout, kOuterRadius, w, 0.0f), color);
        DrawLineV(Project(layout, w, -kOuterRadius, 0.0f),
                  Project(layout, w, kOuterRadius, 0.0f), color);
    }
    constexpr int kCircleSegments = 96;
    for (int i = 0; i < kCircleSegments; ++i) {
        const float a0 = kTwoPi * static_cast<float>(i) / kCircleSegments;
        const float a1 = kTwoPi * static_cast<float>(i + 1) / kCircleSegments;
        DrawLineV(Project(layout, kOuterRadius * std::cos(a0),
                          kOuterRadius * std::sin(a0), 0.0f),
                  Project(layout, kOuterRadius * std::cos(a1),
                          kOuterRadius * std::sin(a1), 0.0f),
                  Color{52, 70, 96, 150});
    }
}

void DrawRing(const RingDiameters& model, const Layout& layout) {
    const std::vector<float>& values = model.Ring();
    for (int pass = 0; pass < 2; ++pass) {
        const float thickness = pass == 0 ? 6.0f : 1.8f;
        const Color color = pass == 0 ? Color{255, 220, 160, 28}
                                      : Color{255, 235, 195, 220};
        for (int k = 0; k < kRingNodes; ++k) {
            const int next = (k + 1) % kRingNodes;
            const float a0 = kTwoPi * static_cast<float>(k) / kRingNodes;
            const float a1 = kTwoPi * static_cast<float>(next) / kRingNodes;
            const Vector2 p0 = Project(layout, kOuterRadius * std::cos(a0),
                                       kOuterRadius * std::sin(a0), values[k]);
            const Vector2 p1 = Project(layout, kOuterRadius * std::cos(a1),
                                       kOuterRadius * std::sin(a1), values[next]);
            DrawLineEx(p0, p1, thickness, color);
        }
    }
    for (int terminal = 0; terminal < 8; ++terminal) {
        const int k = terminal * kRingTerminalStep;
        const float angle = kTwoPi * static_cast<float>(k) / kRingNodes;
        const Vector2 p = Project(layout, kOuterRadius * std::cos(angle),
                                  kOuterRadius * std::sin(angle), values[k]);
        DrawCircleV(p, 3.0f, Color{255, 245, 220, 220});
    }
    for (int fixed = 0; fixed < 4; ++fixed) {
        const int k = fixed * kRingFixedStride + kRingFixedOffset;
        const float angle = kTwoPi * static_cast<float>(k) / kRingNodes;
        const Vector2 p = Project(layout, kOuterRadius * std::cos(angle),
                                  kOuterRadius * std::sin(angle), 0.0f);
        DrawCircleV(p, 4.0f, Color{255, 150, 110, 255});
        DrawCircleLinesV(p, 7.0f, Color{255, 190, 150, 190});
        DrawLineV(Vector2{p.x - 6.0f, p.y}, Vector2{p.x + 6.0f, p.y},
                  Color{255, 190, 150, 150});
        DrawLineV(Vector2{p.x, p.y - 6.0f}, Vector2{p.x, p.y + 6.0f},
                  Color{255, 190, 150, 150});
    }
}

void DrawStrings(const RingDiameters& model, const Layout& layout) {
    for (int s = 0; s < kStringCount; ++s) {
        const float angle = kTwoPi * static_cast<float>(s) / 8.0f;
        const float ux = std::cos(angle);
        const float uy = std::sin(angle);
        const std::vector<float>& values = model.String(s);
        const Color color = StringColor(s);
        for (int pass = 0; pass < 2; ++pass) {
            const float thickness = pass == 0 ? 6.0f : 1.8f;
            const Color draw_color =
                pass == 0 ? Color{color.r, color.g, color.b, 28} : color;
            for (int i = 0; i < kStringNodes - 1; ++i) {
                const float r0 = kOuterRadius - static_cast<float>(i) * kStringSpacing;
                const float r1 = kOuterRadius - static_cast<float>(i + 1) * kStringSpacing;
                const Vector2 p0 = Project(layout, r0 * ux, r0 * uy, values[i]);
                const Vector2 p1 = Project(layout, r1 * ux, r1 * uy, values[i + 1]);
                DrawLineEx(p0, p1, thickness, draw_color);
            }
        }
    }
}

void DrawEndpoints(const Layout& layout) {
    constexpr int kEnds = 8;
    for (int i = 0; i < kEnds; ++i) {
        const float angle = kTwoPi * static_cast<float>(i) / 8.0f;
        DrawCircleV(Project(layout, kOuterRadius * std::cos(angle),
                            kOuterRadius * std::sin(angle), 0.0f),
                    3.0f, Color{150, 172, 200, 255});
    }
    DrawCircleV(Project(layout, 0.0f, 0.0f, 0.0f), 3.0f,
                Color{210, 225, 245, 220});
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

Rectangle ProbeBox(int width) {
    return Rectangle{static_cast<float>(width) - 168.0f, 70.0f, 144.0f, 144.0f};
}

Rectangle ProbeBoxInner(const Rectangle& box) {
    return Rectangle{box.x + 10.0f, box.y + 10.0f, box.width - 20.0f,
                     box.height - 20.0f};
}

void DrawProbeBox(const Rectangle& box, int string_index, int node) {
    const Rectangle inner = ProbeBoxInner(box);
    const float mid_x = inner.x + inner.width * 0.5f;
    const float mid_y = inner.y + inner.height * 0.5f;
    DrawRectangleRec(box, Color{6, 8, 13, 225});
    DrawRectangleLinesEx(box, 1.0f, Color{44, 60, 82, 230});
    DrawText("xy", static_cast<int>(box.x) + 8, static_cast<int>(box.y) + 6, 12,
             Color{100, 120, 145, 220});

    for (int s = 0; s < kStringCount; ++s) {
        const float angle = kTwoPi * static_cast<float>(s) / 8.0f;
        const float dx = std::cos(angle) * inner.width * 0.5f;
        const float dy = -std::sin(angle) * inner.height * 0.5f;
        DrawLineEx(Vector2{mid_x - dx, mid_y - dy},
                   Vector2{mid_x + dx, mid_y + dy}, 1.0f,
                   Color{StringColor(s).r, StringColor(s).g,
                         StringColor(s).b, 120});
    }
    DrawCircleLinesV(Vector2{mid_x, mid_y}, inner.width * 0.5f,
                     Color{255, 225, 175, 170});
    for (int terminal = 0; terminal < 8; ++terminal) {
        const float angle = kTwoPi * static_cast<float>(terminal) / 8.0f;
        DrawCircleV(Vector2{mid_x + std::cos(angle) * inner.width * 0.5f,
                            mid_y - std::sin(angle) * inner.height * 0.5f},
                    2.0f, Color{255, 235, 195, 170});
    }
    for (int fixed = 0; fixed < 4; ++fixed) {
        const int k = fixed * kRingFixedStride + kRingFixedOffset;
        const float angle = kTwoPi * static_cast<float>(k) / kRingNodes;
        DrawCircleV(Vector2{mid_x + std::cos(angle) * inner.width * 0.5f,
                            mid_y - std::sin(angle) * inner.height * 0.5f},
                    3.0f, Color{255, 150, 110, 220});
    }

    float px = mid_x;
    float py = mid_y;
    if (string_index >= 0) {
        const float angle = kTwoPi * static_cast<float>(string_index) / 8.0f;
        const float radius =
            kOuterRadius - static_cast<float>(node) * kStringSpacing;
        px = mid_x + radius / kOuterRadius * std::cos(angle) * inner.width * 0.5f;
        py = mid_y - radius / kOuterRadius * std::sin(angle) * inner.height * 0.5f;
    }
    DrawCircleV(Vector2{px, py}, 4.0f, ProbeColor(string_index));
    DrawCircleLinesV(Vector2{px, py}, 7.0f, Color{255, 255, 255, 170});
}

void DrawProbeScene(const RingDiameters& model, const Layout& layout, int string_index,
                    int node) {
    if (string_index < 0) {
        const float angle = kTwoPi * static_cast<float>(node) / kRingNodes;
        const float z = model.Ring()[node];
        const Vector2 point = Project(layout, kOuterRadius * std::cos(angle),
                                      kOuterRadius * std::sin(angle), z);
        const Vector2 base = Project(layout, kOuterRadius * std::cos(angle),
                                     kOuterRadius * std::sin(angle), 0.0f);
        const float top = std::min(base.y, point.y);
        const float bottom = std::max(base.y, point.y);
        for (float y = top + 4.0f; y + 4.0f < bottom; y += 10.0f) {
            DrawLineV(Vector2{point.x, y}, Vector2{point.x, y + 4.0f},
                      Color{255, 230, 175, 80});
        }
        DrawCircleV(point, 3.0f, ProbeColor(-1));
        DrawCircleLinesV(point, 8.0f, ProbeColor(-1));
        return;
    }

    const float angle = kTwoPi * static_cast<float>(string_index) / 8.0f;
    const float ux = std::cos(angle);
    const float uy = std::sin(angle);
    const float radius = kOuterRadius - static_cast<float>(node) * kStringSpacing;
    const float z = model.String(string_index)[node];
    const Vector2 point = Project(layout, radius * ux, radius * uy, z);
    const Vector2 base = Project(layout, radius * ux, radius * uy, 0.0f);
    const float top = std::min(base.y, point.y);
    const float bottom = std::max(base.y, point.y);
    const Color color = ProbeColor(string_index);
    for (float y = top + 4.0f; y + 4.0f < bottom; y += 10.0f) {
        DrawLineV(Vector2{point.x, y}, Vector2{point.x, y + 4.0f},
                  Color{color.r, color.g, color.b, 80});
    }
    DrawCircleV(point, 3.0f, color);
    DrawCircleLinesV(point, 7.0f, color);
}

void DrawSceneHud(const RingDiameters& model) {
    const std::string line = std::format(
        "t {:.2f} s      R: reset impulse      drag xy: probe      E: export wav",
        model.Time());
    DrawText(line.c_str(), 24, 20, 15, Color{120, 140, 165, 230});
}

std::string ProbeTitle(int string_index, int node) {
    if (string_index < 0) {
        return std::format("PROBE   RING   node {}", node);
    }
    return std::format("PROBE   S{}   T {:.1f} N   node {}", string_index,
                       kStringTensions[string_index], node);
}

void DrawScope(const ProbeHistory& history, const RingDiameters& model,
               int string_index, int node, float& scale) {
    const Layout layout = MakeLayout(GetScreenWidth(), GetScreenHeight());
    const Rectangle panel = layout.scope;
    const Rectangle plot{panel.x + 18.0f, panel.y + 44.0f, panel.width - 36.0f,
                         panel.height - 76.0f};
    const float mid = plot.y + plot.height * 0.5f;
    const int count = std::min(history.Size(), kHistoryWindowSamples);
    const float value = string_index < 0 ? model.Ring()[node]
                                         : model.String(string_index)[node];

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

    const std::string title = ProbeTitle(string_index, node);
    DrawText(title.c_str(), static_cast<int>(panel.x) + 16,
             static_cast<int>(panel.y) + 12, 15, Color{200, 215, 235, 255});

    const std::string live = std::format("z {:+.4f} m      scale {:.4f} m",
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
        const std::string label = std::format("-{}s", 4 - k);
        const int label_width = MeasureText(label.c_str(), 13);
        const float x = plot.x + plot.width * static_cast<float>(k) / 4.0f;
        DrawText(label.c_str(), static_cast<int>(x) - label_width / 2, label_y, 13,
                 Color{100, 120, 145, 220});
    }
    DrawText("now",
             static_cast<int>(plot.x + plot.width - MeasureText("now", 13)),
             label_y, 13, Color{100, 120, 145, 220});

    if (count >= 2) {
        const Color trace_color = ProbeColor(string_index);
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
        RingDiameters model;
        model.SetDamping(0.0f);
        model.ResetImpulse();
        for (int step = 0; step < kStepsPerSecond / 50; ++step) {
            model.Step(kFixedDt);
        }
        const float energy0 = model.TotalEnergy();
        float peak = 0.0f;
        for (int step = 0; step < kStepsPerSecond; ++step) {
            model.Step(kFixedDt);
            peak = std::max(peak, model.PeakDisplacement());
        }
        const float energy1 = model.TotalEnergy();
        const float drift = std::abs(energy1 - energy0) / energy0;
        check(model.IsFinite(), "finite state", "1 s undamped run");
        check(peak > 0.0005f && peak < 0.5f, "impulse response",
              std::format("peak {:.5f} m", peak));
        check(model.MaxStringSpread() > 1e-6f, "detuned string spread",
              std::format("max spread {:.2e} m", model.MaxStringSpread()));
        check(model.MaxRingFixedDrift() < 1e-6f, "ring anchors fixed",
              std::format("max drift {:.2e} m", model.MaxRingFixedDrift()));
        check(drift < 0.02f, "undamped energy drift",
              std::format("relative drift {:.4f}", drift));
        check(model.MaxRingDisplacement() > 0.0001f, "ring vibration",
              std::format("ring peak {:.5f} m", model.MaxRingDisplacement()));
    }

    {
        RingDiameters model;
        model.SetDamping(0.15f);
        model.ResetImpulse();
        for (int step = 0; step < kStepsPerSecond / 50; ++step) {
            model.Step(kFixedDt);
        }
        const float energy0 = model.TotalEnergy();
        for (int step = 0; step < kStepsPerSecond * 4; ++step) {
            model.Step(kFixedDt);
        }
        const float energy1 = model.TotalEnergy();
        check(energy1 > 0.0f && energy1 < 0.8f * energy0, "damped energy decay",
              std::format("E {:.4f} -> {:.4f} mJ", energy0 * 1000.0f,
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
    InitWindow(1280, 900, "Four Strings Ring");
    SetWindowMinSize(960, 720);
    SetTargetFPS(60);

    RingDiameters model;
    ProbeHistory history;
    int probe_string = -1;
    int probe_node = kRingTerminalStep;
    float wave_scale = 0.01f;
    int sample_countdown = kStepsPerSample;
    float accumulator = 0.0f;
    bool probe_dragging = false;
    float speed = kSpeedDefault;
    bool speed_dragging = false;

    while (!WindowShouldClose()) {
        const int width = GetScreenWidth();
        const int height = GetScreenHeight();
        const Layout layout = MakeLayout(width, height);
        const Vector2 mouse{static_cast<float>(GetMouseX()),
                            static_cast<float>(GetMouseY())};

        if (IsKeyPressed(KEY_R)) {
            model.ResetImpulse();
            history.Clear();
            wave_scale = 0.01f;
            sample_countdown = kStepsPerSample;
            accumulator = 0.0f;
        }
        if (IsKeyPressed(KEY_E)) {
            std::cout << std::format("wav export: rendering {:.0f} s...\n",
                                     kExportSeconds)
                      << std::flush;
            RingDiameters export_model;
            export_model.SetDamping(kExportDamping);
            export_model.ResetImpulse();
            const bool ok = ExportProbeWav(
                "vibrate_ring.wav", export_model, kExportStepsPerSecond1D,
                kExportSeconds,
                [](RingDiameters& m, float dt) { m.Step(dt); },
                [&](const RingDiameters& m) {
                    return probe_string < 0
                               ? m.Ring()[probe_node]
                               : m.String(probe_string)[probe_node];
                });
            std::cout << std::format("wav export: {}\n",
                                     ok ? "vibrate_ring.wav" : "failed");
        }

        const Rectangle box = ProbeBox(width);
        const Rectangle inner = ProbeBoxInner(box);
        const Rectangle speed_track = SpeedSliderTrack(width);
        const Rectangle speed_hit = SpeedSliderHit(width);

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) &&
            CheckCollisionPointRec(mouse, speed_hit)) {
            speed_dragging = true;
        } else if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) &&
                   CheckCollisionPointRec(mouse, box)) {
            probe_dragging = true;
        }

        if (speed_dragging) {
            const bool near_slider =
                mouse.y >= 20.0f && mouse.y <= 62.0f &&
                mouse.x >= speed_track.x - 20.0f &&
                mouse.x <= speed_track.x + speed_track.width + 20.0f;
            if (!near_slider) {
                speed_dragging = false;
            } else {
                const float fraction = std::clamp(
                    (mouse.x - speed_track.x) / speed_track.width, 0.0f, 1.0f);
                speed = kSpeedMin + fraction * (kSpeedMax - kSpeedMin);
                if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                    speed_dragging = false;
                }
            }
        }

        if (probe_dragging) {
            const float mid_x = inner.x + inner.width * 0.5f;
            const float mid_y = inner.y + inner.height * 0.5f;
            const float wx = (mouse.x - mid_x) / (inner.width * 0.5f) * kOuterRadius;
            const float wy = (mid_y - mouse.y) / (inner.height * 0.5f) * kOuterRadius;
            const float radius = std::sqrt(wx * wx + wy * wy);
            int new_string = probe_string;
            int new_node = probe_node;
            if (radius > 0.88f * kOuterRadius) {
                float angle = std::atan2(wy, wx);
                if (angle < 0.0f) {
                    angle += kTwoPi;
                }
                new_string = -1;
                new_node = static_cast<int>(std::lround(angle / kTwoPi * kRingNodes)) %
                           kRingNodes;
            } else {
                float best = 1e9f;
                int best_string = 0;
                float best_t = 0.0f;
                for (int s = 0; s < kStringCount; ++s) {
                    const float angle = kTwoPi * static_cast<float>(s) / 8.0f;
                    const float ux = std::cos(angle);
                    const float uy = std::sin(angle);
                    const float t = wx * ux + wy * uy;
                    const float distance = std::abs(-wx * uy + wy * ux);
                    if (distance < best) {
                        best = distance;
                        best_string = s;
                        best_t = t;
                    }
                }
                new_string = best_string;
                new_node = std::clamp(
                    static_cast<int>(std::lround(
                        (kOuterRadius - best_t) / kStringSpacing)),
                    1, kStringNodes - 2);
            }
            if (new_string != probe_string || new_node != probe_node) {
                probe_string = new_string;
                probe_node = new_node;
                history.Clear();
                wave_scale = 0.01f;
                sample_countdown = kStepsPerSample;
            }
            if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                probe_dragging = false;
            }
        }

        accumulator += std::min(GetFrameTime(), kMaxFrameTime) * speed;
        int steps = 0;
        while (accumulator >= kFixedDt && steps < kMaxStepsPerFrame) {
            model.Step(kFixedDt);
            if (--sample_countdown <= 0) {
                const float value =
                    probe_string < 0 ? model.Ring()[probe_node]
                                     : model.String(probe_string)[probe_node];
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
        DrawStrings(model, layout);
        DrawRing(model, layout);
        DrawEndpoints(layout);
        DrawProbeScene(model, layout, probe_string, probe_node);
        EndScissorMode();

        DrawProbeBox(box, probe_string, probe_node);
        DrawSceneHud(model);
        DrawSpeedSlider(width, speed);
        DrawScope(history, model, probe_string, probe_node, wave_scale);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
