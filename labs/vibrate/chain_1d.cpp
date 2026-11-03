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

// 一维张紧弦：N 个质点沿 x 排布，只允许 y 方向位移。
// 采用小振幅弦方程 a = (T / mu) * y_xx；网格足够细，且不引入随振幅变化的波速。
namespace qwqdsp_vibrate {

constexpr int kNodeCount = 401;
constexpr float kRestLength = 3.2f;  // m
constexpr float kSpacing = kRestLength / static_cast<float>(kNodeCount - 1);
constexpr float kLinearDensity = 0.25f;  // kg/m
constexpr float kNodeMass = kLinearDensity * kSpacing;
constexpr float kTensionMin = 10.0f;
constexpr float kTensionMax = 300.0f;
constexpr float kTensionStep = 10.0f;
constexpr float kDampingMin = 0.0f;
constexpr float kDampingMax = 4.0f;
constexpr float kDampingStep = 0.05f;
constexpr int kStepsPerSecond = 16000;
constexpr float kFixedDt = 1.0f / static_cast<float>(kStepsPerSecond);
constexpr float kMaxFrameTime = 0.05f;
constexpr int kMaxStepsPerFrame = 900;
constexpr float kMaxPluck = 0.25f;
constexpr float kDriveAmplitude = 0.02f;
constexpr int kDriveNode = kNodeCount / 3;
constexpr float kHistorySeconds = 4.0f;
constexpr int kHistoryRate = 500;
constexpr int kHistorySize = 2048;
constexpr int kHistoryWindowSamples = static_cast<int>(kHistorySeconds * kHistoryRate);
constexpr int kStepsPerSample = kStepsPerSecond / kHistoryRate;

constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kTwoPi = 2.0f * kPi;

// 波速 c = sqrt(T / mu)，线密度 mu = m / spacing。
float WaveSpeed(float tension) {
    return std::sqrt(tension * kSpacing / kNodeMass);
}

// 两端固定时的一阶固有频率。
float Fundamental(float tension) {
    return WaveSpeed(tension) / (2.0f * kRestLength);
}

class Chain1D {
public:
    Chain1D() : y_(kNodeCount, 0.0f), v_(kNodeCount, 0.0f), accel_(kNodeCount, 0.0f) {
        Reset();
    }

    void Reset() {
        std::fill(y_.begin(), y_.end(), 0.0f);
        std::fill(v_.begin(), v_.end(), 0.0f);
        std::fill(accel_.begin(), accel_.end(), 0.0f);
        time_ = 0.0f;
        drive_phase_ = 0.0f;
    }

    void Step(float dt) {
        // 小振幅近似下的标准弦方程：a = (T / mu) * y_xx。相比逐段张力模型，
        // 它没有随振幅变化的波速，不会在大幅度抓取时产生非物理的高次谐波。
        const float wave_speed_sq = tension_ / kLinearDensity;
        const float inv_spacing_sq = 1.0f / (kSpacing * kSpacing);
        for (int i = 1; i < kNodeCount - 1; ++i) {
            const float curvature = y_[i + 1] - 2.0f * y_[i] + y_[i - 1];
            accel_[i] = wave_speed_sq * curvature * inv_spacing_sq - damping_ * v_[i];
        }
        for (int i = 1; i < kNodeCount - 1; ++i) {
            v_[i] += accel_[i] * dt;
            y_[i] += v_[i] * dt;
        }
        time_ += dt;

        if (driving_) {
            const float omega = kTwoPi * Fundamental(tension_);
            drive_phase_ = std::fmod(drive_phase_ + omega * dt, kTwoPi);
            y_[kDriveNode] = kDriveAmplitude * std::sin(drive_phase_);
            v_[kDriveNode] = kDriveAmplitude * omega * std::cos(drive_phase_);
        }
    }

    // 三角形拨弦：峰位于 length 比例 position 处。
    void PluckTriangle(float height, float position) {
        Reset();
        const int peak = std::clamp(static_cast<int>(std::lround(position * (kNodeCount - 1))),
                                    1, kNodeCount - 2);
        for (int i = 1; i < kNodeCount - 1; ++i) {
            const float t = i <= peak
                                ? static_cast<float>(i) / static_cast<float>(peak)
                                : static_cast<float>(kNodeCount - 1 - i) /
                                      static_cast<float>(kNodeCount - 1 - peak);
            y_[i] = height * t;
        }
    }

    // 鼠标抓取：按固定端到抓取点的静力平衡形状（三角形）拉弦。
    // 只压一个局部窗口会形成很陡的凹坑，释放后把大量能量灌进高次谐波。
    void Grab(int center, float displacement) {
        const float left_span = static_cast<float>(std::max(center, 1));
        const float right_span = static_cast<float>(std::max(kNodeCount - 1 - center, 1));
        for (int i = 1; i < kNodeCount - 1; ++i) {
            const float t = i <= center
                                ? static_cast<float>(i) / left_span
                                : static_cast<float>(kNodeCount - 1 - i) / right_span;
            y_[i] = displacement * t;
            v_[i] = 0.0f;
        }
    }

    float TotalEnergy() const {
        float kinetic = 0.0f;
        for (int i = 1; i < kNodeCount - 1; ++i) {
            kinetic += 0.5f * kNodeMass * v_[i] * v_[i];
        }
        float length = 0.0f;
        for (int i = 0; i < kNodeCount - 1; ++i) {
            const float dy = y_[i + 1] - y_[i];
            length += std::sqrt(kSpacing * kSpacing + dy * dy);
        }
        const float potential = tension_ * (length - kRestLength);
        return kinetic + potential;
    }

    float PeakDisplacement() const {
        float peak = 0.0f;
        for (float value : y_) {
            peak = std::max(peak, std::abs(value));
        }
        return peak;
    }

    bool IsFinite() const {
        for (int i = 0; i < kNodeCount; ++i) {
            if (!std::isfinite(y_[i]) || !std::isfinite(v_[i])) {
                return false;
            }
        }
        return true;
    }

    void SetTension(float tension) {
        tension_ = std::clamp(tension, kTensionMin, kTensionMax);
    }

    void SetDamping(float damping) {
        damping_ = std::clamp(damping, kDampingMin, kDampingMax);
    }

    void SetDriving(bool driving) { driving_ = driving; }

    float Tension() const { return tension_; }
    float Damping() const { return damping_; }
    float Time() const { return time_; }
    const std::vector<float>& Displacement() const { return y_; }
    const std::vector<float>& Velocity() const { return v_; }

private:
    std::vector<float> y_;
    std::vector<float> v_;
    std::vector<float> accel_;
    float tension_ = 80.0f;
    float damping_ = 0.12f;
    float time_ = 0.0f;
    float drive_phase_ = 0.0f;
    bool driving_ = false;
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

    // age 0 是最近一次采样。
    float At(int age) const {
        return samples_[(head_ - 1 - age + kHistorySize * 2) % kHistorySize];
    }

    int Size() const { return filled_; }

private:
    std::array<float, kHistorySize> samples_{};
    int head_ = 0;
    int filled_ = 0;
};

struct View {
    float left = 0.0f;
    float right = 0.0f;
    float center_y = 0.0f;
    float pixels_per_meter = 600.0f;

    float NodeX(int index) const {
        const float t = static_cast<float>(index) / static_cast<float>(kNodeCount - 1);
        return left + (right - left) * t;
    }

    float ScreenY(float displacement) const {
        return center_y - displacement * pixels_per_meter;
    }

    float Displacement(float screen_y) const {
        return (center_y - screen_y) / pixels_per_meter;
    }
};

struct Layout {
    View view;
    float string_top = 185.0f;
    float slider_y = 0.0f;
    Rectangle waveform{};
};

Layout MakeLayout(int width, int height) {
    Layout layout;
    const float panel_height =
        std::clamp(static_cast<float>(height) * 0.26f, 150.0f, 240.0f);
    const float panel_bottom = static_cast<float>(height) - 44.0f;
    const float panel_top = panel_bottom - panel_height;

    layout.string_top = 185.0f;
    layout.slider_y = panel_top - 26.0f;
    layout.view.left = static_cast<float>(width) * 0.08f;
    layout.view.right = static_cast<float>(width) * 0.92f;
    layout.view.center_y = (layout.string_top + layout.slider_y) * 0.5f;
    const float half_span = (layout.slider_y - layout.string_top) * 0.5f * 0.88f;
    layout.view.pixels_per_meter = std::max(half_span / kMaxPluck, 1.0f);

    layout.waveform = Rectangle{static_cast<float>(width) * 0.07f, panel_top,
                                static_cast<float>(width) * 0.86f, panel_height};
    return layout;
}

void DrawBackdrop(const View& view, int width, int height) {
    DrawRectangleGradientV(0, 0, width, height, Color{19, 25, 42, 255}, Color{7, 9, 16, 255});

    const Color grid{28, 38, 56, 80};
    for (int x = 0; x < width; x += 48) {
        DrawLineV(Vector2{static_cast<float>(x), 0.0f},
                  Vector2{static_cast<float>(x), static_cast<float>(height)}, grid);
    }
    for (int y = 0; y < height; y += 48) {
        DrawLineV(Vector2{0.0f, static_cast<float>(y)},
                  Vector2{static_cast<float>(width), static_cast<float>(y)}, grid);
    }

    // 静止位置虚线。
    const float rest_y = view.ScreenY(0.0f);
    for (float x = view.left; x < view.right; x += 18.0f) {
        const float x_end = std::min(x + 9.0f, view.right);
        DrawLineV(Vector2{x, rest_y}, Vector2{x_end, rest_y}, Color{70, 95, 130, 150});
    }

    // 左侧 0.1 m 刻度。
    const Color tick{80, 103, 135, 160};
    for (int k = -2; k <= 2; ++k) {
        if (k == 0) {
            continue;
        }
        const float y = view.ScreenY(0.1f * static_cast<float>(k));
        DrawLineV(Vector2{view.left - 14.0f, y}, Vector2{view.left - 5.0f, y}, tick);
    }
}

void DrawString(const Chain1D& chain, const View& view) {
    const auto& ys = chain.Displacement();
    const std::array<float, 3> thickness{13.0f, 6.0f, 2.2f};
    const std::array<Color, 3> colors{Color{45, 110, 255, 22}, Color{90, 170, 255, 48},
                                      Color{195, 235, 255, 255}};
    for (int pass = 0; pass < 3; ++pass) {
        for (int i = 0; i < kNodeCount - 1; ++i) {
            DrawLineEx(Vector2{view.NodeX(i), view.ScreenY(ys[i])},
                       Vector2{view.NodeX(i + 1), view.ScreenY(ys[i + 1])}, thickness[pass],
                       colors[pass]);
        }
    }
}

void DrawNodes(const Chain1D& chain, const View& view) {
    const auto& ys = chain.Displacement();
    const auto& vs = chain.Velocity();
    for (int i = 3; i < kNodeCount - 3; i += 4) {
        const float speed = std::clamp(std::abs(vs[i]) / 2.5f, 0.0f, 1.0f);
        const Color color = ColorLerp(Color{120, 215, 255, 220}, Color{255, 145, 90, 255}, speed);
        DrawCircleV(Vector2{view.NodeX(i), view.ScreenY(ys[i])}, 2.2f + 1.8f * speed, color);
    }
}

void DrawAnchors(const View& view) {
    const float y = view.ScreenY(0.0f);
    for (int index : {0, kNodeCount - 1}) {
        const float x = view.NodeX(index);
        DrawRectangle(static_cast<int>(x) - 2, static_cast<int>(y) - 18, 4, 36,
                      Color{55, 75, 105, 255});
        DrawCircleV(Vector2{x, y}, 7.0f, Color{35, 55, 85, 255});
        DrawCircleV(Vector2{x, y}, 3.5f, Color{175, 215, 255, 255});
    }
}

void DrawDriveMarker(const View& view) {
    const float x = view.NodeX(kDriveNode);
    const float y = view.ScreenY(0.0f);
    DrawCircleV(Vector2{x, y}, 9.0f, Color{90, 230, 160, 40});
    DrawCircleLinesV(Vector2{x, y}, 9.0f, Color{110, 240, 175, 180});
}

void DrawProbeGuide(const Chain1D& chain, const View& view, float slider_y, int probe) {
    const float x = view.NodeX(probe);
    const float node_y = view.ScreenY(chain.Displacement()[probe]);
    for (float y = node_y + 12.0f; y < slider_y - 10.0f; y += 16.0f) {
        DrawLineV(Vector2{x, y}, Vector2{x, std::min(y + 8.0f, slider_y - 10.0f)},
                  Color{255, 185, 100, 45});
    }
    DrawCircleLinesV(Vector2{x, node_y}, 8.0f, Color{255, 195, 120, 90});
    DrawCircleLinesV(Vector2{x, node_y}, 5.5f, Color{255, 200, 130, 220});
}

void DrawProbeSlider(const Layout& layout, int probe) {
    const View& view = layout.view;
    const float y = layout.slider_y;
    DrawLineEx(Vector2{view.left, y}, Vector2{view.right, y}, 2.0f, Color{38, 54, 78, 255});
    for (int i = 0; i <= kNodeCount - 1; i += 20) {
        const float x = view.NodeX(i);
        DrawLineV(Vector2{x, y + 3.0f}, Vector2{x, y + 8.0f}, Color{46, 64, 92, 200});
    }
    const float x = view.NodeX(probe);
    DrawCircleV(Vector2{x, y}, 11.0f, Color{255, 180, 90, 30});
    DrawCircleV(Vector2{x, y}, 6.5f, Color{255, 192, 108, 255});
    DrawCircleV(Vector2{x, y}, 2.4f, Color{48, 34, 20, 255});
}

void DrawProbeWaveform(const ProbeHistory& history, const Chain1D& chain, const Layout& layout,
                       int probe, float& scale, bool running) {
    const Rectangle panel = layout.waveform;
    const Rectangle plot{panel.x + 18.0f, panel.y + 44.0f, panel.width - 36.0f,
                         panel.height - 76.0f};
    const float mid = plot.y + plot.height * 0.5f;
    const int count = std::min(history.Size(), kHistoryWindowSamples);

    float peak = 0.0f;
    for (int age = 0; age < count; ++age) {
        peak = std::max(peak, std::abs(history.At(age)));
    }
    if (running && count > 0) {
        const float target = std::max(peak * 1.15f, 0.004f);
        if (target > scale) {
            scale = target;
        } else {
            scale += (target - scale) * 0.04f;
        }
    }

    DrawRectangleRounded(panel, 0.05f, 8, Color{8, 12, 21, 225});
    DrawRectangleRoundedLinesEx(panel, 0.05f, 8, 1.0f, Color{50, 76, 110, 170});

    const std::string title = std::format("PROBE WAVEFORM      node {}   x {:.2f} m", probe,
                                          static_cast<float>(probe) * kSpacing);
    DrawText(title.c_str(), static_cast<int>(panel.x) + 18, static_cast<int>(panel.y) + 13, 18,
             Color{255, 200, 140, 255});

    const std::string live = std::format("y {:+.3f} m      scale {:.3f} m",
                                         chain.Displacement()[probe], scale);
    const int live_width = MeasureText(live.c_str(), 16);
    DrawText(live.c_str(), static_cast<int>(panel.x + panel.width) - 18 - live_width,
             static_cast<int>(panel.y) + 14, 16, Color{150, 170, 195, 255});

    DrawRectangleRec(plot, Color{5, 8, 14, 220});
    DrawRectangleLinesEx(plot, 1.0f, Color{38, 56, 82, 200});

    for (float x = plot.x; x < plot.x + plot.width; x += 12.0f) {
        DrawLineV(Vector2{x, mid}, Vector2{std::min(x + 6.0f, plot.x + plot.width), mid},
                  Color{60, 80, 110, 120});
    }
    for (int i = 1; i <= 3; ++i) {
        const float x = plot.x + plot.width * static_cast<float>(i) / 4.0f;
        DrawLineV(Vector2{x, plot.y + 2.0f}, Vector2{x, plot.y + plot.height - 2.0f},
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
    const char* now_label = "now";
    DrawText(now_label, static_cast<int>(plot.x + plot.width - MeasureText(now_label, 13)),
             label_y, 13, Color{100, 120, 145, 220});

    if (count >= 2) {
        const Color color{255, 190, 110, 235};
        Vector2 previous{};
        bool has_previous = false;
        for (int age = count - 1; age >= 0; --age) {
            const float x =
                plot.x + plot.width *
                             (1.0f - static_cast<float>(age) /
                                         static_cast<float>(kHistoryWindowSamples - 1));
            const float value = history.At(age);
            const float y = std::clamp(mid - value / scale * (plot.height * 0.5f - 4.0f),
                                       plot.y + 1.0f, plot.y + plot.height - 1.0f);
            const Vector2 current{x, y};
            if (has_previous) {
                DrawLineEx(previous, current, 1.7f, color);
            }
            previous = current;
            has_previous = true;
        }
        const float latest = std::clamp(mid - history.At(0) / scale * (plot.height * 0.5f - 4.0f),
                                        plot.y + 1.0f, plot.y + plot.height - 1.0f);
        DrawCircleV(Vector2{plot.x + plot.width, latest}, 3.0f, color);
    } else {
        const char* waiting = "accumulating...";
        const int waiting_width = MeasureText(waiting, 14);
        DrawText(waiting, static_cast<int>(plot.x + plot.width * 0.5f) - waiting_width / 2,
                 static_cast<int>(mid) - 7, 14, Color{90, 110, 135, 255});
    }
}

std::string FormatMeters(float value) {
    return std::format("{:.3f}", value);
}

void DrawHud(const Chain1D& chain, bool running, bool driving) {
    const int width = GetScreenWidth();
    const int height = GetScreenHeight();

    const std::string title = "1D VIBRATION / TENSIONED STRING";
    const std::string line1 = std::format("tension {:.0f} N      damping {:.2f} 1/s",
                                          chain.Tension(), chain.Damping());
    const std::string line2 =
        std::format("f1 {:.2f} Hz      wave speed {:.1f} m/s      nodes {}", Fundamental(chain.Tension()),
                    WaveSpeed(chain.Tension()), kNodeCount);
    const std::string line3 = std::format("energy {:.1f} mJ      sim {:.2f} s      peak {} m",
                                          chain.TotalEnergy() * 1000.0f, chain.Time(),
                                          FormatMeters(chain.PeakDisplacement()));

    const int panel_width = std::max({MeasureText(title.c_str(), 19),
                                      MeasureText(line2.c_str(), 17),
                                      MeasureText(line3.c_str(), 17)}) +
                            44;
    const Rectangle panel{24.0f, 20.0f, static_cast<float>(panel_width), 142.0f};
    DrawRectangleRounded(panel, 0.14f, 8, Color{8, 12, 20, 205});
    DrawRectangleRoundedLinesEx(panel, 0.14f, 8, 1.0f, Color{55, 86, 122, 170});

    const int text_x = 44;
    DrawText(title.c_str(), text_x, 34, 19, Color{155, 210, 255, 255});
    DrawText(line1.c_str(), text_x, 66, 17, Color{150, 166, 190, 255});
    DrawText(line2.c_str(), text_x, 92, 17, Color{150, 166, 190, 255});
    DrawText(line3.c_str(), text_x, 118, 17, Color{150, 166, 190, 255});

    int badge_y = 34;
    if (!running) {
        const char* text = "PAUSED";
        const int text_width = MeasureText(text, 18);
        DrawText(text, width - 40 - text_width, badge_y, 18, Color{255, 190, 95, 255});
        badge_y += 26;
    }
    if (driving) {
        const char* text = "DRIVE";
        const int text_width = MeasureText(text, 18);
        DrawText(text, width - 40 - text_width, badge_y, 18, Color{120, 235, 170, 255});
    }

    DrawRectangle(0, height - 38, width, 38, Color{6, 9, 15, 210});
    DrawText("drag pluck | slider probe | SPACE pause | R reset | UP/DOWN tension | "
             "LEFT/RIGHT damping | A drive | E wav | F12 shot",
             24, height - 27, 15, Color{115, 135, 160, 255});
}

int RunSelfTest() {
    int failures = 0;
    const auto check = [&failures](bool ok, const char* name, const std::string& detail) {
        std::cout << std::format("[{}] {:<26} {}\n", ok ? "PASS" : "FAIL", name, detail);
        if (!ok) {
            ++failures;
        }
    };

    {
        Chain1D chain;
        chain.SetTension(80.0f);
        chain.SetDamping(0.35f);
        chain.PluckTriangle(0.18f, 0.5f);
        const float energy0 = chain.TotalEnergy();
        float peak = 0.0f;
        for (int i = 0; i < kStepsPerSecond * 4; ++i) {
            chain.Step(kFixedDt);
            peak = std::max(peak, chain.PeakDisplacement());
        }
        const float energy1 = chain.TotalEnergy();
        check(chain.IsFinite(), "finite state", "damped 4 s run");
        check(energy1 > 0.0f && energy1 < 0.5f * energy0, "damped energy decay",
              std::format("E {:.2f} -> {:.2f} mJ", energy0 * 1000.0f, energy1 * 1000.0f));
        check(peak < 0.19f, "bounded pluck amplitude", std::format("peak {:.4f} m", peak));
    }

    {
        Chain1D chain;
        chain.SetTension(80.0f);
        chain.SetDamping(0.0f);
        chain.PluckTriangle(0.10f, 0.5f);
        const float energy0 = chain.TotalEnergy();
        for (int i = 0; i < kStepsPerSecond * 2; ++i) {
            chain.Step(kFixedDt);
        }
        const float energy1 = chain.TotalEnergy();
        const float drift = std::abs(energy1 - energy0) / energy0;
        check(drift < 0.02f, "undamped energy drift",
              std::format("relative drift {:.4f}", drift));
    }

    {
        Chain1D chain;
        chain.SetTension(kTensionMax);
        chain.SetDamping(0.0f);
        chain.PluckTriangle(0.20f, 0.4f);
        float peak = 0.0f;
        for (int i = 0; i < kStepsPerSecond * 2; ++i) {
            chain.Step(kFixedDt);
            peak = std::max(peak, chain.PeakDisplacement());
        }
        check(chain.IsFinite() && peak < 0.21f, "max tension stable",
              std::format("T {:.0f} N, peak {:.4f} m", kTensionMax, peak));
    }

    {
        Chain1D chain;
        chain.SetTension(80.0f);
        chain.SetDamping(0.35f);
        chain.Reset();
        chain.SetDriving(true);
        float peak = 0.0f;
        for (int i = 0; i < kStepsPerSecond * 5; ++i) {
            chain.Step(kFixedDt);
            peak = std::max(peak, chain.PeakDisplacement());
        }
        const float energy = chain.TotalEnergy();
        check(chain.IsFinite() && peak > 0.01f, "driven response",
              std::format("peak {:.4f} m", peak));
        check(energy < 5.0f, "driven energy bounded",
              std::format("E {:.2f} mJ", energy * 1000.0f));
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
    InitWindow(1280, 840, "1D Vibration - Tensioned String");
    SetWindowMinSize(960, 680);
    SetTargetFPS(60);

    Chain1D chain;
    chain.PluckTriangle(0.18f, 0.42f);

    bool running = true;
    bool driving = false;
    int grab_center = -1;
    float grab_value = 0.0f;
    float accumulator = 0.0f;
    int probe_node = kNodeCount / 2;
    bool probe_dragging = false;
    ProbeHistory probe_history;
    float wave_scale = 0.1f;
    int sample_countdown = kStepsPerSample;

    while (!WindowShouldClose()) {
        const Layout layout = MakeLayout(GetScreenWidth(), GetScreenHeight());
        const View& view = layout.view;

        if (IsKeyPressed(KEY_R)) {
            chain.Reset();
            accumulator = 0.0f;
            probe_history.Clear();
            wave_scale = 0.1f;
            sample_countdown = kStepsPerSample;
        }
        if (IsKeyPressed(KEY_E)) {
            std::cout << std::format("wav export: rendering {:.0f} s...\n",
                                     kExportSeconds)
                      << std::flush;
            Chain1D export_model = chain;
            export_model.SetDriving(false);
            export_model.SetDamping(kExportDamping);
            export_model.PluckTriangle(0.18f, 0.42f);
            const bool ok = ExportProbeWav(
                "vibrate_chain.wav", export_model, kExportStepsPerSecond1D,
                kExportSeconds,
                [](Chain1D& m, float dt) { m.Step(dt); },
                [&](const Chain1D& m) { return m.Displacement()[probe_node]; });
            std::cout << std::format("wav export: {}\n",
                                     ok ? "vibrate_chain.wav" : "failed");
        }
        if (IsKeyPressed(KEY_SPACE)) {
            running = !running;
        }
        if (IsKeyPressed(KEY_A)) {
            driving = !driving;
        }
        chain.SetDriving(driving);

        if (IsKeyPressed(KEY_UP) || IsKeyPressedRepeat(KEY_UP)) {
            chain.SetTension(chain.Tension() + kTensionStep);
        }
        if (IsKeyPressed(KEY_DOWN) || IsKeyPressedRepeat(KEY_DOWN)) {
            chain.SetTension(chain.Tension() - kTensionStep);
        }
        if (IsKeyPressed(KEY_RIGHT) || IsKeyPressedRepeat(KEY_RIGHT)) {
            chain.SetDamping(chain.Damping() + kDampingStep);
        }
        if (IsKeyPressed(KEY_LEFT) || IsKeyPressedRepeat(KEY_LEFT)) {
            chain.SetDamping(chain.Damping() - kDampingStep);
        }

        const float mouse_x = static_cast<float>(GetMouseX());
        const float mouse_y = static_cast<float>(GetMouseY());

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            const bool on_slider = std::abs(mouse_y - layout.slider_y) <= 14.0f &&
                                   mouse_x >= view.left - 12.0f &&
                                   mouse_x <= view.right + 12.0f;
            if (on_slider) {
                probe_dragging = true;
            } else if (view.right > view.left && mouse_y > layout.string_top - 20.0f &&
                       mouse_y < layout.slider_y - 20.0f && mouse_x >= view.left &&
                       mouse_x <= view.right) {
                const float t = (mouse_x - view.left) / (view.right - view.left);
                grab_center = std::clamp(
                    static_cast<int>(std::lround(t * static_cast<float>(kNodeCount - 1))), 1,
                    kNodeCount - 2);
                grab_value = chain.Displacement()[grab_center];
            }
        }
        if (probe_dragging && IsMouseButtonDown(MOUSE_BUTTON_LEFT) && view.right > view.left) {
            const float t =
                std::clamp((mouse_x - view.left) / (view.right - view.left), 0.0f, 1.0f);
            const int new_probe = std::clamp(
                static_cast<int>(std::lround(t * static_cast<float>(kNodeCount - 1))), 1,
                kNodeCount - 2);
            if (new_probe != probe_node) {
                probe_node = new_probe;
                probe_history.Clear();
                wave_scale = 0.05f;
                sample_countdown = kStepsPerSample;
            }
        }
        if (grab_center >= 0 && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            const float target =
                std::clamp(view.Displacement(mouse_y), -kMaxPluck, kMaxPluck);
            const float blend =
                1.0f - std::exp(-std::min(GetFrameTime(), kMaxFrameTime) / 0.045f);
            grab_value += (target - grab_value) * blend;
            chain.Grab(grab_center, grab_value);
        }
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) {
            grab_center = -1;
            probe_dragging = false;
        }

        if (IsKeyPressed(KEY_F12)) {
            TakeScreenshot("vibrate-chain.png");
        }

        if (running) {
            accumulator += std::min(GetFrameTime(), kMaxFrameTime);
            int steps = 0;
            while (accumulator >= kFixedDt && steps < kMaxStepsPerFrame) {
                chain.Step(kFixedDt);
                if (--sample_countdown <= 0) {
                    probe_history.Push(chain.Displacement()[probe_node]);
                    sample_countdown = kStepsPerSample;
                }
                accumulator -= kFixedDt;
                ++steps;
            }
        }

        BeginDrawing();
        DrawBackdrop(view, GetScreenWidth(), GetScreenHeight());
        DrawString(chain, view);
        DrawNodes(chain, view);
        DrawAnchors(view);
        DrawProbeGuide(chain, view, layout.slider_y, probe_node);
        if (driving) {
            DrawDriveMarker(view);
        }
        DrawProbeSlider(layout, probe_node);
        DrawProbeWaveform(probe_history, chain, layout, probe_node, wave_scale, running);
        DrawHud(chain, running, driving);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
