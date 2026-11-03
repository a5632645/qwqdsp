#include <raylib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <string>
#include <vector>

#include "wav_export.hpp"

// 正方形薄膜：二维波动方程，四边固定，顶视图用发散 colormap 显示位移。
namespace qwqdsp_vibrate {

constexpr int kSize = 161;
constexpr float kLength = 2.0f;
constexpr float kDx = kLength / static_cast<float>(kSize - 1);
constexpr float kWaveSpeed = 8.0f;
constexpr float kDensity = 0.25f;  // kg/m^2
constexpr float kDefaultDamping = 0.15f;

constexpr int kStepsPerSecond = 2000;
constexpr float kFixedDt = 1.0f / static_cast<float>(kStepsPerSecond);
constexpr float kMaxFrameTime = 0.05f;
constexpr int kMaxStepsPerFrame = 120;
constexpr float kSpeedMin = 0.10f;
constexpr float kSpeedMax = 1.00f;
constexpr float kSpeedDefault = 0.25f;

// 中心瞬时速度冲量；空间上收窄到几个网格单元，避免单单元 delta 激起满屏网格噪声。
constexpr float kImpulseVelocity = 100.0f;  // m/s
constexpr float kImpulseSigmaCells = 2.0f;

constexpr float kHistorySeconds = 4.0f;
constexpr int kHistoryRate = 500;
constexpr int kHistorySize = 2048;
constexpr int kHistoryWindowSamples = static_cast<int>(kHistorySeconds * kHistoryRate);
constexpr int kStepsPerSample = kStepsPerSecond / kHistoryRate;

class Membrane {
public:
    Membrane() : u_(kSize * kSize, 0.0f), v_(kSize * kSize, 0.0f),
                 a_(kSize * kSize, 0.0f) {
        BuildImpulseWeights();
        ResetImpulse();
    }

    void ResetImpulse() {
        std::fill(u_.begin(), u_.end(), 0.0f);
        std::fill(v_.begin(), v_.end(), 0.0f);
        std::fill(a_.begin(), a_.end(), 0.0f);
        for (const ImpulseCell& cell : impulse_cells_) {
            v_[cell.index] = kImpulseVelocity * cell.weight;
        }
        time_ = 0.0f;
    }

    void Step(float dt) {
        ComputeAccelerations();
        for (int j = 1; j < kSize - 1; ++j) {
            for (int i = 1; i < kSize - 1; ++i) {
                const int index = j * kSize + i;
                v_[index] += 0.5f * a_[index] * dt;
                u_[index] += v_[index] * dt;
            }
        }

        ComputeAccelerations();
        for (int j = 1; j < kSize - 1; ++j) {
            for (int i = 1; i < kSize - 1; ++i) {
                const int index = j * kSize + i;
                v_[index] += 0.5f * a_[index] * dt;
            }
        }
        time_ += dt;
    }

    float TotalEnergy() const {
        const float cell_mass = kDensity * kDx * kDx;
        float kinetic = 0.0f;
        for (float value : v_) {
            kinetic += 0.5f * cell_mass * value * value;
        }

        float potential = 0.0f;
        for (int j = 1; j < kSize - 1; ++j) {
            for (int i = 1; i < kSize - 1; ++i) {
                const int index = j * kSize + i;
                const float neighbors = u_[index + 1] + u_[index - 1] +
                                        u_[index + kSize] + u_[index - kSize];
                const float diagonals =
                    u_[index + kSize + 1] + u_[index + kSize - 1] +
                    u_[index - kSize + 1] + u_[index - kSize - 1];
                const float laplacian =
                    (4.0f * neighbors + diagonals - 20.0f * u_[index]) / 6.0f;
                potential -= 0.5f * kDensity * kWaveSpeed * kWaveSpeed *
                             u_[index] * laplacian;
            }
        }
        return kinetic + potential;
    }

    float PeakDisplacement() const {
        float peak = 0.0f;
        for (float value : u_) {
            peak = std::max(peak, std::abs(value));
        }
        return peak;
    }

    float MaxSymmetryError() const {
        float error = 0.0f;
        for (int j = 0; j < kSize; ++j) {
            for (int i = 0; i < kSize; ++i) {
                error = std::max(
                    error, std::abs(u_[j * kSize + i] - u_[i * kSize + j]));
            }
        }
        return error;
    }

    bool IsFinite() const {
        for (int index = 0; index < kSize * kSize; ++index) {
            if (!std::isfinite(u_[index]) || !std::isfinite(v_[index])) {
                return false;
            }
        }
        return true;
    }

    void SetDamping(float damping) { damping_ = damping; }
    float Time() const { return time_; }
    float At(int i, int j) const { return u_[j * kSize + i]; }
    const std::vector<float>& U() const { return u_; }

private:
    struct ImpulseCell {
        int index;
        float weight;
    };

    void BuildImpulseWeights() {
        const int center = (kSize - 1) / 2;
        const int radius = static_cast<int>(std::ceil(3.0f * kImpulseSigmaCells));
        for (int dj = -radius; dj <= radius; ++dj) {
            for (int di = -radius; di <= radius; ++di) {
                const float distance_sq =
                    static_cast<float>(di * di + dj * dj) /
                    (kImpulseSigmaCells * kImpulseSigmaCells);
                if (distance_sq > 9.0f) {
                    continue;
                }
                const float weight = std::exp(-0.5f * distance_sq);
                impulse_cells_.push_back(
                    ImpulseCell{(center + dj) * kSize + (center + di), weight});
            }
        }
    }

    void ComputeAccelerations() {
        const float wave_speed_sq = kWaveSpeed * kWaveSpeed;
        const float inv_spacing_sq = 1.0f / (kDx * kDx);
        for (int j = 1; j < kSize - 1; ++j) {
            for (int i = 1; i < kSize - 1; ++i) {
                const int index = j * kSize + i;
                // 九点各向同性拉普拉斯，抑制方形网格的方向性色散。
                const float neighbors = u_[index + 1] + u_[index - 1] +
                                        u_[index + kSize] + u_[index - kSize];
                const float diagonals =
                    u_[index + kSize + 1] + u_[index + kSize - 1] +
                    u_[index - kSize + 1] + u_[index - kSize - 1];
                const float laplacian =
                    (4.0f * neighbors + diagonals - 20.0f * u_[index]) / 6.0f;
                a_[index] = wave_speed_sq * inv_spacing_sq * laplacian -
                            damping_ * v_[index];
            }
        }
    }

    std::vector<float> u_;
    std::vector<float> v_;
    std::vector<float> a_;
    std::vector<ImpulseCell> impulse_cells_;
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
    Rectangle membrane{};
    Rectangle scope{};
};

Layout MakeLayout(int width, int height) {
    Layout layout;
    const float w = static_cast<float>(width);
    const float h = static_cast<float>(height);
    const float scene_height = h * 0.63f;
    const float size = std::min(w * 0.50f, scene_height * 0.82f);
    const float cx = w * 0.5f;
    const float cy = scene_height * 0.54f;
    layout.scene = Rectangle{0.0f, 0.0f, w, scene_height};
    layout.membrane = Rectangle{cx - size * 0.5f, cy - size * 0.5f, size, size};
    const float scope_top = scene_height + 22.0f;
    layout.scope = Rectangle{w * 0.07f, scope_top, w * 0.86f,
                             h - 24.0f - scope_top};
    return layout;
}

Color ColorFromValue(float t) {
    t = std::clamp(t, -1.0f, 1.0f);
    const Color zero{18, 22, 32, 255};
    if (t < 0.0f) {
        return ColorLerp(zero, Color{60, 130, 255, 255}, -t);
    }
    return ColorLerp(zero, Color{255, 150, 70, 255}, t);
}

void FillColormap(const Membrane& model, float scale, Image& image) {
    auto* pixels = reinterpret_cast<Color*>(image.data);
    const std::vector<float>& values = model.U();
    for (int index = 0; index < kSize * kSize; ++index) {
        pixels[index] = ColorFromValue(values[index] / scale);
    }
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

void DrawSceneHud(const Membrane& model) {
    const std::string line = std::format(
        "t {:.2f} s      R: reset impulse      click/drag: probe      E: export wav",
        model.Time());
    DrawText(line.c_str(), 24, 20, 15, Color{120, 140, 165, 230});
}

void DrawProbeMarker(const Rectangle& membrane, int i, int j) {
    const float x = membrane.x + (static_cast<float>(i) + 0.5f) / kSize * membrane.width;
    const float y = membrane.y + (static_cast<float>(j) + 0.5f) / kSize * membrane.height;
    const Vector2 p{x, y};
    DrawCircleLinesV(p, 8.0f, Color{0, 0, 0, 180});
    DrawCircleLinesV(p, 7.0f, Color{255, 255, 255, 220});
    DrawLineV(Vector2{p.x - 10.0f, p.y}, Vector2{p.x + 10.0f, p.y},
              Color{0, 0, 0, 150});
    DrawLineV(Vector2{p.x, p.y - 10.0f}, Vector2{p.x, p.y + 10.0f},
              Color{0, 0, 0, 150});
    DrawLineV(Vector2{p.x - 9.0f, p.y}, Vector2{p.x + 9.0f, p.y},
              Color{255, 255, 255, 230});
    DrawLineV(Vector2{p.x, p.y - 9.0f}, Vector2{p.x, p.y + 9.0f},
              Color{255, 255, 255, 230});
}

std::string ProbeTitle(int i, int j) {
    const float x = ((static_cast<float>(i) + 0.5f) / kSize - 0.5f) * kLength;
    const float y = (0.5f - (static_cast<float>(j) + 0.5f) / kSize) * kLength;
    return std::format("PROBE   x {:+.2f} m   y {:+.2f} m", x, y);
}

void DrawScope(const ProbeHistory& history, const Membrane& model, int i, int j,
               float& scale) {
    const Layout layout = MakeLayout(GetScreenWidth(), GetScreenHeight());
    const Rectangle panel = layout.scope;
    const Rectangle plot{panel.x + 18.0f, panel.y + 44.0f, panel.width - 36.0f,
                         panel.height - 76.0f};
    const float mid = plot.y + plot.height * 0.5f;
    const int count = std::min(history.Size(), kHistoryWindowSamples);
    const float value = model.At(i, j);

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

    const std::string title = ProbeTitle(i, j);
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
        const Color trace_color{190, 220, 255, 235};
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
        Membrane model;
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
        check(model.MaxSymmetryError() < 1e-6f, "x/y symmetry",
              std::format("max error {:.2e} m", model.MaxSymmetryError()));
        check(drift < 0.02f, "undamped energy drift",
              std::format("relative drift {:.4f}", drift));
    }

    {
        Membrane model;
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
    InitWindow(1280, 900, "Square Membrane");
    SetWindowMinSize(960, 720);
    SetTargetFPS(60);

    Membrane model;
    ProbeHistory history;
    int probe_i = (kSize - 1) / 2;
    int probe_j = (kSize - 1) / 2;
    float wave_scale = 0.01f;
    float color_scale = 0.001f;
    int sample_countdown = kStepsPerSample;
    float accumulator = 0.0f;
    bool probe_dragging = false;
    float speed = kSpeedDefault;
    bool speed_dragging = false;

    Image image = GenImageColor(kSize, kSize, BLANK);
    Texture2D texture = LoadTextureFromImage(image);
    SetTextureFilter(texture, TEXTURE_FILTER_BILINEAR);

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
            color_scale = 0.001f;
            sample_countdown = kStepsPerSample;
            accumulator = 0.0f;
        }
        if (IsKeyPressed(KEY_E)) {
            std::cout << std::format("wav export: rendering {:.0f} s...\n",
                                     kExportSeconds)
                      << std::flush;
            Membrane export_model;
            export_model.SetDamping(kExportDamping);
            export_model.ResetImpulse();
            const bool ok = ExportProbeWav(
                "vibrate_membrane.wav", export_model, kExportStepsPerSecond2D,
                kExportSeconds,
                [](Membrane& m, float dt) { m.Step(dt); },
                [&](const Membrane& m) { return m.At(probe_i, probe_j); });
            std::cout << std::format("wav export: {}\n",
                                     ok ? "vibrate_membrane.wav" : "failed");
        }

        const Rectangle speed_hit = SpeedSliderHit(width);
        const Rectangle speed_track = SpeedSliderTrack(width);
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) &&
            CheckCollisionPointRec(mouse, speed_hit)) {
            speed_dragging = true;
        } else if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) &&
                   CheckCollisionPointRec(mouse, layout.membrane)) {
            probe_dragging = true;
        }
        if (speed_dragging && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
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
            }
        }
        if (probe_dragging && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            const float nx = std::clamp(
                (mouse.x - layout.membrane.x) / layout.membrane.width, 0.0f, 1.0f);
            const float ny = std::clamp(
                (mouse.y - layout.membrane.y) / layout.membrane.height, 0.0f, 1.0f);
            const int i = std::clamp(static_cast<int>(nx * kSize), 1, kSize - 2);
            const int j = std::clamp(static_cast<int>(ny * kSize), 1, kSize - 2);
            if (i != probe_i || j != probe_j) {
                probe_i = i;
                probe_j = j;
                history.Clear();
                wave_scale = 0.01f;
                sample_countdown = kStepsPerSample;
            }
        }
        if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            probe_dragging = false;
            speed_dragging = false;
        }

        accumulator += std::min(GetFrameTime(), kMaxFrameTime) * speed;
        int steps = 0;
        while (accumulator >= kFixedDt && steps < kMaxStepsPerFrame) {
            model.Step(kFixedDt);
            if (--sample_countdown <= 0) {
                history.Push(model.At(probe_i, probe_j));
                sample_countdown = kStepsPerSample;
            }
            accumulator -= kFixedDt;
            ++steps;
        }

        float peak = 0.0f;
        for (float value : model.U()) {
            peak = std::max(peak, std::abs(value));
        }
        const float target_scale = std::max(peak * 1.08f, 0.000001f);
        if (target_scale > color_scale) {
            color_scale = target_scale;
        } else {
            color_scale += (target_scale - color_scale) * 0.02f;
        }
        FillColormap(model, color_scale, image);
        UpdateTexture(texture, image.data);

        BeginDrawing();
        ClearBackground(Color{8, 10, 16, 255});
        DrawTexturePro(texture, Rectangle{0.0f, 0.0f, static_cast<float>(kSize),
                                          static_cast<float>(kSize)},
                       layout.membrane, Vector2{0.0f, 0.0f}, 0.0f, WHITE);
        DrawRectangleLinesEx(layout.membrane, 1.0f, Color{52, 68, 90, 255});
        DrawProbeMarker(layout.membrane, probe_i, probe_j);
        DrawSceneHud(model);
        DrawSpeedSlider(width, speed);
        DrawScope(history, model, probe_i, probe_j, wave_scale);
        EndDrawing();
    }

    UnloadTexture(texture);
    UnloadImage(image);
    CloseWindow();
    return 0;
}
