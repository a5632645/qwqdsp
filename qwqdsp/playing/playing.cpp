#include <raylib.h>
#include <numbers>
#include <algorithm>

#include "qwqdsp/filter/fir_design.hpp"
#include "qwqdsp/window/helper.hpp"
#include "qwqdsp/spectral/real_fft.hpp"
#include "qwqdsp/convert.hpp"

#include "slider.hpp"

int main() {
    InitWindow(800, 600, "fir");
    SetTargetFPS(15);

    Rectangle button{0,0,100,30};
    Knob w;
    w.set_bound(0, 30, 80, 80);
    w.set_range(3.14f / 512, 3.14f / 2, 3.14f / 512, 0.314f);
    w.set_bg_color(Color{255, 251, 232, 0xff});
    w.set_fore_color(ORANGE);
    w.set_title("width");
    Knob k;
    k.set_bound(0, 110, 80, 80);
    k.set_range(1, 1000, 10, 50);
    k.set_bg_color(Color{255, 251, 232, 0xff});
    k.set_fore_color(ORANGE);
    k.set_title("k");

    float h[64]{};
    float pad[2048]{};
    float gain[1025]{};
    float last_gain[1025]{};

    qwqdsp::spectral::RealFFT fft;
    fft.Init(2048);

    while (!WindowShouldClose()) {
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            auto pos = GetMousePosition();
            if (CheckCollisionPointRec(pos, button)) {
                qwqdsp::filter::FirDesign::LowpassLeastSquare(h, std::numbers::pi_v<float> / 2, w.get_value(), k.get_value());
                qwqdsp::window::Helper::ZeroPad(pad, h);
                std::copy(gain, std::end(gain), last_gain);
                fft.FFTGainPhase(pad, gain);
                for (auto& s : gain) {
                    s = qwqdsp::convert::Gain2Db(s);
                    s = std::max(s, -120.0f);
                    s = (s + 120.0f) / 120.0f;
                }
            }
        }

        BeginDrawing();
            ClearBackground(Color{255, 251, 232, 0xff});
            DrawRectangleLines(button.x, button.y, button.width, button.height, BLACK);
            DrawText("draw", button.x, button.y, button.height, ORANGE);

            float lastx = 0;
            float lasty = 0;
            for (int i = 0; i < 800; ++i) {
                int idx = i * 1023.0f / 800.0f;
                auto x = gain[idx];
                auto y = std::lerp(600, 50, std::clamp(x, 0.0f, 1.0f));
                DrawLine(lastx, lasty, i, y, ORANGE);
                lastx = i;
                lasty = y;
            }

            lastx = 0;
            lasty = 0;
            for (int i = 0; i < 800; ++i) {
                int idx = i * 1023.0f / 800.0f;
                auto x = last_gain[idx];
                auto y = std::lerp(600, 50, std::clamp(x, 0.0f, 1.0f));
                DrawLine(lastx, lasty, i, y, BLUE);
                lastx = i;
                lasty = y;
            }

            w.display();
            k.display();
        EndDrawing();
    }
    CloseWindow();
}