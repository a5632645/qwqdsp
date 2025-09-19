#include <vector>

#include "qwqdsp/filter/window_fir.hpp"
#include "qwqdsp/window/kaiser.hpp"

int main() {
    size_t num_subspan = 4;
    size_t n = 3;
    size_t size = n + (n - 1) * num_subspan;
    size_t toal = n * (num_subspan + 2);

    std::vector<float> kenerl;
    kenerl.resize(size);
    qwqdsp::filter::WindowFIR::Lowpass(kenerl, std::numbers::pi_v<float> / 4 / (num_subspan + 1));
    qwqdsp::window::Kaiser::ApplyWindow(kenerl, qwqdsp::window::Kaiser::Beta(60), false);

    std::vector<float> lut;
    lut.resize(toal);
    for (size_t i = 0; i < num_subspan + 1; ++i) {
        size_t const aa = i == 0 ? n : n - 1;
        for (size_t j = 0; j < aa; ++j) {
            lut[i * n + j] = kenerl[i + (num_subspan + 1) * j];
        }
    }
    for (size_t i = 0; i < n - 1; ++i) {
        lut[(num_subspan + 1) * n + i] = lut[i + 1];
    }
}
