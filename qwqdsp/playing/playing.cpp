#include "qwqdsp/filter/window_fir.hpp"

int main() {
    float test[7];
    qwqdsp::filter::WindowFIR::Lowpass(test, std::numbers::pi_v<float> / 2);
}