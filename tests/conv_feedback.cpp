#include <qwqdsp/fx/uniform_convolution.hpp>

int main() {
    qwqdsp_fx::UniformConvolution conv;
    conv.Init(32);
    conv.Reset();

    float ir[512]{1.0f};
    conv.SetIR(ir);

    float input[2048]{1.0f};

    float lag = 0.0f;
    for (int i = 0; i < 2048; ++i) {
        float x = input[i];

        float temp[1]{x + lag};
        conv.Process(temp);
        lag = temp[0] * 0.99f;

        input[i] = temp[0];
    }
}
