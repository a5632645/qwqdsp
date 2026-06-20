#include "AudioFile.h"
#include "rnnoise.h"

int main() {
    AudioFile<float> file;
    if (file.load(R"(C:\Users\Kawai\Desktop\noise.wav)")) {
        auto& io = file.samples.front();

        std::vector<int16_t> int16;
        int16.reserve(io.size());
        for (auto& s : io) {
            int16.push_back(s * INT16_MAX);
        }

        auto int16bck = int16;

        RnnoiseState s;
        RnnoiseState_Init(&s);
        RnnoiseState_Reset(&s);
        RnnoiseState_Process(&s, int16bck.data(), int16.data(), int16.size());

        for (size_t i = 0; i < io.size(); ++i) {
            io[i] = int16[i] / (float)INT16_MAX;
        }

        file.setNumChannels(1);
        file.save(R"(C:\Users\Kawai\Desktop\noise_test.wav)");
    }
}
