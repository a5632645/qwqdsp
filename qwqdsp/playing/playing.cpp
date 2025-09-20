#include "qwqdsp/filter/parallel_allpass.hpp"
#include "qwqdsp/convert.hpp"


int main() {
    qwqdsp::filter::ParallelAllpass apf;
    qwqdsp::filter::ParallelAllpass apf2;
    apf.BuildButterworth(9, std::numbers::pi_v<float> / 4);
    apf2.BuildButterworth(9, std::numbers::pi_v<float> / 4);

    float gains[1024];
    for (size_t i = 0; i < 1024; ++i) {
        float const w = std::numbers::pi_v<float> * i / 1024.0f;
        auto[up, down] = apf.GetResponce(std::polar(1.0f, w));
        up *= up;
        down *= down;
        up *= up;
        down *= down;

        gains[i] = std::abs(up + down);
        gains[i] = qwqdsp::convert::Gain2Db(gains[i]);
        gains[i] = std::max(-80.0f, gains[i]);
    }
}