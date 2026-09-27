#include "work_dir.hpp"
#include <AudioFile.h>
#include <cmath>
#include <cstdio>
#include <format>
#include <iostream>
#include <qwqdsp/filter/any_hilbert.hpp>
#include <vector>

// ===== 主程序 =====
int main() {
    using namespace qwqdsp_filter;

    constexpr int64_t N = 1 << 18; // 262144 样本
    constexpr int kSampleRate = 48000;

    AnyHilbert<float> hilbert;
    hilbert.Reset();

    // 生成脉冲响应（立体声: ch0=实部, ch1=虚部）
    AudioFile<float>::AudioBuffer buf(2);
    buf[0].resize(N);
    buf[1].resize(N);

    // 第一个样本: 单位冲激
    auto out = hilbert.Tick(1.0f);
    buf[0][0] = out.real();
    buf[1][0] = out.imag();

    // 后续样本: 输入为 0
    for (int64_t n = 1; n < N; ++n) {
        out = hilbert.Tick(0.0f);
        buf[0][n] = out.real();
        buf[1][n] = out.imag();
    }

    // 导出为 WAV
    AudioFile<float> file;
    file.setAudioBuffer(buf);
    file.setSampleRate(kSampleRate);
    file.setBitDepth(32);

    auto outPath = qwqdsp_support::OutputFile("cheby_hilbert_impulse.wav");
    if (file.save(outPath)) {
        std::cout << std::format("已导出: {} ({} 样本, 立体声 float32, {} Hz)\n", outPath, N, kSampleRate);
        std::cout << "  左声道 = 实部, 右声道 = 虚部\n";
    }
    else {
        std::cerr << std::format("写入失败: {}\n", outPath);
        return 1;
    }

    return 0;
}
