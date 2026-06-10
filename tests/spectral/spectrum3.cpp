#include <array>

#include <qwqdsp/spectral/real_fft.hpp>
#include <qwqdsp/window/hann.hpp>
#include <qwqdsp/window/helper.hpp>

#include "../../playing/miniaudio.h"
#include "raylib.h"

// ── 窗口与 UI 常量 ──
static constexpr int kWindowWidth  = 640;
static constexpr int kWindowHeight = 320;

// ── 画布布局 ──
static constexpr int kCanvasX = 50;
static constexpr int kCanvasY = 30;
static constexpr int kCanvasW = kWindowWidth - kCanvasX - 20;
static constexpr int kCanvasH = kWindowHeight - kCanvasY - 60;

// ── 网格与刻度 ──
static constexpr float kDbFloor       = -80.0f;
static constexpr float kYTicks[]      = {0.0f, -20.0f, -40.0f, -60.0f, -80.0f};
static constexpr float kXTickLabels[] = {20, 200, 2000, 20000};
static const Color     kGridColor     = {50, 50, 50, 255};
static const Color     kTextColor     = {180, 180, 180, 255};
static const Color     kBgColor       = {20, 20, 20, 255};

// ── Lanczos ──
static constexpr int   kLanczosA = 2;  // 6-tap
static constexpr float kPi       = 3.14159265358979323846f;

// ── 曲线与散点 ──
static const Color     kCurveColor = {120, 200, 255, 255};
static constexpr float kScatterSize = 3.0f;

// ----------------------------------------
// fft
// ----------------------------------------

constexpr int kSampleRate = 48000;
constexpr int kFftSize = 4096;
constexpr int kBinSize = kFftSize / 2 + 1;
constexpr int kHopSize = kFftSize / 4;

static std::array<float, kFftSize> hann_window_{};
static std::array<float, kFftSize> in_buffer_{};
static std::array<float, kFftSize> fft_in_buffer_{};
static qwqdsp_spectral::RealFFT fft_;
static int in_count_{};

// ----------------------------------------
// octave
// ----------------------------------------

constexpr float kRefFreq = 1000;
constexpr float kCenterFreqFloor = kRefFreq / 64;
constexpr int   kOctave     = 12;
constexpr float kAttackMs   = 1.0f;
constexpr float kReleaseMs  = 150.0f;

struct OctaveSmoothIndex {
  int begin_idx;
  int end_idx;
  float center_freq;
};
std::vector<OctaveSmoothIndex> octave_idx_;
std::vector<float> octave_avg_mul_;
std::vector<float> octave_tilt_mul_;
std::vector<float> octave_gain_;
std::vector<float> octave_log_freq_;
std::vector<float> octave_gain_dB_;
static std::array<float, kBinSize> fft_envelope_dB_{};
static float kAttackAlpha  = 0.0f;
static float kReleaseAlpha = 0.0f;

static void Octave_Init() {
  float fcenter = kCenterFreqFloor;
  float fceil_mul = std::exp2(1.0f / (2.0f * kOctave));
  int i = 1;
  do {
    fcenter = kCenterFreqFloor *
              std::exp2(static_cast<float>(i++) / static_cast<float>(kOctave));
    float fceil = fcenter * fceil_mul;
    float ffloor = fcenter / fceil_mul;

    float ceil_bin = fceil / kSampleRate * kFftSize;
    float floor_bin = ffloor / kSampleRate * kFftSize;

    int i32_ceil_bin = std::round(ceil_bin);
    int i32_floor_bin = std::round(floor_bin);
    if (i32_ceil_bin == i32_floor_bin) {
      // ++i32_ceil_bin;
      continue;
    }

    i32_ceil_bin = std::min(i32_ceil_bin, kBinSize - 1);
    i32_floor_bin = std::min(i32_floor_bin, kBinSize - 1);

    octave_idx_.push_back(
        OctaveSmoothIndex{i32_floor_bin, i32_ceil_bin, fcenter});
    // octave_avg_mul_.push_back(1.0f / (ceil_bin - floor_bin));
    octave_avg_mul_.push_back(1.0f / (i32_ceil_bin - i32_floor_bin));

    float tilt_dB = 4.5f * std::log2(fcenter / kRefFreq);
    octave_tilt_mul_.push_back(std::pow(10.0f, tilt_dB / 20.0f));
  } while (fcenter < 20000);

  octave_gain_.resize(octave_idx_.size());
  octave_gain_dB_.resize(octave_idx_.size());
  octave_log_freq_.resize(octave_idx_.size());
  for (size_t j = 0; j < octave_idx_.size(); ++j) {
    octave_log_freq_[j] = std::log10(octave_idx_[j].center_freq);
  }
}

// ----------------------------------------
// function
// ----------------------------------------

static void Dsp_Init() {
  fft_.Init(kFftSize);
  qwqdsp_window::Hann::Window(hann_window_, false);
  qwqdsp_window::Helper::Normalize(hann_window_);

  // 预计算一阶平滑系数
  float T = static_cast<float>(kHopSize) / static_cast<float>(kSampleRate);
  kAttackAlpha  = 1.0f - std::exp(-T / (kAttackMs  * 0.001f));
  kReleaseAlpha = 1.0f - std::exp(-T / (kReleaseMs * 0.001f));

  fft_envelope_dB_.fill(kDbFloor);
}

static void Dsp_Process() {
  std::array<float, kBinSize> gain_buffer{};
  fft_.FFTGainPhase(fft_in_buffer_, gain_buffer);

  // 功率前缀和 (平方)
  std::array<float, kBinSize + 1> power_cum{};
  for (size_t j = 0; j < kBinSize; ++j) {
    float m = gain_buffer[j];
    power_cum[j + 1] = power_cum[j] + m * m;
  }

  int octave_count = octave_idx_.size();
  for (int i = 0; i < octave_count; ++i) {
    auto &d = octave_idx_[i];
    float sum_power = power_cum[d.end_idx + 1] - power_cum[d.begin_idx];
    float n = static_cast<float>(d.end_idx + 1 - d.begin_idx);
    float raw    = std::sqrt(sum_power / n) * octave_tilt_mul_[i];
    float raw_dB = 20.0f * log10f(raw + 1e-12f);
    float prev   = octave_gain_[i];  // dB
    float alpha  = (raw_dB > prev) ? kAttackAlpha : kReleaseAlpha;
    octave_gain_[i]    = prev + alpha * (raw_dB - prev);
    octave_gain_dB_[i] = octave_gain_[i];
  }

  // ── FFT bin 包络 (一阶平滑) ──
  for (int j = 0; j < kBinSize; ++j) {
    float dB  = 20.0f * log10f(gain_buffer[j] + 1e-12f);
    float prev = fft_envelope_dB_[j];
    float alpha = (dB > prev) ? kAttackAlpha : kReleaseAlpha;
    fft_envelope_dB_[j] = prev + alpha * (dB - prev);
  }
}

extern "C" void MaCaptureCallback(ma_device *pDevice, void *pOutput,
                                  const void *pInput, ma_uint32 frameCount) {
  float const *src = reinterpret_cast<float const *>(pInput);
  float* dst = reinterpret_cast<float*>(pOutput);
  int num_frame = frameCount;
  std::copy_n(src, num_frame, dst);

  while (num_frame != 0) {
    int need = kFftSize - in_count_;
    int can_do = std::min(need, num_frame);
    std::copy_n(src, can_do, in_buffer_.begin() + in_count_);

    in_count_ += can_do;
    num_frame -= can_do;
    src += can_do;

    if (in_count_ == kFftSize) {
      for (int i = 0; i < kFftSize; ++i) {
        fft_in_buffer_[i] = hann_window_[i] * in_buffer_[i];
      }
      for (int i = 0; i < (kFftSize - kHopSize); ++i) {
        in_buffer_[i] = in_buffer_[i + kHopSize];
      }
      in_count_ -= kHopSize;
      Dsp_Process();
    }
  }
}

// ── Lanczos 查找表 ──
struct LanczosTap {
  int   idx;
  float weight;
};
static std::vector<std::vector<LanczosTap>> lanczos_table_;

static void InitLanczosTable() {
  const int octaveCnt = static_cast<int>(octave_idx_.size());
  const float logMin = std::log10(octave_idx_[0].center_freq);
  const float logMax = std::log10(octave_idx_[octaveCnt - 1].center_freq);
  lanczos_table_.resize(kCanvasW);
  for (int px = 0; px < kCanvasW; ++px) {
    float normPx = static_cast<float>(px) / static_cast<float>(kCanvasW - 1);
    float logF  = logMin + normPx * (logMax - logMin);
    float center = (logF - logMin) / (logMax - logMin) * static_cast<float>(octaveCnt - 1);
    int start = std::max(0, static_cast<int>(std::ceil(center - kLanczosA)));
    int end   = std::min(octaveCnt - 1,
                         static_cast<int>(std::floor(center + kLanczosA)));
    float wsum = 0.0f;
    for (int j = start; j <= end; ++j) {
      float x = center - static_cast<float>(j);
      float w;
      if (x == 0.0f) {
        w = 1.0f;
      } else {
        float pix = kPi * x;
        w = static_cast<float>(kLanczosA) * std::sin(pix) *
            std::sin(pix / static_cast<float>(kLanczosA)) / (pix * pix);
      }
      lanczos_table_[px].push_back({j, w});
      wsum += w;
    }
    if (wsum > 0.0f) {
      for (auto &tap : lanczos_table_[px]) {
        tap.weight /= wsum;
      }
    }
  }
}

// ── 抛物线插值表 (非等距) ──
struct ParaTap {
  int   idx_l, idx_m, idx_r;
  float w_l, w_m, w_r;
};
static std::vector<ParaTap> para_table_;

static void InitParaTable() {
  const int octaveCnt = static_cast<int>(octave_idx_.size());
  const float logMin = octave_log_freq_[0];
  const float logMax = octave_log_freq_[octaveCnt - 1];
  para_table_.resize(kCanvasW);

  for (int px = 0; px < kCanvasW; ++px) {
    float normPx = static_cast<float>(px) / static_cast<float>(kCanvasW - 1);
    float logF  = logMin + normPx * (logMax - logMin);

    // 找到区间 [seg, seg+1]
    int seg = 0;
    if (logF >= octave_log_freq_[octaveCnt - 1]) {
      seg = octaveCnt - 2;
    } else if (logF <= octave_log_freq_[0]) {
      seg = 0;
    } else {
      for (int j = 0; j < octaveCnt - 1; ++j) {
        if (logF >= octave_log_freq_[j] && logF < octave_log_freq_[j + 1]) { seg = j; break; }
      }
    }

    // 对区间 [seg, seg+1] 使用 3 个 band (seg-1, seg, seg+1)
    int left  = std::max(0, seg - 1);
    int mid   = seg;
    int right = std::min(octaveCnt - 1, seg + 1);
    if (left == mid) right = std::min(octaveCnt - 1, mid + 2);
    if (right == mid) left  = std::max(0, mid - 2);
    mid = left + 1;

    float x0 = octave_log_freq_[left];
    float x1 = octave_log_freq_[mid];
    float x2 = octave_log_freq_[right];

    // 二次 Lagrange 基函数
    float w0 = (logF - x1) * (logF - x2) / ((x0 - x1) * (x0 - x2));
    float w1 = (logF - x0) * (logF - x2) / ((x1 - x0) * (x1 - x2));
    float w2 = (logF - x0) * (logF - x1) / ((x2 - x0) * (x2 - x1));

    para_table_[px] = {left, mid, right, w0, w1, w2};
  }
}

// ── 绘制函数 ──

static void DrawBackground() {
  DrawRectangle(kCanvasX, kCanvasY, kCanvasW, kCanvasH, kBgColor);
}

static void DrawYAxisGrid() {
  for (float db : kYTicks) {
    float norm = (db - kDbFloor) / (-kDbFloor);
    int y = kCanvasY + kCanvasH - static_cast<int>(norm * kCanvasH);
    DrawLine(kCanvasX, y, kCanvasX + kCanvasW, y, kGridColor);
    char label[16];
    snprintf(label, sizeof(label), "%.0f", db);
    int tw = MeasureText(label, 10);
    DrawText(label, kCanvasX - tw - 6, y - 5, 10, kTextColor);
  }
}

static void DrawXAxisGrid(int octaveCnt) {
  const float logMin = octave_log_freq_[0];
  const float logMax = octave_log_freq_[octaveCnt - 1];

  // 三个频率 decade: 20→200, 200→2k, 2k→20k，每个等分 10 格，按 log 定位
  constexpr float kDecadeStart[] = {20, 200, 2000};
  for (float start : kDecadeStart) {
    float step = start;
    for (float freq = start; freq <= start * 10.0f; freq += step) {
      float norm = (std::log10(freq) - logMin) / (logMax - logMin);
      float sx = kCanvasX + norm * kCanvasW;
      DrawLine(static_cast<int>(sx), kCanvasY, static_cast<int>(sx),
               kCanvasY + kCanvasH, kGridColor);
    }
  }

  // 标签仅 20 / 200 / 2k / 20k
  for (float freq : kXTickLabels) {
    float norm = (std::log10(freq) - logMin) / (logMax - logMin);
    float sx = kCanvasX + norm * kCanvasW;
    char label[16];
    if (freq >= 1000.0f)
      snprintf(label, sizeof(label), "%.0fk", freq / 1000.0f);
    else
      snprintf(label, sizeof(label), "%.0f", freq);
    int tw = MeasureText(label, 10);
    DrawText(label, static_cast<int>(sx) - tw / 2, kCanvasY + kCanvasH + 4,
             10, kTextColor);
  }
}

static void DrawFFTEnvelope(int octaveCnt) {
  const float logMin = octave_log_freq_[0];
  const float logMax = octave_log_freq_[octaveCnt - 1];
  const float binToFreq = static_cast<float>(kSampleRate) / static_cast<float>(kFftSize);

  // 收集可见范围内的 bin 点
  std::vector<Vector2> pts;
  pts.reserve(kBinSize);
  for (int j = 0; j < kBinSize; ++j) {
    float freqHz = static_cast<float>(j) * binToFreq;
    float logF = std::log10(freqHz);
    if (logF < logMin || logF > logMax) continue;

    float normX = (logF - logMin) / (logMax - logMin);
    float sx = kCanvasX + normX * kCanvasW;
    float dbVal = fft_envelope_dB_[j];
    float normY = (dbVal - kDbFloor) / (-kDbFloor);
    normY = std::clamp(normY, 0.0f, 1.0f);
    float sy = kCanvasY + kCanvasH - normY * kCanvasH;
    pts.push_back({sx, sy});
  }
  if (pts.size() > 1)
    DrawLineStrip(pts.data(), static_cast<int>(pts.size()),
                  Color{80, 220, 120, 180});  // 浅绿半透明
}

static void DrawOctaveCurve(int /*octaveCnt*/) {
  std::vector<Vector2> curvePoints(kCanvasW);
  for (int px = 0; px < kCanvasW; ++px) {
    float dbSum = 0.0f;
    for (const auto &tap : lanczos_table_[px]) {
      dbSum += tap.weight * (20.0f * log10f(octave_gain_[tap.idx] + 1e-12f));
    }
    float dbVal = std::max(dbSum, kDbFloor);
    float norm  = (dbVal - kDbFloor) / (-kDbFloor);
    float y = kCanvasY + kCanvasH - norm * kCanvasH;
    curvePoints[px] = {static_cast<float>(kCanvasX + px), y};
  }
  DrawLineStrip(curvePoints.data(), kCanvasW, kCurveColor);
}

static void DrawOctaveCurve2(int /*octaveCnt*/) {
  std::vector<Vector2> curvePoints(kCanvasW);
  for (int px = 0; px < kCanvasW; ++px) {
    const auto &tap = para_table_[px];
    float dbVal = std::max(
        tap.w_l * octave_gain_dB_[tap.idx_l] +
        tap.w_m * octave_gain_dB_[tap.idx_m] +
        tap.w_r * octave_gain_dB_[tap.idx_r], kDbFloor);
    float norm  = (dbVal - kDbFloor) / (-kDbFloor);
    float y = kCanvasY + kCanvasH - norm * kCanvasH;
    curvePoints[px] = {static_cast<float>(kCanvasX + px), y};
  }
  DrawLineStrip(curvePoints.data(), kCanvasW, kCurveColor);
}

static void DrawOctaveScatter(int octaveCnt) {
  const float logMin = octave_log_freq_[0];
  const float logMax = octave_log_freq_[octaveCnt - 1];
  for (int i = 0; i < octaveCnt; ++i) {
    float norm = (octave_gain_dB_[i] - kDbFloor) / (-kDbFloor);
    norm       = std::clamp(norm, 0.0f, 1.0f);
    float sx = kCanvasX + (octave_log_freq_[i] - logMin) /
                              (logMax - logMin) * kCanvasW;
    float sy = kCanvasY + kCanvasH - norm * kCanvasH;
    DrawCircleV({sx, sy}, kScatterSize, YELLOW);
  }
}

int main(void) {
  SetConfigFlags(FLAG_MSAA_4X_HINT);
  InitWindow(kWindowWidth, kWindowHeight,
             "Spectrogram - miniaudio + qwqdsp + raylib");
  SetTargetFPS(60);

  // ── miniaudio 全双工 (捕获 + 播放) ──
  ma_device_config config = ma_device_config_init(ma_device_type_duplex);
  config.capture.format = ma_format_f32;
  config.capture.channels = 1;
  config.playback.format = ma_format_f32;
  config.playback.channels = 1;
  config.sampleRate = static_cast<ma_uint32>(kSampleRate);
  config.dataCallback = MaCaptureCallback;
  config.pUserData = nullptr;
  config.periodSizeInMilliseconds = 10;

  ma_device device;
  ma_result result = ma_device_init(nullptr, &config, &device);
  if (result == MA_SUCCESS) {
    ma_device_start(&device);
  } else {
    TraceLog(LOG_WARNING, "miniaudio 捕获设备初始化失败，以静默模式运行");
  }

  // ── 初始化 DSP ──
  Octave_Init();
  Dsp_Init();
  InitLanczosTable();
  InitParaTable();

  // ── 主循环 ──
  while (!WindowShouldClose()) {
    BeginDrawing();
    ClearBackground(BLACK);

    DrawBackground();
    DrawYAxisGrid();

    const int octaveCnt = static_cast<int>(octave_idx_.size());
    if (octaveCnt > 0) {
      DrawFFTEnvelope(octaveCnt);
      DrawXAxisGrid(octaveCnt);
      // DrawOctaveCurve(octaveCnt);
      DrawOctaveCurve2(octaveCnt);
      // DrawOctaveScatter(octaveCnt);
    }

    DrawFPS(kWindowWidth - 80, 10);
    EndDrawing();
  }

  // ── 清理 ──
  if (result == MA_SUCCESS) {
    ma_device_stop(&device);
    ma_device_uninit(&device);
  }
  CloseWindow();
  return 0;
}
