#!/usr/bin/env python3
"""
PSD Spectrum — 实时音频 25×7 LED 频谱显示

从麦克风采集 16-bit PCM → 加汉宁窗 → FFT → 峰值 PSD → 25 个 Bin
→ 与 7 级滑条阈值比较 → 点亮对应 LED

所有处理均在 16-bit 整数域模拟，与 C++ 算法行为一致。
"""

import sys
import math
from typing import Optional

import numpy as np
import sounddevice as sd
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QSlider, QLabel, QFrame, QGridLayout, QSizePolicy, QStatusBar,
    QPushButton, QGroupBox, QSpinBox, QComboBox,
)
from PyQt6.QtCore import Qt, QTimer, pyqtSignal, pyqtSlot
from PyQt6.QtGui import QPainter, QPainterPath, QColor, QBrush, QPen, QFont, QLinearGradient

# ════════════════════════════════════════════════════════════════
#  常量 (与 spectrum_driver.h / .cpp 一致)
# ════════════════════════════════════════════════════════════════
SPECTRUM_BIN_SIZE = 25
MAX_AMP = 32768
SAMPLE_RATE = 44100

LED_ROW_COUNT = 7

# 基础频率步进表 (用于 SPECTRUM_SIZE=256)
# 实际使用时按 kMul = size/256 缩放
_BASE_STEP_TABLE = [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5]


def build_step_table(fft_size: int) -> list:
    """根据 FFT 大小生成对应的频率步进表"""
    k = max(1, fft_size // 256)
    return [v * k for v in _BASE_STEP_TABLE]


def compute_bin_frequencies(sample_rate: float, fft_size: int,
                            step_table: list) -> list:
    """计算每个 bin 的中心频率 (Hz), 对应 C++ bin_frq"""
    frq_one_unit = sample_rate / fft_size
    freqs = []
    s = 0
    for st in step_table:
        freqs.append(frq_one_unit * (st + s))
        s += st
    return freqs


def fmt_freq(hz: float) -> str:
    """频率格式: ≥1000 Hz → x.xk, 否则整数 Hz"""
    if hz >= 1000:
        return f"{hz / 1000:.1f}k"
    return f"{int(hz)}"


# ════════════════════════════════════════════════════════════════
#  算法初始化 — 窗函数 + 增益
# ════════════════════════════════════════════════════════════════

WINDOW_HANN     = "Hann"
WINDOW_HAMMING  = "Hamming"
WINDOW_BLACKMAN = "Blackman"
WINDOW_BHARRIS  = "Blackman-Harris"
WINDOW_FLATTOP  = "Flat top"
WINDOW_RECT     = "Rectangle"

WINDOW_NAMES = [WINDOW_HANN, WINDOW_HAMMING, WINDOW_BLACKMAN,
                WINDOW_BHARRIS, WINDOW_FLATTOP, WINDOW_RECT]


def init_window(size: int, window_type: str = WINDOW_HANN):
    """生成指定窗系数 (int32) 与增益

    返回值 (hm, window_gain):
        hm[i]          = int(窗函数值 × 32767)
        window_gain    = sum(窗函数值) / 2
    """
    hm = [0] * size
    fsum = 0.0
    twopi = 2.0 * math.pi
    N = size
    for i in range(size):
        if window_type == WINDOW_HANN:
            w = 0.5 - 0.5 * math.cos(twopi * i / (N - 1))
        elif window_type == WINDOW_HAMMING:
            w = 0.54 - 0.46 * math.cos(twopi * i / (N - 1))
        elif window_type == WINDOW_BLACKMAN:
            w = (0.42 - 0.5 * math.cos(twopi * i / (N - 1))
                 + 0.08 * math.cos(2 * twopi * i / (N - 1)))
        elif window_type == WINDOW_BHARRIS:
            w = (0.35875 - 0.48829 * math.cos(twopi * i / (N - 1))
                 + 0.14128 * math.cos(2 * twopi * i / (N - 1))
                 - 0.01168 * math.cos(3 * twopi * i / (N - 1)))
        elif window_type == WINDOW_FLATTOP:
            w = (0.2156 - 0.4160 * math.cos(twopi * i / (N - 1))
                 + 0.2781 * math.cos(2 * twopi * i / (N - 1))
                 - 0.0836 * math.cos(3 * twopi * i / (N - 1))
                 + 0.0069 * math.cos(4 * twopi * i / (N - 1)))
        else:  # Rectangle
            w = 1.0
        hm[i] = int(w * 32767)
        fsum += w
    window_gain = fsum / 2.0 if fsum > 0 else 1.0
    return hm, window_gain


# ════════════════════════════════════════════════════════════════
#  算法核心 — 处理一帧 int16 PCM, 返回 25-bin PSD (0~MAX_PSD)
# ════════════════════════════════════════════════════════════════
def compute_bin_psd(pcm_int16: np.ndarray, hm: list, wgain: float,
                    step_table: list) -> np.ndarray:
    """全程模拟 16-bit 整数计算: 加窗 → FFT → 峰值 PSD

    不做任何 [-1,1] 归一化。模拟 C++ int32 定点运算:
        fft_buf[i] = (pcm[i] * hm[i] / window_gain) >> 15
    """
    n = len(pcm_int16)

    # ── 1. 加窗 (纯 int32, wgain 移至 FFT 后) ──
    fft_in = np.zeros(n, dtype=np.float64)
    for i in range(n):
        product = int(pcm_int16[i]) * hm[i]        # int32 × int32
        val = int(product) >> 15                    # 仅右移, 不含 wgain
        fft_in[i] = float(val)

    # ── 2. 实数 FFT ──
    spectrum = np.fft.rfft(fft_in, n=n)

    # ── 3. 原始 PSD (未归一化): re² + im² ──
    re_raw = np.real(spectrum)
    im_raw = np.imag(spectrum)
    psd_raw = re_raw * re_raw + im_raw * im_raw

    # ── 4. 溢出检测 ──
    INT64_MAX = 9_223_372_036_854_775_807
    if np.any(psd_raw > INT64_MAX):
        raise OverflowError(
            f"PSD 值 ({psd_raw.max():.0f}) 超过 int64 最大值 ({INT64_MAX})"
        )

    # ── 5. wgain² 归一化 ──
    w2 = wgain * wgain
    psd = psd_raw / w2

    # ── 6. 置零 DC / Nyquist (与 C++ 一致) ──
    psd[0] = 0.0
    if len(psd) > 1:
        psd[-1] = 0.0

    # ── 7. 按步进表分组, 取组内峰值 ──
    bin_out = np.zeros(SPECTRUM_BIN_SIZE, dtype=np.float64)
    j = k = 0
    for i in range(len(psd)):
        if k >= SPECTRUM_BIN_SIZE:
            break
        if psd[i] > bin_out[k]:
            bin_out[k] = psd[i]
        j += 1
        if j >= step_table[k]:
            j -= step_table[k]
            k += 1

    # ── 8. 原始 FFT 幅度 (dB), 供折线频谱使用 ──
    mag = np.sqrt(psd_raw) / wgain  # 等效于 sqrt(psd_raw / w2)
    mag_db = np.full_like(mag, -100.0)
    ok = mag > 0
    mag_db[ok] = 20.0 * np.log10(mag[ok] / MAX_AMP)
    mag_db[0] = -100.0          # 移除 DC
    if len(mag_db) > 1:
        mag_db[-1] = -100.0     # 移除 Nyquist

    return bin_out, mag_db


def db_to_power(db: float) -> float:
    """幅值 dB → PSD 阈值: (32768 × 10^(dB/20))²"""
    amp = MAX_AMP * (10.0 ** (db / 20.0))
    return amp * amp


def power_to_db(power: float) -> float:
    """功率谱 → 幅值 dB: 20·log10(√power / 32768), 钳位 [-100, 0]"""
    if power <= 0.0:
        return -100.0
    amp = math.sqrt(power)
    db = 20.0 * math.log10(amp / MAX_AMP)
    return max(-100.0, min(0.0, db))


def psd_to_level(bin_psd: np.ndarray, thresholds_db: list) -> np.ndarray:
    """根据 7 个滑条的 dB 阈值, 返回 (7, 25) bool 矩阵

    led_lit[row][col] = True 表示该 LED 亮起:
        row 0 = Lv.7 (最强), row 6 = Lv.1 (最弱)
    每个 bin 的 PSD 与 7 个阈值逐一比较, 超过则亮。
    """
    lit = np.zeros((LED_ROW_COUNT, SPECTRUM_BIN_SIZE), dtype=bool)
    for row in range(LED_ROW_COUNT):
        threshold = db_to_power(thresholds_db[row])
        lit[row, :] = bin_psd > threshold
    return lit


def format_power(power: float) -> str:
    if power >= 1e9:
        return f"{power / 1e9:.3f} G"
    if power >= 1e6:
        return f"{power / 1e6:.3f} M"
    if power >= 1e3:
        return f"{power / 1e3:.3f} K"
    return f"{power:.1f}"


# ════════════════════════════════════════════════════════════════
#  LED 指示灯控件
# ════════════════════════════════════════════════════════════════
class LedWidget(QWidget):
    """长方形 LED 灯条 — 间隙小, 排列紧密, 仿均衡器矩阵"""

    def __init__(self, color: QColor, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._on = False
        self._color = color
        self.setFixedSize(32, 14)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)

    def set_on(self, state: bool) -> None:
        if self._on != state:
            self._on = state
            self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        margin = 1  # 仅 1px 暗边 → 间隙极小
        pw = self.width() - margin * 2    # 30
        ph = self.height() - margin * 2   # 12
        radius = 1  # 几乎直角, 更有 LED 灯条感

        if self._on:
            c = self._color
            painter.setPen(QPen(c.lighter(130), 1))
            painter.setBrush(QBrush(c))
            painter.drawRoundedRect(margin, margin, pw, ph, radius, radius)
            # 细长高光 (顶部)
            hl = QColor(255, 255, 255, 45)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QBrush(hl))
            painter.drawRect(margin + 2, margin + 1, pw - 4, ph // 2 - 1)
        else:
            painter.setPen(QPen(QColor(30, 30, 30), 1))
            painter.setBrush(QBrush(QColor(18, 18, 18)))
            painter.drawRoundedRect(margin, margin, pw, ph, radius, radius)


# ════════════════════════════════════════════════════════════════
#  实时频谱线图 (左侧对比)
# ════════════════════════════════════════════════════════════════
class SpectrumLineWidget(QWidget):
    """绘制一条实时的 PSD 折线频谱, 作为 LED 矩阵的直观对比"""

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._mag_db = np.zeros(1)
        self._palette: list[QColor] = []
        self.setMinimumWidth(220)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)

    def update_data(self, mag_db: np.ndarray, palette: list):
        self._mag_db = mag_db
        self._palette = palette
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        if not self._palette:
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        W = self.width()
        H = self.height()
        margin_l, margin_r, margin_t, margin_b = 46, 12, 12, 28
        pw = W - margin_l - margin_r
        ph = H - margin_t - margin_b

        # ── 背景 ──
        painter.fillRect(0, 0, W, H, QColor(12, 12, 12))
        painter.fillRect(margin_l, margin_t, pw, ph, QColor(18, 18, 22))

        mag = self._mag_db
        bin_count = len(mag)

        # ── Y 轴网格 & 刻度 ──
        db_refs = [0, -12, -24, -36, -48, -60, -80]
        painter.setFont(QFont("Consolas", 7))
        for db in db_refs:
            # db → 归一化 y: 0dB 在顶部, -80dB 在底部
            ratio = (db - (-80)) / (0 - (-80))  # 0→1, 0dB→1, -80dB→0
            y = margin_t + int((1.0 - ratio) * ph)
            painter.setPen(QColor(40, 40, 45))
            painter.drawLine(margin_l, y, margin_l + pw, y)
            painter.setPen(QColor(160, 160, 160))
            painter.drawText(2, y + 3, f"{db}")

        # ── X 轴刻度 ──
        step = max(1, bin_count // 10)
        painter.setPen(QColor(160, 160, 160))
        painter.setFont(QFont("Consolas", 6))
        for i in range(0, bin_count, step):
            x = margin_l + int((i + 0.5) / bin_count * pw)
            painter.drawText(x - 6, margin_t + ph + 12, str(i))

        # ── 单位标注 ──
        painter.setPen(QColor(120, 120, 120))
        painter.setFont(QFont("Consolas", 6))
        painter.drawText(margin_l, margin_t + ph + 22, "FFT bin")
        painter.drawText(2, margin_t + 8, "dB")

        # ── 折线 (原始 FFT 幅度 dB) ──
        pts = []
        for i in range(bin_count):
            db = mag[i]
            ratio = (db - (-80)) / (0 - (-80))
            ratio = max(0.0, min(1.0, ratio))
            x = margin_l + int((i + 0.5) / bin_count * pw)
            y = margin_t + int((1.0 - ratio) * ph)
            pts.append((x, y))

        # 填充区域
        if len(pts) > 1:
            path = QPainterPath()
            path.moveTo(pts[0][0], margin_t + ph)
            for x, y in pts:
                path.lineTo(x, y)
            path.lineTo(pts[-1][0], margin_t + ph)
            path.closeSubpath()
            grad = QLinearGradient(0, margin_t, 0, margin_t + ph)
            grad.setColorAt(0.0, QColor(0, 200, 255, 45))
            grad.setColorAt(0.5, QColor(0, 120, 200, 18))
            grad.setColorAt(1.0, QColor(0, 60, 150, 4))
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QBrush(grad))
            painter.drawPath(path)

        # 折线 (青色)
        painter.setPen(QPen(QColor(0, 210, 255), 1.5))
        for i in range(len(pts) - 1):
            painter.drawLine(pts[i][0], pts[i][1],
                             pts[i + 1][0], pts[i + 1][1])

        # 顶点圆点
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QBrush(QColor(0, 230, 255)))
        for x, y in pts:
            painter.drawEllipse(x - 2, y - 2, 4, 4)


class PowerOfTwoSpinBox(QSpinBox):
    """步进按 ×2 / ÷2 的 SpinBox"""

    def stepBy(self, steps: int) -> None:
        v = self.value()
        if steps > 0:
            while v <= self.value() and v < self.maximum():
                v <<= 1
        elif steps < 0:
            while v >= self.value() and v > self.minimum():
                v >>= 1
        self.setValue(max(self.minimum(), min(self.maximum(), v)))


# ════════════════════════════════════════════════════════════════
#  主窗口 — 25×7 LED 矩阵 + 底部滑条
# ════════════════════════════════════════════════════════════════
class MainWindow(QMainWindow):
    """实时频谱 LED 矩阵显示器"""

    # 从音频线程发射新数据
    psd_ready = pyqtSignal(object, object, object)  # bin_psd, mag_db, thresholds_db

    def __init__(self):
        super().__init__()
        self.setWindowTitle("PSD Spectrum — 25×7 LED 实时频谱")
        self.setMinimumSize(940, 600)
        self.resize(940, 620)

        # ── FFT 大小 / 窗类型 ──
        self.fft_size = 256
        self.window_type = WINDOW_HANN
        self.step_table = build_step_table(self.fft_size)

        # ── 算法状态 ──
        self.hm, self.window_gain = init_window(self.fft_size, self.window_type)

        # ── 音频环形暂存 (大小跟随 fft_size) ──
        self._pcm_buf = np.zeros(self.fft_size, dtype=np.int16)
        self._float_ring = np.zeros(self.fft_size, dtype=np.float64)
        self._ring_pos = 0

        # ── 7 个滑条的 dB 阈值 (初始值) ──
        self.threshold_db = [-12, -24, -36, -48, -60, -72, -80]

        # ── 最新 PSD / 原始幅度 / 平滑状态 ──
        self.latest_psd = np.zeros(SPECTRUM_BIN_SIZE)
        self.latest_mag_db = np.zeros(1)
        self._smooth_psd = np.zeros(SPECTRUM_BIN_SIZE)
        self._smooth_mag_db = np.zeros(1)
        self.smooth_tau_ms = 0.0  # 0 = 无平滑

        # ── 音频流 ──
        self.audio_stream: Optional[sd.InputStream] = None

        self._setup_ui()
        self._setup_audio()

        # ── 定时刷新 UI (30 FPS) ──
        self.timer = QTimer(self)
        self.timer.timeout.connect(self._refresh_display)
        self.timer.start(33)

    def _setup_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_vbox = QVBoxLayout(central)
        main_vbox.setSpacing(4)

        # ══════════ 标题 + FFT 尺寸控制 ══════════
        title_row = QHBoxLayout()
        title_row.addStretch()
        title = QLabel("25 频段 × 7 级 LED 实时频谱")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        tf = QFont("Microsoft YaHei", 11)
        tf.setBold(True)
        title.setFont(tf)
        title_row.addWidget(title)
        title_row.addStretch()

        # FFT 大小选择 (×2 / ÷2 步进)
        title_row.addWidget(QLabel("  FFT:"))
        self.fft_spin = PowerOfTwoSpinBox()
        self.fft_spin.setRange(64, 2048)
        self.fft_spin.setValue(self.fft_size)
        self.fft_spin.setSuffix(" pt")
        self.fft_spin.valueChanged.connect(self._on_fft_size_changed)
        title_row.addWidget(self.fft_spin)

        # 窗选择 + 平滑滑条
        title_row.addWidget(QLabel("  窗:"))
        self.window_combo = QComboBox()
        self.window_combo.addItems(WINDOW_NAMES)
        self.window_combo.setCurrentText(self.window_type)
        self.window_combo.currentTextChanged.connect(self._on_window_changed)
        title_row.addWidget(self.window_combo)

        title_row.addWidget(QLabel("  平滑:"))
        self.smooth_slider = QSlider(Qt.Orientation.Horizontal)
        self.smooth_slider.setRange(0, 500)
        self.smooth_slider.setValue(0)
        self.smooth_slider.setFixedWidth(100)
        self.smooth_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.smooth_slider.setTickInterval(100)
        self.smooth_slider.valueChanged.connect(self._on_smooth_changed)
        title_row.addWidget(self.smooth_slider)
        self.smooth_label = QLabel("0 ms")
        self.smooth_label.setFont(QFont("Consolas", 8))
        self.smooth_label.setFixedWidth(40)
        title_row.addWidget(self.smooth_label)

        main_vbox.addLayout(title_row)

        # ══════════ 频谱对比区: 左侧折线 + 右侧 LED 矩阵 ══════════
        spectrum_row = QHBoxLayout()
        spectrum_row.setSpacing(6)

        # ── 左侧: 实时折线频谱 ──
        self.line_spectrum = SpectrumLineWidget()
        spectrum_row.addWidget(self.line_spectrum, stretch=1)

        # ── 右侧: LED 矩阵 (25列 × 7行) ──
        led_group = QGroupBox("LED 矩阵")
        led_group.setSizePolicy(QSizePolicy.Policy.Expanding,
                                QSizePolicy.Policy.Expanding)
        led_group.setFixedWidth(750)
        matrix_layout = QVBoxLayout(led_group)
        matrix_layout.setSpacing(3)

        self.led_grid = QGridLayout()
        self.led_grid.setSpacing(1)

        # ── 生成彩虹调色板 ──
        hue_start, hue_end = 0, 270
        self._palette = []
        for c in range(SPECTRUM_BIN_SIZE):
            hue = int(hue_start + (hue_end - hue_start) * c / (SPECTRUM_BIN_SIZE - 1))
            sat = 200 + 55 * (c % 2)
            self._palette.append(QColor.fromHsv(hue, min(sat, 255), 230))

        self.led_widgets: list[list[LedWidget]] = []
        level_labels = ["Lv.7", "Lv.6", "Lv.5", "Lv.4",
                        "Lv.3", "Lv.2", "Lv.1"]
        for row in range(LED_ROW_COUNT):
            lbl = QLabel(level_labels[row])
            lbl.setFont(QFont("Consolas", 8))
            lbl.setFixedWidth(32)
            self.led_grid.addWidget(lbl, row, 0, Qt.AlignmentFlag.AlignCenter)

            row_leds = []
            for col in range(SPECTRUM_BIN_SIZE):
                led = LedWidget(self._palette[col])
                row_leds.append(led)
                self.led_grid.addWidget(led, row, col + 1,
                                        Qt.AlignmentFlag.AlignCenter)
            self.led_widgets.append(row_leds)

        # 底部频率标签 (黑色)
        bin_freqs = compute_bin_frequencies(SAMPLE_RATE, self.fft_size,
                                            self.step_table)
        for col in range(SPECTRUM_BIN_SIZE):
            lbl = QLabel(fmt_freq(bin_freqs[col]))
            lbl.setFont(QFont("Consolas", 6))
            lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            lbl.setStyleSheet("color: #000000")
            self.led_grid.addWidget(lbl, LED_ROW_COUNT + 1, col + 1)

        matrix_layout.addLayout(self.led_grid)
        spectrum_row.addWidget(led_group, stretch=3)

        main_vbox.addLayout(spectrum_row, stretch=1)

        # ══════════ 分隔 ══════════
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        main_vbox.addWidget(line)

        # ══════════ 7 级阈值滑条 (底部) ══════════
        slider_group = QGroupBox("电平阈值滑条 (dB)")
        slider_layout = QHBoxLayout(slider_group)

        self.sliders: list[QSlider] = []
        self.slider_db_labels: list[QLabel] = []
        self.slider_power_labels: list[QLabel] = []

        default_dbs = [-12, -24, -36, -48, -60, -72, -80]

        for row in range(LED_ROW_COUNT):
            col_layout = QVBoxLayout()

            # 标题: Lv.X
            lbl = QLabel(f"Lv.{LED_ROW_COUNT - row}")
            lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            lbl.setFont(QFont("Consolas", 9, QFont.Weight.Bold))
            col_layout.addWidget(lbl)

            # 垂直滑条
            sl = QSlider(Qt.Orientation.Vertical)
            sl.setRange(-80, 0)
            sl.setValue(default_dbs[row])
            sl.setTickPosition(QSlider.TickPosition.TicksRight)
            sl.setTickInterval(10)
            sl.setFixedHeight(100)
            sl.valueChanged.connect(self._on_threshold_changed)
            col_layout.addWidget(sl, alignment=Qt.AlignmentFlag.AlignCenter)
            self.sliders.append(sl)

            # dB 数值标签
            db_lbl = QLabel(f"{default_dbs[row]} dB")
            db_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            db_lbl.setFont(QFont("Consolas", 8))
            col_layout.addWidget(db_lbl)
            self.slider_db_labels.append(db_lbl)

            # 功率阈值标签
            p_lbl = QLabel(format_power(db_to_power(default_dbs[row])))
            p_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            p_lbl.setFont(QFont("Consolas", 7))
            col_layout.addWidget(p_lbl)
            self.slider_power_labels.append(p_lbl)

            slider_layout.addLayout(col_layout)

        main_vbox.addWidget(slider_group)

        # ══════════ 状态栏 ══════════
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_label = QLabel("就绪")
        self.status_bar.addWidget(self.status_label)

        # ── 信号连接 ──
        self.psd_ready.connect(self._on_psd_ready)

        # ── 初始刷新阈值 → LED ──
        self._enforce_monotonic()
        self._update_labels()

    # ══════════════════════════════════════════════
    #  音频
    # ══════════════════════════════════════════════

    def _setup_audio(self):
        """启动 sounddevice 输入流"""
        try:
            self.audio_stream = sd.InputStream(
                samplerate=SAMPLE_RATE,
                channels=1,
                blocksize=self.fft_size,
                dtype='float32',
                callback=self._audio_callback,
            )
            self.audio_stream.start()
            self.status_label.setText("音频已启动")
        except Exception as e:
            self.status_label.setText(f"音频启动失败: {e}")

    def _audio_callback(self, indata: np.ndarray, frames: int,
                        time_info, status):
        """sounddevice 回调 (在独立线程中运行)"""
        if status:
            print(f"音频状态: {status}", file=sys.stderr)

        # float32 → int16
        pcm_int16 = np.clip(indata[:, 0] * MAX_AMP,
                            -MAX_AMP, MAX_AMP - 1).astype(np.int16)

        # 算法: 计算 25-bin PSD (全程 16-bit 模拟, 无归一化)
        bin_psd, mag_db = compute_bin_psd(pcm_int16, self.hm, self.window_gain,
                                  self.step_table)

        # 发送到 GUI 线程
        self.psd_ready.emit(bin_psd, mag_db, self.threshold_db)

    # ══════════════════════════════════════════════
    #  FFT 尺寸切换
    # ══════════════════════════════════════════════

    def _on_fft_size_changed(self, new_size: int):
        """用户改变了 FFT 尺寸 → 重建窗、步进表、重启音频流"""
        # 对齐到有效的 2^N
        size = 1
        while size < new_size:
            size <<= 1
        if size > 2048:
            size = 2048

        # 避免重复触发
        if size == self.fft_size:
            return
        self.fft_size = size

        # 更新步进表
        self.step_table = build_step_table(size)

        # 重建窗
        self.hm, self.window_gain = init_window(size, self.window_type)

        # 重置音频暂存
        self._pcm_buf = np.zeros(size, dtype=np.int16)
        self._float_ring = np.zeros(size, dtype=np.float64)
        self._ring_pos = 0

        # 重启音频流
        if self.audio_stream:
            self.audio_stream.stop()
            self.audio_stream.close()
            self.audio_stream = None
        self._setup_audio()

        # 同步 spin 显示
        self.fft_spin.blockSignals(True)
        self.fft_spin.setValue(size)
        self.fft_spin.blockSignals(False)

        self.status_label.setText(f"FFT 切换至 {size} pt")

    # ══════════════════════════════════════════════
    #  窗函数切换
    # ══════════════════════════════════════════════

    def _on_window_changed(self, win: str):
        """用户切换窗函数 → 重建窗系数, 重启音频流"""
        if win == self.window_type:
            return
        self.window_type = win
        self.hm, self.window_gain = init_window(self.fft_size, win)

        # 重启音频流
        if self.audio_stream:
            self.audio_stream.stop()
            self.audio_stream.close()
            self.audio_stream = None
        self._setup_audio()

        self.status_label.setText(f"窗: {win}")

    # ══════════════════════════════════════════════
    #  平滑控制
    # ══════════════════════════════════════════════

    def _on_smooth_changed(self, val: int):
        self.smooth_tau_ms = float(val)
        self.smooth_label.setText(f"{val} ms")

    # ══════════════════════════════════════════════
    #  滑条回调
    # ══════════════════════════════════════════════

    def _on_threshold_changed(self):
        self._enforce_monotonic()
        self._update_labels()

    def _enforce_monotonic(self):
        """保证阈值严格递减: 0≥1≥2≥...≥6

        当用户拖动某个滑条时, 相邻滑条会被自动推挤以维持递减顺序。
        """
        vals = [sl.value() for sl in self.sliders]

        # 正向: 每个 ≤ 前一个
        for i in range(1, LED_ROW_COUNT):
            if vals[i] > vals[i - 1]:
                vals[i] = vals[i - 1]

        # 反向: 每个 ≥ 后一个
        for i in range(LED_ROW_COUNT - 2, -1, -1):
            if vals[i] < vals[i + 1]:
                vals[i] = vals[i + 1]

        # 回写到滑块 (阻断信号以避免递归)
        for i in range(LED_ROW_COUNT):
            changed = self.sliders[i].value() != vals[i]
            if changed:
                self.sliders[i].blockSignals(True)
                self.sliders[i].setValue(vals[i])
                self.sliders[i].blockSignals(False)
            self.threshold_db[i] = vals[i]

    def _update_labels(self):
        """更新滑条下方的 dB / 功率标签"""
        for row in range(LED_ROW_COUNT):
            db_val = self.threshold_db[row]
            self.slider_db_labels[row].setText(f"{db_val} dB")
            p = db_to_power(db_val)
            self.slider_power_labels[row].setText(format_power(p))

    # ══════════════════════════════════════════════
    #  显示刷新
    # ══════════════════════════════════════════════

    @pyqtSlot(object, object, object)
    def _on_psd_ready(self, bin_psd, mag_db, thresholds_db):
        """接收音频线程的 PSD 及原始幅度数据"""
        # 从时间常数 (ms) 计算 EMA 因子 α = 1 - exp(-Δt / τ)
        tau = self.smooth_tau_ms
        if tau <= 0.0:
            a = 1.0
        else:
            dt = self.fft_size / SAMPLE_RATE          # 每帧间隔 (秒)
            a = 1.0 - math.exp(-dt / (tau / 1000.0))

        if a >= 1.0:
            self.latest_psd = bin_psd
            self.latest_mag_db = mag_db
        else:
            # 确保平滑缓冲区尺寸匹配
            if len(bin_psd) != len(self._smooth_psd):
                self._smooth_psd = np.zeros_like(bin_psd)
            if len(mag_db) != len(self._smooth_mag_db):
                self._smooth_mag_db = np.zeros_like(mag_db)
            # 指数平滑: y = α·x + (1-α)·y_prev
            self._smooth_psd = a * bin_psd + (1.0 - a) * self._smooth_psd
            self._smooth_mag_db = a * mag_db + (1.0 - a) * self._smooth_mag_db
            self.latest_psd = self._smooth_psd
            self.latest_mag_db = self._smooth_mag_db
        self.threshold_db = thresholds_db[:]

    def _refresh_display(self):
        """定时器: 更新 LED 矩阵 + 折线频谱"""
        bin_psd = self.latest_psd

        # 更新折线频谱 (原始 FFT 幅度, 非 PSD)
        self.line_spectrum.update_data(self.latest_mag_db, self._palette)

        # 更新 LED 矩阵
        for col in range(SPECTRUM_BIN_SIZE):
            psd_val = bin_psd[col]
            for row in range(LED_ROW_COUNT):
                threshold_power = db_to_power(self.threshold_db[row])
                self.led_widgets[row][col].set_on(psd_val > threshold_power)

        # 状态栏
        max_db = power_to_db(np.max(bin_psd)) if np.max(bin_psd) > 0 else -100
        self.status_label.setText(
            f"峰值: {format_power(np.max(bin_psd))}  "
            f"({max_db:.1f} dB)  |  "
            f"活跃 bin: {np.sum(bin_psd > 0)}/{SPECTRUM_BIN_SIZE}"
        )

    # ══════════════════════════════════════════════
    #  清理
    # ══════════════════════════════════════════════

    def closeEvent(self, event):
        self.timer.stop()
        if self.audio_stream:
            self.audio_stream.stop()
            self.audio_stream.close()
        super().closeEvent(event)


# ════════════════════════════════════════════════════════════════
#  入口
# ════════════════════════════════════════════════════════════════
def main():
    # 预检: 检查音频设备
    try:
        sd.check_input_settings(samplerate=SAMPLE_RATE, channels=1)
    except Exception as e:
        print(f"音频设备检查失败: {e}", file=sys.stderr)
        print("将以离线模式运行 (可使用测试滑条手动测试)", file=sys.stderr)

    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
