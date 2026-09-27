#!/usr/bin/env python3
"""
PSD to LED — 功率谱 → 类频谱刻度可视化工具

功率谱由硬件 FFT 计算而来, 始终位于 0 ~ 32768² 范围。
每个 LED 通道:
  ① LED 指示灯 (测试功率 > LED 阈值功率时亮起)
  ② 幅值滑条 (-80 ~ 0 dB, 10^(dB/20) 换算功率)
  ③ dB 读数
  ④ 阈值功率文本
底部测试滑条 -100 ~ 0 dB, 拉动时计算测试功率并刷新所有 LED。
"""

import sys
import math
from typing import Optional

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QSlider, QLabel, QFrame, QSizePolicy, QStatusBar,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QPainter, QColor, QBrush, QPen, QFont

# ──────────────────────────────────────────────
#  常量
# ──────────────────────────────────────────────
MAX_AMP = 32768
MAX_POWER = MAX_AMP * MAX_AMP          # 1 073 741 824
LED_COUNT = 7
LED_DB_MIN = -80
LED_DB_MAX = 0
TEST_DB_MIN = -100
TEST_DB_MAX = 0


def amp_db_to_power(db: float) -> float:
    """幅值 dB → 功率谱

    公式: power = 32768² × 10^(dB / 20)
    """
    return MAX_POWER * (10.0 ** (db / 20.0))


def power_to_amp_db(power: float) -> float:
    """功率谱 → 幅值 dB, 钳位至 [-100, 0]"""
    if power <= 0.0:
        return float(TEST_DB_MIN)
    db = 20.0 * math.log10(power / MAX_POWER)
    return max(float(TEST_DB_MIN), min(0.0, db))


def format_power(power: float) -> str:
    """友好格式化功率值"""
    if power >= 1e9:
        return f"{power / 1e9:.3f} G"
    if power >= 1e6:
        return f"{power / 1e6:.3f} M"
    if power >= 1e3:
        return f"{power / 1e3:.3f} K"
    return f"{power:.1f}"


# ──────────────────────────────────────────────
#  LED 指示灯
# ──────────────────────────────────────────────
class LedWidget(QWidget):
    """圆形 LED, 亮: 绿色, 灭: 深灰"""

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._on = False
        self.setFixedSize(36, 36)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)

    def set_on(self, state: bool) -> None:
        if self._on != state:
            self._on = state
            self.update()

    def is_on(self) -> bool:
        return self._on

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        if self._on:
            fill = QColor(0, 220, 0)
            border = QColor(0, 180, 0)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QBrush(QColor(0, 255, 0, 50)))
            painter.drawEllipse(2, 2, 32, 32)
        else:
            fill = QColor(55, 55, 55)
            border = QColor(35, 35, 35)

        painter.setPen(QPen(border, 2))
        painter.setBrush(QBrush(fill))
        painter.drawEllipse(4, 4, 28, 28)

        if self._on:
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QBrush(QColor(180, 255, 180, 100)))
            painter.drawEllipse(9, 9, 18, 18)


# ──────────────────────────────────────────────
#  单个 LED 通道: LED + 幅值滑条 + 标签
# ──────────────────────────────────────────────
class LedChannel(QFrame):
    """一个 LED 的阈值标定控件"""

    def __init__(self, index: int, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._index = index
        self._test_power: float = 0.0

        self.setFrameStyle(QFrame.Shape.StyledPanel | QFrame.Shadow.Raised)
        self.setFixedWidth(100)

        layout = QVBoxLayout(self)
        layout.setSpacing(4)
        layout.setContentsMargins(8, 8, 8, 10)

        # ── 编号 ──
        idx_lbl = QLabel(f"LED {index}")
        idx_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        f = QFont("Consolas", 10)
        f.setBold(True)
        idx_lbl.setFont(f)
        layout.addWidget(idx_lbl)

        # ── LED ──
        self.led = LedWidget()
        layout.addWidget(self.led, alignment=Qt.AlignmentFlag.AlignCenter)

        # ── 幅值滑条 (-80 ~ 0 dB) ──
        self.slider = QSlider(Qt.Orientation.Vertical)
        self.slider.setRange(LED_DB_MIN, LED_DB_MAX)
        self.slider.setValue(-40)
        self.slider.setTickPosition(QSlider.TickPosition.TicksRight)
        self.slider.setTickInterval(10)
        self.slider.setFixedHeight(140)
        self.slider.valueChanged.connect(self._on_slider_changed)
        layout.addWidget(self.slider, alignment=Qt.AlignmentFlag.AlignCenter)

        # ── dB 读数 ──
        self.db_label = QLabel("-40 dB")
        self.db_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        f_db = QFont("Consolas", 10)
        self.db_label.setFont(f_db)
        layout.addWidget(self.db_label)

        # ── 功率阈值文本 ──
        self.power_label = QLabel()
        self.power_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.power_label.setWordWrap(True)
        f_pw = QFont("Consolas", 9)
        self.power_label.setFont(f_pw)
        layout.addWidget(self.power_label)

        self._update_power_label()
        self._update_led()

    # ── public ──

    def set_test_power(self, power: float) -> None:
        self._test_power = power
        self._update_led()

    @property
    def threshold_power(self) -> float:
        return amp_db_to_power(float(self.slider.value()))

    @property
    def threshold_db(self) -> int:
        return self.slider.value()

    # ── internal ──

    def _on_slider_changed(self, value: int) -> None:
        self.db_label.setText(f"{value} dB")
        self._update_power_label()
        self._update_led()

    def _update_power_label(self) -> None:
        p = amp_db_to_power(float(self.slider.value()))
        self.power_label.setText(format_power(p))

    def _update_led(self) -> None:
        self.led.set_on(self._test_power > self.threshold_power)


# ──────────────────────────────────────────────
#  主窗口
# ──────────────────────────────────────────────
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PSD to LED — 功率谱阈值标定")
        self.setMinimumSize(800, 520)
        self.setMaximumHeight(580)
        self._setup_ui()

    def _setup_ui(self) -> None:
        central = QWidget()
        self.setCentralWidget(central)
        main_vbox = QVBoxLayout(central)
        main_vbox.setSpacing(8)

        # ═══════════════ 7 个 LED 通道行 ═══════════════
        led_row = QHBoxLayout()
        led_row.setSpacing(6)
        led_row.setContentsMargins(10, 8, 10, 4)

        self.channels: list[LedChannel] = []
        for i in range(LED_COUNT):
            ch = LedChannel(i)
            self.channels.append(ch)
            led_row.addWidget(ch)

        led_row.addStretch()
        main_vbox.addLayout(led_row)

        # ═══════════════ 分隔线 ═══════════════
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        main_vbox.addWidget(line)

        # ═══════════════ 测试幅度滑条 (-100 ~ 0 dB) ═══════════════
        test_group = QVBoxLayout()
        test_group.setContentsMargins(20, 4, 40, 8)

        # 标题行
        title_row = QHBoxLayout()
        title_row.addWidget(QLabel("测试幅度"))
        title_row.addStretch()
        self.test_db_value_label = QLabel("0 dB")
        self.test_db_value_label.setFont(QFont("Consolas", 12, QFont.Weight.Bold))
        title_row.addWidget(self.test_db_value_label)
        test_group.addLayout(title_row)

        # 水平滑条
        slider_row = QHBoxLayout()
        slider_row.addWidget(QLabel("-100 dB"))
        self.test_slider = QSlider(Qt.Orientation.Horizontal)
        self.test_slider.setRange(TEST_DB_MIN, TEST_DB_MAX)
        self.test_slider.setValue(0)
        self.test_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.test_slider.setTickInterval(10)
        self.test_slider.valueChanged.connect(self._on_test_slider_changed)
        slider_row.addWidget(self.test_slider, stretch=1)
        slider_row.addWidget(QLabel("0 dB"))
        test_group.addLayout(slider_row)

        # 信息行
        info_row = QHBoxLayout()
        info_row.addStretch()
        self.test_power_label = QLabel("测试功率: 1.074 G")
        self.test_power_label.setFont(QFont("Consolas", 11))
        info_row.addWidget(self.test_power_label)
        info_row.addStretch()
        test_group.addLayout(info_row)

        main_vbox.addLayout(test_group)

        # ═══════════════ 状态栏 ═══════════════
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_label = QLabel("就绪")
        self.status_bar.addWidget(self.status_label)

        # ── 初始刷新 ──
        self._update_test(0)

    # ══════════════════════════════════════════
    #  测试幅度变化
    # ══════════════════════════════════════════

    def _on_test_slider_changed(self, value: int) -> None:
        self._update_test(value)

    def _update_test(self, db_value: int) -> None:
        test_power = amp_db_to_power(float(db_value))

        # 更新标签
        self.test_db_value_label.setText(f"{db_value} dB")
        self.test_power_label.setText(
            f"测试功率: {format_power(test_power)}"
        )

        # 刷新所有 LED
        for ch in self.channels:
            ch.set_test_power(test_power)

        # 更新状态栏
        lit = sum(1 for ch in self.channels if ch.led.is_on())
        self.status_label.setText(
            f"LED 亮起: {lit}/{LED_COUNT}  |  "
            f"测试幅度: {db_value} dB  |  "
            f"测试功率: {format_power(test_power)}"
        )


# ──────────────────────────────────────────────
#  入口
# ──────────────────────────────────────────────
def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
