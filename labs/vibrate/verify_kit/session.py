"""Launch a window, drive the mouse/keyboard, and capture its client area.

The module uses only the Python standard library (Pillow accelerates PNG
encoding when installed). Mouse positions are fractions of the client
rectangle so scripts survive window size and DPI differences.
"""

from __future__ import annotations

import ctypes
import ctypes.wintypes as wt
import struct
import subprocess
import threading
import time
import zlib
from dataclasses import dataclass
from pathlib import Path

try:
    from PIL import Image as _PIL_IMAGE
except ImportError:
    _PIL_IMAGE = None

_user32 = ctypes.WinDLL("user32", use_last_error=True)
_gdi32 = ctypes.WinDLL("gdi32", use_last_error=True)
_kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)

SW_RESTORE = 9
WM_CLOSE = 0x0010
SRCCOPY = 0x00CC0020
DIB_RGB_COLORS = 0
BI_RGB = 0
SWP_NOMOVE = 0x0002
SWP_NOSIZE = 0x0001
KEYEVENTF_KEYUP = 0x0002
MOUSEEVENTF_LEFTDOWN, MOUSEEVENTF_LEFTUP = 0x0002, 0x0004
MOUSEEVENTF_RIGHTDOWN, MOUSEEVENTF_RIGHTUP = 0x0008, 0x0010
MOUSEEVENTF_MIDDLEDOWN, MOUSEEVENTF_MIDDLEUP = 0x0020, 0x0040

_BUTTONS = {
    "left": (MOUSEEVENTF_LEFTDOWN, MOUSEEVENTF_LEFTUP),
    "right": (MOUSEEVENTF_RIGHTDOWN, MOUSEEVENTF_RIGHTUP),
    "middle": (MOUSEEVENTF_MIDDLEDOWN, MOUSEEVENTF_MIDDLEUP),
}

_VK = {
    "space": 0x20,
    "enter": 0x0D,
    "esc": 0x1B,
    "escape": 0x1B,
    "tab": 0x09,
    "backspace": 0x08,
    "shift": 0x10,
    "ctrl": 0x11,
    "alt": 0x12,
    "up": 0x26,
    "down": 0x28,
    "left": 0x25,
    "right": 0x27,
}


class RECT(ctypes.Structure):
    _fields_ = [
        ("left", wt.LONG),
        ("top", wt.LONG),
        ("right", wt.LONG),
        ("bottom", wt.LONG),
    ]


class POINT(ctypes.Structure):
    _fields_ = [("x", wt.LONG), ("y", wt.LONG)]


class BITMAPINFOHEADER(ctypes.Structure):
    _fields_ = [
        ("biSize", wt.DWORD),
        ("biWidth", wt.LONG),
        ("biHeight", wt.LONG),
        ("biPlanes", wt.WORD),
        ("biBitCount", wt.WORD),
        ("biCompression", wt.DWORD),
        ("biSizeImage", wt.DWORD),
        ("biXPelsPerMeter", wt.LONG),
        ("biYPelsPerMeter", wt.LONG),
        ("biClrUsed", wt.DWORD),
        ("biClrImportant", wt.DWORD),
    ]


class BITMAPINFO(ctypes.Structure):
    _fields_ = [("bmiHeader", BITMAPINFOHEADER), ("bmiColors", wt.DWORD * 3)]


_ENUM_PROC = ctypes.WINFUNCTYPE(wt.BOOL, wt.HWND, wt.LPARAM)

_user32.EnumWindows.argtypes = [_ENUM_PROC, wt.LPARAM]
_user32.EnumWindows.restype = wt.BOOL
_user32.GetWindowThreadProcessId.argtypes = [wt.HWND, ctypes.POINTER(wt.DWORD)]
_user32.GetWindowThreadProcessId.restype = wt.DWORD
_user32.GetWindowTextLengthW.argtypes = [wt.HWND]
_user32.GetWindowTextLengthW.restype = ctypes.c_int
_user32.GetWindowTextW.argtypes = [wt.HWND, wt.LPWSTR, ctypes.c_int]
_user32.GetWindowTextW.restype = ctypes.c_int
_user32.IsWindowVisible.argtypes = [wt.HWND]
_user32.IsWindowVisible.restype = wt.BOOL
_user32.GetWindowRect.argtypes = [wt.HWND, ctypes.POINTER(RECT)]
_user32.GetWindowRect.restype = wt.BOOL
_user32.GetClientRect.argtypes = [wt.HWND, ctypes.POINTER(RECT)]
_user32.GetClientRect.restype = wt.BOOL
_user32.ClientToScreen.argtypes = [wt.HWND, ctypes.POINTER(POINT)]
_user32.ClientToScreen.restype = wt.BOOL
_user32.SetForegroundWindow.argtypes = [wt.HWND]
_user32.SetForegroundWindow.restype = wt.BOOL
_user32.GetForegroundWindow.restype = wt.HWND
_user32.BringWindowToTop.argtypes = [wt.HWND]
_user32.BringWindowToTop.restype = wt.BOOL
_user32.SetFocus.argtypes = [wt.HWND]
_user32.SetFocus.restype = wt.HWND
_user32.ShowWindow.argtypes = [wt.HWND, ctypes.c_int]
_user32.ShowWindow.restype = wt.BOOL
_user32.AttachThreadInput.argtypes = [wt.DWORD, wt.DWORD, wt.BOOL]
_user32.AttachThreadInput.restype = wt.BOOL
_user32.SetWindowPos.argtypes = [
    wt.HWND,
    wt.HWND,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    wt.UINT,
]
_user32.SetWindowPos.restype = wt.BOOL
_user32.SetCursorPos.argtypes = [ctypes.c_int, ctypes.c_int]
_user32.SetCursorPos.restype = wt.BOOL
_user32.mouse_event.argtypes = [wt.DWORD, wt.DWORD, wt.DWORD, wt.DWORD, ctypes.c_size_t]
_user32.mouse_event.restype = None
_user32.keybd_event.argtypes = [wt.BYTE, wt.BYTE, wt.DWORD, ctypes.c_size_t]
_user32.keybd_event.restype = None
_user32.PostMessageW.argtypes = [wt.HWND, wt.UINT, wt.WPARAM, wt.LPARAM]
_user32.PostMessageW.restype = wt.BOOL
_user32.GetDC.argtypes = [wt.HWND]
_user32.GetDC.restype = wt.HDC
_user32.ReleaseDC.argtypes = [wt.HWND, wt.HDC]
_user32.ReleaseDC.restype = ctypes.c_int
_gdi32.CreateCompatibleDC.argtypes = [wt.HDC]
_gdi32.CreateCompatibleDC.restype = wt.HDC
_gdi32.CreateCompatibleBitmap.argtypes = [wt.HDC, ctypes.c_int, ctypes.c_int]
_gdi32.CreateCompatibleBitmap.restype = wt.HBITMAP
_gdi32.SelectObject.argtypes = [wt.HDC, wt.HGDIOBJ]
_gdi32.SelectObject.restype = wt.HGDIOBJ
_gdi32.BitBlt.argtypes = [
    wt.HDC,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    wt.HDC,
    ctypes.c_int,
    ctypes.c_int,
    wt.DWORD,
]
_gdi32.BitBlt.restype = wt.BOOL
_gdi32.GetDIBits.argtypes = [
    wt.HDC,
    wt.HBITMAP,
    wt.UINT,
    wt.UINT,
    ctypes.c_void_p,
    ctypes.POINTER(BITMAPINFO),
    wt.UINT,
]
_gdi32.GetDIBits.restype = ctypes.c_int
_gdi32.DeleteObject.argtypes = [wt.HGDIOBJ]
_gdi32.DeleteObject.restype = wt.BOOL
_gdi32.DeleteDC.argtypes = [wt.HDC]
_gdi32.DeleteDC.restype = wt.BOOL
_kernel32.GetCurrentThreadId.restype = wt.DWORD

try:
    _user32.SetProcessDpiAwarenessContext.argtypes = [ctypes.c_void_p]
    _user32.SetProcessDpiAwarenessContext(ctypes.c_void_p(-4))
except AttributeError:
    try:
        _user32.SetProcessDPIAware()
    except AttributeError:
        pass


@dataclass(frozen=True)
class Rect:
    left: int
    top: int
    width: int
    height: int

    @property
    def right(self) -> int:
        return self.left + self.width

    @property
    def bottom(self) -> int:
        return self.top + self.height


def _window_text(hwnd: int) -> str:
    length = _user32.GetWindowTextLengthW(wt.HWND(hwnd))
    if length <= 0:
        return ""
    buffer = ctypes.create_unicode_buffer(length + 1)
    _user32.GetWindowTextW(wt.HWND(hwnd), buffer, length + 1)
    return buffer.value


def _find_window(pid: int, title: str | None) -> int | None:
    candidates: list[tuple[int, int]] = []

    @_ENUM_PROC
    def callback(hwnd, _lparam):
        window_pid = wt.DWORD()
        _user32.GetWindowThreadProcessId(hwnd, ctypes.byref(window_pid))
        if window_pid.value != pid or not _user32.IsWindowVisible(hwnd):
            return True
        text = _window_text(hwnd)
        if not text:
            return True
        if title is not None and title.lower() not in text.lower():
            return True
        rect = RECT()
        _user32.GetWindowRect(hwnd, ctypes.byref(rect))
        area = (rect.right - rect.left) * (rect.bottom - rect.top)
        candidates.append((area, int(hwnd)))
        return True

    _user32.EnumWindows(callback, 0)
    if not candidates:
        return None
    return max(candidates)[1]


def _wait_for_window(process: subprocess.Popen, title: str | None, timeout: float) -> int:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise RuntimeError(f"process exited before opening a window (code {process.returncode})")
        hwnd = _find_window(process.pid, title)
        if hwnd is not None:
            return hwnd
        time.sleep(0.05)
    raise TimeoutError(f"window not found within {timeout:.1f}s (title={title!r})")


def _bgra_to_rgb(bgra: bytes) -> bytes:
    rgb = bytearray((len(bgra) // 4) * 3)
    for index in range(len(bgra) // 4):
        source = index * 4
        target = index * 3
        rgb[target] = bgra[source + 2]
        rgb[target + 1] = bgra[source + 1]
        rgb[target + 2] = bgra[source]
    return bytes(rgb)


def _write_png_bgra(path: Path, width: int, height: int, bgra: bytes) -> None:
    if _PIL_IMAGE is not None:
        try:
            image = _PIL_IMAGE.frombytes("RGB", (width, height), bgra, "raw", "BGRX", 0, 1)
            image.save(path, format="PNG")
            return
        except Exception:
            pass
    _write_png(path, width, height, _bgra_to_rgb(bgra))


def _write_png(path: Path, width: int, height: int, rgb: bytes) -> None:
    stride = width * 3
    rows = bytearray((stride + 1) * height)
    for y in range(height):
        source = y * stride
        target = y * (stride + 1)
        rows[target] = 0  # filter: none
        rows[target + 1 : target + 1 + stride] = rgb[source : source + stride]

    def chunk(tag: bytes, data: bytes) -> bytes:
        return (
            struct.pack(">I", len(data))
            + tag
            + data
            + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)
        )

    payload = b"\x89PNG\r\n\x1a\n"
    payload += chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
    payload += chunk(b"IDAT", zlib.compress(bytes(rows), 6))
    payload += chunk(b"IEND", b"")
    path.write_bytes(payload)


class Session:
    """A launched process plus its main window."""

    def __init__(
        self,
        exe: str | Path,
        *,
        cwd: str | Path | None = None,
        title: str | None = None,
        args: tuple[str, ...] = (),
        startup_timeout: float = 10.0,
    ) -> None:
        self.exe = Path(exe)
        self.cwd = Path(cwd) if cwd is not None else self.exe.parent
        self.title = title
        self.args = tuple(str(arg) for arg in args)
        self.startup_timeout = startup_timeout
        self._process: subprocess.Popen | None = None
        self._hwnd: int | None = None
        self._screen_threads: list[threading.Thread] = []
        self._screen_errors: list[BaseException] = []

    @property
    def hwnd(self) -> int:
        if self._hwnd is None:
            raise RuntimeError("session is not started")
        return self._hwnd

    @property
    def process(self) -> subprocess.Popen:
        if self._process is None:
            raise RuntimeError("session is not started")
        return self._process

    def start(self) -> Session:
        if self._process is not None:
            return self
        if not self.exe.exists():
            raise FileNotFoundError(self.exe)
        self._process = subprocess.Popen([str(self.exe), *self.args], cwd=str(self.cwd))
        self._hwnd = _wait_for_window(self._process, self.title, self.startup_timeout)
        self.focus()
        return self

    def focus(self) -> None:
        hwnd = wt.HWND(self.hwnd)
        _user32.ShowWindow(hwnd, SW_RESTORE)
        _user32.BringWindowToTop(hwnd)
        if _user32.GetForegroundWindow() == self.hwnd:
            return

        foreground = _user32.GetForegroundWindow()
        current_thread = _kernel32.GetCurrentThreadId()
        foreground_thread = wt.DWORD()
        attached = False
        if foreground:
            _user32.GetWindowThreadProcessId(foreground, ctypes.byref(foreground_thread))
        if foreground_thread.value and foreground_thread.value != current_thread:
            attached = bool(
                _user32.AttachThreadInput(current_thread, foreground_thread.value, True)
            )
        _user32.SetForegroundWindow(hwnd)
        _user32.SetFocus(hwnd)
        if attached:
            _user32.AttachThreadInput(current_thread, foreground_thread.value, False)

        if _user32.GetForegroundWindow() != self.hwnd:
            topmost = ctypes.c_void_p(-1 & ((1 << (8 * ctypes.sizeof(ctypes.c_void_p))) - 1))
            not_topmost = ctypes.c_void_p(-2 & ((1 << (8 * ctypes.sizeof(ctypes.c_void_p))) - 1))
            _user32.SetWindowPos(hwnd, topmost, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE)
            _user32.SetWindowPos(hwnd, not_topmost, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE)
        time.sleep(0.1)

    def client_rect(self) -> Rect:
        rect = RECT()
        if not _user32.GetClientRect(wt.HWND(self.hwnd), ctypes.byref(rect)):
            raise OSError("GetClientRect failed")
        origin = POINT(0, 0)
        if not _user32.ClientToScreen(wt.HWND(self.hwnd), ctypes.byref(origin)):
            raise OSError("ClientToScreen failed")
        return Rect(origin.x, origin.y, rect.right - rect.left, rect.bottom - rect.top)

    def point(self, fx: float, fy: float) -> tuple[int, int]:
        rect = self.client_rect()
        x = rect.left + min(max(int(round(fx * rect.width)), 0), rect.width - 1)
        y = rect.top + min(max(int(round(fy * rect.height)), 0), rect.height - 1)
        return x, y

    def point_px(self, x: float, y: float) -> tuple[int, int]:
        rect = self.client_rect()
        px = rect.left + min(max(int(round(x)), 0), rect.width - 1)
        py = rect.top + min(max(int(round(y)), 0), rect.height - 1)
        return px, py

    def move(self, fx: float, fy: float) -> None:
        x, y = self.point(fx, fy)
        if not _user32.SetCursorPos(x, y):
            raise OSError("SetCursorPos failed")

    def move_px(self, x: float, y: float) -> None:
        px, py = self.point_px(x, y)
        if not _user32.SetCursorPos(px, py):
            raise OSError("SetCursorPos failed")

    def mouse_down(self, button: str = "left") -> None:
        if button not in _BUTTONS:
            raise ValueError(f"unknown button: {button}")
        _user32.mouse_event(_BUTTONS[button][0], 0, 0, 0, 0)

    def mouse_up(self, button: str = "left") -> None:
        if button not in _BUTTONS:
            raise ValueError(f"unknown button: {button}")
        _user32.mouse_event(_BUTTONS[button][1], 0, 0, 0, 0)

    def click(self, fx: float, fy: float, button: str = "left") -> None:
        self.move(fx, fy)
        time.sleep(0.02)
        self.mouse_down(button)
        time.sleep(0.03)
        self.mouse_up(button)

    def key_down(self, name: str) -> None:
        _user32.keybd_event(_vk(name), 0, 0, 0)

    def key_up(self, name: str) -> None:
        _user32.keybd_event(_vk(name), 0, KEYEVENTF_KEYUP, 0)

    def key(self, name: str, hold: float = 0.03) -> None:
        self.key_down(name)
        time.sleep(hold)
        self.key_up(name)

    def screenshot(self, path: str | Path) -> Path:
        target = self._prepare_path(path)
        width, height, bgra = self.capture_bgra()
        _write_png_bgra(target, width, height, bgra)
        return target

    def screenshot_async(self, path: str | Path) -> Path:
        """Capture now, encode/write the PNG on a worker thread."""
        target = self._prepare_path(path)
        width, height, bgra = self.capture_bgra()

        def write() -> None:
            try:
                _write_png_bgra(target, width, height, bgra)
            except BaseException as exc:  # keep the error for wait_screenshots()
                self._screen_errors.append(exc)

        thread = threading.Thread(target=write, name=f"png:{target.name}", daemon=True)
        self._screen_threads.append(thread)
        thread.start()
        return target

    def wait_screenshots(self) -> None:
        for thread in self._screen_threads:
            thread.join()
        self._screen_threads.clear()
        if self._screen_errors:
            error = self._screen_errors[0]
            self._screen_errors.clear()
            raise error

    @staticmethod
    def _prepare_path(path: str | Path) -> Path:
        target = Path(path)
        if not target.is_absolute():
            target = Path.cwd() / target
        target.parent.mkdir(parents=True, exist_ok=True)
        return target

    def capture(self) -> tuple[int, int, bytes]:
        width, height, bgra = self.capture_bgra()
        return width, height, _bgra_to_rgb(bgra)

    def capture_bgra(self) -> tuple[int, int, bytes]:
        self.focus()
        rect = self.client_rect()
        width, height = rect.width, rect.height
        if width <= 0 or height <= 0:
            raise RuntimeError(f"invalid client rect: {rect}")

        screen_dc = _user32.GetDC(None)
        memory_dc = _gdi32.CreateCompatibleDC(screen_dc)
        bitmap = _gdi32.CreateCompatibleBitmap(screen_dc, width, height)
        previous = _gdi32.SelectObject(memory_dc, bitmap)
        try:
            if not _gdi32.BitBlt(
                memory_dc, 0, 0, width, height, screen_dc, rect.left, rect.top, SRCCOPY
            ):
                raise OSError("BitBlt failed")
            _gdi32.SelectObject(memory_dc, previous)

            info = BITMAPINFO()
            info.bmiHeader.biSize = ctypes.sizeof(BITMAPINFOHEADER)
            info.bmiHeader.biWidth = width
            info.bmiHeader.biHeight = -height
            info.bmiHeader.biPlanes = 1
            info.bmiHeader.biBitCount = 32
            info.bmiHeader.biCompression = BI_RGB
            pixels = ctypes.create_string_buffer(width * height * 4)
            if not _gdi32.GetDIBits(
                memory_dc, bitmap, 0, height, pixels, ctypes.byref(info), DIB_RGB_COLORS
            ):
                raise OSError("GetDIBits failed")
        finally:
            _gdi32.SelectObject(memory_dc, previous)
            _gdi32.DeleteObject(bitmap)
            _gdi32.DeleteDC(memory_dc)
            _user32.ReleaseDC(None, screen_dc)
        return width, height, pixels.raw

    def close(self, force: bool = False) -> None:
        if self._process is None:
            return
        self.wait_screenshots()
        if self._process.poll() is None:
            if self._hwnd is not None and not force:
                _user32.PostMessageW(wt.HWND(self._hwnd), WM_CLOSE, 0, 0)
                try:
                    self._process.wait(timeout=2.0)
                except subprocess.TimeoutExpired:
                    self._process.terminate()
            else:
                self._process.terminate()
            try:
                self._process.wait(timeout=2.0)
            except subprocess.TimeoutExpired:
                self._process.kill()
                self._process.wait()
        self._process = None
        self._hwnd = None

    def __enter__(self) -> Session:
        return self.start()

    def __exit__(self, *_exc) -> None:
        self.close()


def _vk(name: str) -> int:
    key = name.lower()
    if key in _VK:
        return _VK[key]
    if len(key) == 1:
        if "a" <= key <= "z" or "0" <= key <= "9":
            return ord(key.upper())
    if key.startswith("f") and key[1:].isdigit():
        number = int(key[1:])
        if 1 <= number <= 12:
            return 0x70 + number - 1
    raise ValueError(f"unknown key: {name!r}")
