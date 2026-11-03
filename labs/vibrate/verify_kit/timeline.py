"""Absolute-time action scheduling for a Session."""

from __future__ import annotations

import time
from collections.abc import Callable
from pathlib import Path

from .session import Session

Action = tuple[float, int, Callable[[], None], str]


class Timeline:
    """Schedule actions at absolute seconds, measured from window focus."""

    def __init__(self, session: Session, *, close: bool = True, verbose: bool = True) -> None:
        self.session = session
        self.close = close
        self.verbose = verbose
        self._actions: list[Action] = []
        self._sequence = 0

    def at(self, seconds: float) -> Step:
        return Step(self, float(seconds))

    def _add(self, seconds: float, action: Callable[[], None], label: str) -> None:
        self._actions.append((seconds, self._sequence, action, label))
        self._sequence += 1

    def drag(
        self,
        start_time: float,
        start: tuple[float, float],
        end: tuple[float, float],
        *,
        duration: float = 0.3,
        steps: int = 12,
        button: str = "left",
        release_delay: float = 0.05,
    ) -> Timeline:
        if steps < 1:
            raise ValueError("steps must be >= 1")
        self.at(start_time).move(*start).mouse_down(button)
        x0, y0 = start
        x1, y1 = end
        for index in range(1, steps + 1):
            fraction = index / steps
            self.at(start_time + duration * fraction).move(
                x0 + (x1 - x0) * fraction,
                y0 + (y1 - y0) * fraction,
            )
        self.at(start_time + duration + release_delay).mouse_up(button)
        return self

    def click(self, seconds: float, point: tuple[float, float], button: str = "left") -> Timeline:
        self.at(seconds).move(*point).mouse_down(button)
        self.at(seconds + 0.04).mouse_up(button)
        return self

    def key(self, seconds: float, name: str) -> Timeline:
        self.at(seconds).key(name)
        return self

    def run(self, *, until: float | None = None) -> Timeline:
        self.session.start()
        self.session.focus()
        started = time.monotonic()
        try:
            for seconds, _sequence, action, label in sorted(
                self._actions, key=lambda item: (item[0], item[1])
            ):
                delay = started + seconds - time.monotonic()
                if delay > 0:
                    time.sleep(delay)
                if self.verbose:
                    elapsed = time.monotonic() - started
                    print(f"[{elapsed:6.3f}s] {label}", flush=True)
                action()
            if until is not None:
                delay = started + until - time.monotonic()
                if delay > 0:
                    time.sleep(delay)
        finally:
            self.session.wait_screenshots()
            if self.close:
                self.session.close()
        return self


class Step:
    """Actions sharing one timestamp. Methods return self for chaining."""

    def __init__(self, timeline: Timeline, seconds: float) -> None:
        self._timeline = timeline
        self._seconds = seconds

    def move(self, fx: float, fy: float) -> Step:
        self._timeline._add(
            self._seconds,
            lambda: self._timeline.session.move(fx, fy),
            f"move ({fx:.3f}, {fy:.3f})",
        )
        return self

    def move_px(self, x: float, y: float) -> Step:
        self._timeline._add(
            self._seconds,
            lambda: self._timeline.session.move_px(x, y),
            f"move_px ({x:.0f}, {y:.0f})",
        )
        return self

    def mouse_down(self, button: str = "left") -> Step:
        self._timeline._add(
            self._seconds,
            lambda: self._timeline.session.mouse_down(button),
            f"mouse_down {button}",
        )
        return self

    def mouse_up(self, button: str = "left") -> Step:
        self._timeline._add(
            self._seconds,
            lambda: self._timeline.session.mouse_up(button),
            f"mouse_up {button}",
        )
        return self

    def key(self, name: str) -> Step:
        self._timeline._add(
            self._seconds,
            lambda: self._timeline.session.key(name),
            f"key {name}",
        )
        return self

    def screenshot(self, path: str | Path) -> Step:
        target = Path(path)
        self._timeline._add(
            self._seconds,
            lambda: self._timeline.session.screenshot_async(target),
            f"screenshot {target}",
        )
        return self
