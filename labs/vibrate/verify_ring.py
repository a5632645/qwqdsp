"""Drive the four-diameter plus outer-ring lab and capture its scene/scope."""

from __future__ import annotations

import math
from pathlib import Path

from verify_kit import Session, Timeline

ROOT = Path(__file__).resolve().parents[3]
EXE = (
    ROOT / "build" / "qwqdsp" / "labs" / "vibrate" / "Release"
    / "qwqdsp-lab-vibrate-ring.exe"
)
OUT = ROOT / "temp" / "ring_verify"


def points(session: Session):
    width = session.client_rect().width
    height = session.client_rect().height
    box_x = width - 168.0
    box_y = 70.0
    inner_x = box_x + 10.0
    inner_y = box_y + 10.0
    inner = 124.0
    mid_x = inner_x + inner * 0.5
    mid_y = inner_y + inner * 0.5

    ring = (
        mid_x + 0.98 * inner * 0.5 * math.cos(math.pi / 4),
        mid_y - 0.98 * inner * 0.5 * math.sin(math.pi / 4),
    )
    string0 = (mid_x + 0.50 * inner * 0.5, mid_y)
    string1 = (
        mid_x + 0.70 * inner * 0.5 * math.cos(math.pi / 4),
        mid_y - 0.70 * inner * 0.5 * math.sin(math.pi / 4),
    )

    track_x = width - 250.0
    track_w = 200.0
    speed_a = track_x + (0.25 - 0.10) / 0.90 * track_w
    speed_b = track_x + (0.70 - 0.10) / 0.90 * track_w

    def frac(p: tuple[float, float]) -> tuple[float, float]:
        return (p[0] / width, p[1] / height)

    return (
        frac(ring),
        frac(string0),
        frac(string1),
        frac((speed_a, 40.0)),
        frac((speed_b, 40.0)),
    )


def main() -> None:
    session = Session(EXE, cwd=ROOT / "build", title="Four Strings Ring")
    session.start()
    ring, string0, string1, speed_a, speed_b = points(session)

    timeline = Timeline(session, close=True)
    timeline.at(0.3).screenshot(OUT / "00_early.png")
    timeline.at(0.6).screenshot(OUT / "01_scene.png")
    timeline.drag(1.0, speed_a, speed_b, duration=0.3, steps=8)
    timeline.at(1.5).screenshot(OUT / "02_speed.png")
    timeline.click(2.0, string0)
    timeline.at(2.6).screenshot(OUT / "03_string0.png")
    timeline.click(3.2, string1)
    timeline.at(3.8).screenshot(OUT / "04_string1.png")
    timeline.click(4.4, ring)
    timeline.at(5.0).screenshot(OUT / "05_ring.png")
    timeline.at(5.4).key("r")
    timeline.at(5.7).screenshot(OUT / "06_reset.png")
    timeline.at(6.4).screenshot(OUT / "07_after.png")
    timeline.run(until=6.6)
    print(f"screenshots: {OUT}")


if __name__ == "__main__":
    main()
