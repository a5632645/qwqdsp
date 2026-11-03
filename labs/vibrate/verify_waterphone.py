"""Drive the waterphone lab: strike rods, bow a rod, probe the rim."""

from __future__ import annotations

import math
from pathlib import Path

from verify_kit import Session, Timeline

ROOT = Path(__file__).resolve().parents[3]
EXE = (
    ROOT / "build" / "qwqdsp" / "labs" / "vibrate" / "Release"
    / "qwqdsp-lab-vibrate-waterphone.exe"
)
OUT = ROOT / "temp" / "waterphone_verify"


def points(session: Session):
    width = session.client_rect().width
    height = session.client_rect().height
    scene_height = height * 0.63
    cx = width * 0.5
    cy = scene_height * 0.74
    axis = min(width * 0.26, scene_height * 0.45)
    depth_x = axis * 0.55
    depth_y = -axis * 0.30
    z_scale = scene_height * 0.85
    radius = 0.55

    def rod_mid(degree: float, length: float) -> tuple[float, float]:
        angle = math.radians(degree)
        base_x = cx + radius * axis * math.cos(angle) + \
            radius * depth_x * math.sin(angle)
        base_y = cy + radius * depth_y * math.sin(angle)
        return (base_x, base_y - 0.5 * length * z_scale)

    # rod 18: front, length 0.54 m; rod 6: back, length 0.28 m; rod 0: right.
    front_long = rod_mid(270.0, 0.54)
    back_bow = rod_mid(90.0, 0.28)

    pad_x = width - 168.0 + 10.0 + (148.0 - 20.0) * 0.5
    pad_y = 70.0 + 10.0 + (148.0 - 20.0) * 0.5
    pad_rim = (148.0 - 20.0) * 0.36
    pad_front = (pad_x, pad_y + pad_rim)

    track_x = width - 250.0
    track_w = 200.0
    speed_a = track_x + (0.25 - 0.10) / 0.90 * track_w
    speed_b = track_x + (0.70 - 0.10) / 0.90 * track_w

    def frac(p: tuple[float, float]) -> tuple[float, float]:
        return (p[0] / width, p[1] / height)

    # rod 6 tangent at 90 degrees projects to screen -x; bow by dragging left.
    bow_end = (back_bow[0] - 55.0, back_bow[1])

    return (
        frac(front_long),
        frac(back_bow),
        frac(bow_end),
        frac(pad_front),
        frac((speed_a, 40.0)),
        frac((speed_b, 40.0)),
    )


def main() -> None:
    session = Session(EXE, cwd=ROOT / "build", title="Waterphone")
    session.start()
    strike_rod, bow_rod, bow_end, pad_front, speed_a, speed_b = points(
        session)

    timeline = Timeline(session, close=True)
    timeline.at(0.3).screenshot(OUT / "00_early.png")
    timeline.at(0.6).screenshot(OUT / "01_scene.png")
    timeline.drag(1.0, speed_a, speed_b, duration=0.3, steps=8)
    timeline.at(1.5).screenshot(OUT / "02_speed.png")
    timeline.click(2.0, strike_rod)
    timeline.at(2.6).screenshot(OUT / "03_strike.png")
    timeline.drag(3.2, bow_rod, bow_end, duration=0.6, steps=12)
    timeline.at(3.7).screenshot(OUT / "04_bow.png")
    timeline.at(4.2).screenshot(OUT / "05_bow_after.png")
    timeline.click(4.8, pad_front)
    timeline.at(5.3).screenshot(OUT / "06_probe.png")
    timeline.at(5.5).key("r")
    timeline.at(5.8).screenshot(OUT / "07_reset.png")
    timeline.at(6.4).screenshot(OUT / "08_after.png")
    timeline.run(until=6.6)
    print(f"screenshots: {OUT}")


if __name__ == "__main__":
    main()
