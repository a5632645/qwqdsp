"""Drive the circular-membrane lab and capture its colormap and oscilloscope."""

from __future__ import annotations

from pathlib import Path

from verify_kit import Session, Timeline

ROOT = Path(__file__).resolve().parents[3]
EXE = (
    ROOT / "build" / "qwqdsp" / "labs" / "vibrate" / "Release"
    / "qwqdsp-lab-vibrate-disk.exe"
)
OUT = ROOT / "temp" / "disk_verify"


def points(session: Session):
    width = session.client_rect().width
    height = session.client_rect().height
    scene_height = height * 0.63
    size = min(width * 0.50, scene_height * 0.82)
    left = width * 0.5 - size * 0.5
    top = scene_height * 0.54 - size * 0.5

    a = (left + 0.35 * size, top + 0.35 * size)
    b = (left + 0.65 * size, top + 0.60 * size)

    track_x = width - 250.0
    track_w = 200.0
    speed_a = track_x + (0.25 - 0.10) / 0.90 * track_w
    speed_b = track_x + (0.70 - 0.10) / 0.90 * track_w

    def frac(p: tuple[float, float]) -> tuple[float, float]:
        return (p[0] / width, p[1] / height)

    return frac(a), frac(b), frac((speed_a, 40.0)), frac((speed_b, 40.0))


def main() -> None:
    session = Session(EXE, cwd=ROOT / "build", title="Circular Membrane")
    session.start()
    probe_a, probe_b, speed_a, speed_b = points(session)

    timeline = Timeline(session, close=True)
    timeline.at(0.3).screenshot(OUT / "00_early.png")
    timeline.at(0.6).screenshot(OUT / "01_scene.png")
    timeline.drag(1.0, speed_a, speed_b, duration=0.3, steps=8)
    timeline.at(1.5).screenshot(OUT / "02_speed.png")
    timeline.click(2.0, probe_a)
    timeline.at(2.6).screenshot(OUT / "03_probe_a.png")
    timeline.drag(3.2, probe_a, probe_b, duration=0.4, steps=12)
    timeline.at(4.0).screenshot(OUT / "04_probe_b.png")
    timeline.at(4.4).key("r")
    timeline.at(4.7).screenshot(OUT / "05_reset.png")
    timeline.at(5.4).screenshot(OUT / "06_after.png")
    timeline.run(until=5.6)
    print(f"screenshots: {OUT}")


if __name__ == "__main__":
    main()
