"""Drive the crossed-strings lab and capture its scene and oscilloscope."""

from __future__ import annotations

from pathlib import Path

from verify_kit import Session, Timeline

ROOT = Path(__file__).resolve().parents[3]
EXE = ROOT / "build" / "qwqdsp" / "labs" / "vibrate" / "Release" / "qwqdsp-lab-vibrate-cross.exe"
OUT = ROOT / "temp" / "cross_verify"


def pad_points(session: Session) -> tuple[tuple[float, float], tuple[float, float]]:
    width = session.client_rect().width
    height = session.client_rect().height
    box_x = width - 168.0
    box_y = 24.0
    inner_x = box_x + 10.0
    inner_y = box_y + 10.0
    inner_size = 124.0

    # x=-0.75 在 X 轴上；y=+0.30 在 Y 轴上。
    x_point = (inner_x + 0.25 * inner_size, inner_y + 0.5 * inner_size)
    y_point = (inner_x + 0.5 * inner_size, inner_y + 0.4 * inner_size)
    return (
        (x_point[0] / width, x_point[1] / height),
        (y_point[0] / width, y_point[1] / height),
    )


def main() -> None:
    session = Session(EXE, cwd=ROOT / "build", title="Two Coupled Strings")
    session.start()
    x_point, y_point = pad_points(session)

    timeline = Timeline(session, close=True)
    timeline.at(0.6).screenshot(OUT / "01_scene.png")
    timeline.click(1.0, x_point)
    timeline.at(1.2).screenshot(OUT / "02_probe_x.png")
    timeline.drag(1.8, x_point, y_point, duration=0.5, steps=16)
    timeline.at(2.4).screenshot(OUT / "03_probe_y.png")
    timeline.at(2.8).key("r")
    timeline.at(3.0).screenshot(OUT / "04_reset.png")
    timeline.at(3.8).screenshot(OUT / "05_after.png")
    timeline.run(until=4.0)
    print(f"screenshots: {OUT}")


if __name__ == "__main__":
    main()
