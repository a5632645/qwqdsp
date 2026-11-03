"""Drive the 1D string lab and capture a deterministic screenshot sequence.

Run from anywhere:
    python qwqdsp/labs/vibrate/verify_chain_1d.py
"""

from __future__ import annotations

from pathlib import Path

from verify_kit import Session, Timeline

ROOT = Path(__file__).resolve().parents[3]
EXE = ROOT / "build" / "qwqdsp" / "labs" / "vibrate" / "Release" / "qwqdsp-lab-vibrate-chain1d.exe"
OUT = ROOT / "temp" / "vibrate_verify"


def layout_fractions(session: Session) -> tuple[float, float]:
    height = session.client_rect().height
    panel_height = min(max(height * 0.26, 150.0), 240.0)
    slider_y = height - 44.0 - panel_height - 26.0
    rest_y = (185.0 + slider_y) * 0.5
    return rest_y / height, slider_y / height


def main() -> None:
    session = Session(EXE, cwd=ROOT / "build", title="1D Vibration")
    session.start()
    rest, slider = layout_fractions(session)

    timeline = Timeline(session, close=True)
    timeline.at(0.8).screenshot(OUT / "01_initial.png")
    timeline.drag(1.2, (0.60, rest), (0.60, rest - 0.14), duration=0.35, steps=16)
    timeline.at(1.8).screenshot(OUT / "02_pluck.png")
    timeline.at(2.8).screenshot(OUT / "03_ring.png")
    timeline.drag(3.2, (0.50, slider), (0.75, slider), duration=0.4, steps=12)
    timeline.at(4.0).screenshot(OUT / "04_probe.png")
    timeline.at(4.3).key("a")
    timeline.at(5.4).screenshot(OUT / "05_drive.png")
    timeline.run(until=5.6)
    print(f"screenshots: {OUT}")


if __name__ == "__main__":
    main()
