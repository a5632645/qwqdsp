"""一键跑完 AAIIR-Osc 论文复现的全部实验。

用法: python qwqdsp/labs/adaa_iir/run_all.py
"""
from __future__ import annotations

import runpy
import sys
from pathlib import Path

HERE = Path(__file__).parent
SCRIPTS = ("check_impl.py", "exp_saw_snr.py", "exp_wavetable_sweep.py",
           "exp_time_domain.py", "exp_cost.py",
           # ADAA-IIR waveshaper 探索（见 ADAA-IIR-waveshaper.md）
           "check_ws.py", "exp_ws_snr.py", "exp_ws_cost.py",
           "exp_ws_phase.py", "exp_ws_tanh.py")


def main() -> int:
    sys.path.insert(0, str(HERE))
    for name in SCRIPTS:
        print(f"\n{'=' * 70}\n== {name}\n{'=' * 70}")
        try:
            runpy.run_path(str(HERE / name), run_name="__main__")
        except SystemExit as exc:
            if exc.code:
                print(f"!! {name} 以退出码 {exc.code} 结束")
                return int(exc.code)
    print("\n全部实验完成。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
