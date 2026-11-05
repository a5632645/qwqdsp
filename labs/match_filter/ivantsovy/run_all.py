"""一键运行本 lab：正确性自检 + 两张幅度对照图。

用法: python qwqdsp/labs/match_filter/ivantsovy/run_all.py
"""
from __future__ import annotations

import runpy
import sys
from pathlib import Path

HERE = Path(__file__).parent
SCRIPTS = ("check_impl.py", "check_shelf.py", "check_two_point.py",
           "exp_magnitude.py", "exp_four_types.py", "exp_high_shelf.py",
           "exp_two_point.py", "exp_mapping.py")


def main() -> int:
    sys.path.insert(0, str(HERE))
    for name in SCRIPTS:
        print(f"\n{'=' * 72}\n== {name}\n{'=' * 72}")
        try:
            runpy.run_path(str(HERE / name), run_name="__main__")
        except SystemExit as exc:
            if exc.code:
                print(f"!! {name} 以退出码 {exc.code} 结束")
                return int(exc.code)
    print("\n完成。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
