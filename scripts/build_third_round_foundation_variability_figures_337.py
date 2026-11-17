#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "scripts/build_redesign_figures_337.py"
SPEC = importlib.util.spec_from_file_location("redesign_figures_v08", SOURCE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not load {SOURCE}")
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
BASE.OUT = PROJECT / "figures/current_v08"


if __name__ == "__main__":
    BASE.style()
    BASE.figure_1()
    BASE.figure_2()
    print(BASE.OUT)
