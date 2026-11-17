#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "scripts/build_ancestral_state_figure_337.py"
SPEC = importlib.util.spec_from_file_location("ancestral_state_v06", SOURCE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not load {SOURCE}")
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)

OUTDIR = PROJECT / "figures/current_v07"
BASE.OUT_PNG = OUTDIR / "Figure_4_ancestral_state_context_caveated.png"
BASE.OUT_PDF = OUTDIR / "Figure_4_ancestral_state_context_caveated.pdf"
ORIGINAL_DRAW_TREE = BASE.draw_tree


def draw_tree_with_mapping_caveat(ax, nodes) -> None:
    ORIGINAL_DRAW_TREE(ax, nodes)
    ax.text(
        0.05,
        0.01,
        "Mapping caveat: all four amphibian NP349 calls fail the local-context identity gate; "
        "after omitting Geotrypetes, Amniota K posterior = 0.518 / 0.499.",
        transform=ax.transAxes,
        fontsize=7.7,
        color="#991B1B",
        fontweight="bold",
    )


BASE.draw_tree = draw_tree_with_mapping_caveat


if __name__ == "__main__":
    OUTDIR.mkdir(parents=True, exist_ok=True)
    BASE.main()
