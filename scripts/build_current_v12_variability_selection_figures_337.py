#!/usr/bin/env python3
"""Regenerate the variability and selection figures with current terminology."""

from pathlib import Path

import build_redesign_figures_337 as figures


def main() -> None:
    figures.OUT = Path(__file__).resolve().parents[1] / "figures/current_v12"
    figures.style()
    figures.figure_1()
    figures.figure_2()
    figures.figure_4()

    canonical_stems = {
        "Figure_1_updated_data_foundation": "Figure_4_data_foundation",
        "Figure_2_five_aligner_variability": "Figure_3_five_aligner_variability",
        "Figure_4_selection_sensitivity": "Supplementary_Figure_S1_selection_sensitivity",
    }
    for source_stem, destination_stem in canonical_stems.items():
        for suffix in (".png", ".pdf"):
            source = figures.OUT / f"{source_stem}{suffix}"
            destination = figures.OUT / f"{destination_stem}{suffix}"
            source.replace(destination)
    print(figures.OUT)


if __name__ == "__main__":
    main()
