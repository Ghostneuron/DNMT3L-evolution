from __future__ import annotations

import importlib.util
import shutil
from tempfile import TemporaryDirectory
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "analysis" / "selection_reanalysis_337"
INPUTS = ROOT / "inputs"
OUTPUT = ROOT / "saturation_diagnostics"
RESULT = PROJECT / "results" / "selection_337_codon_saturation_summary.tsv"
ARCHIVED_SCRIPT = (
    PROJECT / "../../../source_package/scripts/assess_codon_saturation.py"
).resolve()
ANALYSES = ["CONSENSUS", "MACSE", "MAFFT", "PRANK", "MUSCLE5", "COBALT"]


def main() -> None:
    with TemporaryDirectory(prefix="dnmt3l_saturation_337_") as temporary:
        staging = Path(temporary)
        for analysis in ANALYSES:
            source = (
                INPUTS
                / f"Placentalia_{analysis}_human_coordinate_occupancy80_deduplicated.fasta"
            )
            destination = staging / f"Placentalia_{analysis}_selection_input.fasta"
            shutil.copyfile(source, destination)

        spec = importlib.util.spec_from_file_location("archived_saturation", ARCHIVED_SCRIPT)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Cannot load {ARCHIVED_SCRIPT}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.ALIGNMENT_DIR = staging
        module.OUT = OUTPUT
        module.SEED = 20261102
        module.MAX_PAIRS = 3000
        module.MIN_VALID_CODONS = 100
        module.main()

    summary = pd.read_csv(OUTPUT / "codon_saturation_summary.csv")
    summary.to_csv(RESULT, sep="\t", index=False)
    print(f"Wrote {RESULT.relative_to(PROJECT)}")


if __name__ == "__main__":
    main()
