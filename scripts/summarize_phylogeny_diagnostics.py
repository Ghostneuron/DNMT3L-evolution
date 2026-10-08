from __future__ import annotations

import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo


PROJECT = Path(__file__).resolve().parents[1]
V2 = PROJECT / "analysis_v2"
PHYLOGENIES = V2 / "phylogenies"
INPUTS = V2 / "selection_inputs"
SUMMARIES = V2 / "summaries"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def capture(pattern: str, text: str, cast):
    match = re.search(pattern, text, flags=re.MULTILINE)
    if not match:
        raise ValueError(f"Pattern not found: {pattern}")
    return cast(match.group(1))


def main() -> None:
    rows: list[dict[str, object]] = []
    iqtree_paths = sorted(
        path for path in PHYLOGENIES.glob("*/*.iqtree") if not path.name.startswith("._")
    )
    if len(iqtree_paths) != 12:
        raise ValueError(f"Expected 12 IQ-TREE reports, found {len(iqtree_paths)}")

    for report in iqtree_paths:
        dataset = report.parent.name
        scope, aligner = dataset.rsplit("_", 1)
        text = report.read_text()
        input_match = re.search(r"^Input data: (\d+) sequences with (\d+) codon sites$", text, flags=re.MULTILINE)
        if not input_match:
            raise ValueError(f"Could not parse dimensions from {report}")
        sequence_count = int(input_match.group(1))
        codon_sites = int(input_match.group(2))
        free_parameters = capture(r"^Number of free parameters \(#branches \+ #model parameters\): (\d+)$", text, int)
        near_zero_match = re.search(r"^WARNING: (\d+) near-zero internal branches", text, flags=re.MULTILINE)
        near_zero = int(near_zero_match.group(1)) if near_zero_match else 0
        internal_percent = capture(r"^Sum of internal branch lengths: [0-9.eE+-]+ \(([0-9.eE+-]+)% of tree length\)$", text, float)
        tree_path = report.with_suffix(".treefile")
        alignment_path = INPUTS / f"{dataset}_selection_input.fasta"
        tree = Phylo.read(tree_path, "newick")
        supports = np.array(
            [float(clade.confidence) for clade in tree.get_nonterminals() if clade.confidence is not None],
            dtype=float,
        )
        rows.append(
            {
                "Scope": scope,
                "Aligner": aligner,
                "Dataset": dataset,
                "Sequence_Count": sequence_count,
                "Codon_Sites": codon_sites,
                "Free_Parameters": free_parameters,
                "Parameter_to_Site_Ratio": free_parameters / codon_sites,
                "Parameters_At_Least_Sites": free_parameters >= codon_sites,
                "Near_Zero_Internal_Branches": near_zero,
                "Internal_Branch_Length_Pct": internal_percent,
                "Internal_Support_Values_N": len(supports),
                "Median_SH_aLRT": float(np.median(supports)) if len(supports) else np.nan,
                "Internal_SH_aLRT_ge_80_Pct": float(100 * np.mean(supports >= 80)) if len(supports) else np.nan,
                "Internal_SH_aLRT_ge_95_Pct": float(100 * np.mean(supports >= 95)) if len(supports) else np.nan,
                "Alignment_SHA256": sha256(alignment_path),
                "Tree_SHA256": sha256(tree_path),
            }
        )

    result = pd.DataFrame(rows).sort_values(["Scope", "Aligner"])
    SUMMARIES.mkdir(parents=True, exist_ok=True)
    result.to_csv(SUMMARIES / "phylogeny_diagnostics.csv", index=False)
    print(
        result[
            [
                "Dataset",
                "Sequence_Count",
                "Codon_Sites",
                "Free_Parameters",
                "Parameter_to_Site_Ratio",
                "Near_Zero_Internal_Branches",
                "Median_SH_aLRT",
                "Internal_SH_aLRT_ge_80_Pct",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()

