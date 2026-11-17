#!/usr/bin/env python3
"""Summarize the lungfish DNMT3-family placement sensitivities."""

from pathlib import Path

import pandas as pd
from Bio import Phylo


PROJECT = Path(__file__).resolve().parents[1]
TREE_DIR = PROJECT / "analysis/locus_history_337"
OUTPUT = PROJECT / "results/lungfish_dnmt3_family_tree_summary.tsv"
CANDIDATE = "LUNGFISH_candidate_translated_segment"

TREES = [
    ("targeted", "full segment; internal stop masked as X", "lungfish_targeted_candidate_panel_full_masked_stop_human_columns_LG_F_R4.treefile"),
    ("targeted", "pre-stop segment", "lungfish_targeted_candidate_panel_prestop_human_columns_LG_F_R4.treefile"),
    ("expanded", "full segment; internal stop masked as X", "lungfish_expanded_family_full_masked_stop_LG_F_R4_fast.treefile"),
    ("expanded", "pre-stop segment", "lungfish_expanded_family_prestop_LG_F_R4_fast.treefile"),
]


def main() -> None:
    rows = []
    for panel, treatment, filename in TREES:
        tree = Phylo.read(TREE_DIR / filename, "newick")
        candidate = next(tip for tip in tree.get_terminals() if tip.name == CANDIDATE)
        path = tree.get_path(candidate)
        parent = path[-2]
        sisters = [
            tip.name
            for child in parent.clades
            if child is not candidate
            for tip in child.get_terminals()
        ]
        sister_families = sorted({name.split("__", 1)[0] for name in sisters})
        dnmt3l_only_ancestors = []
        for ancestor in path[:-1]:
            families = {
                tip.name.split("__", 1)[0]
                for tip in ancestor.get_terminals()
                if tip.name != CANDIDATE
            }
            if families == {"DNMT3L"}:
                dnmt3l_only_ancestors.append(ancestor)
        smallest_dnmt3l_only = min(
            (len(node.get_terminals()) for node in dnmt3l_only_ancestors),
            default=0,
        )
        if panel == "targeted":
            interpretation = "DNMT3L-side placement next to sampled squamate DNMT3L references"
        elif treatment == "pre-stop segment":
            interpretation = "DNMT3L-side placement, but the immediate attachment is unsupported and long-branched"
        else:
            interpretation = "unresolved in the expanded panel; no DNMT3L-only attachment is recovered"
        rows.append(
            {
                "panel": panel,
                "candidate_treatment": treatment,
                "candidate_branch_length": candidate.branch_length,
                "immediate_parent_sh_alrt": parent.confidence,
                "sister_clade_size": len(sisters),
                "sister_family_labels": ";".join(sister_families),
                "sister_tips": ";".join(sisters),
                "smallest_dnmt3l_only_clade_size_including_candidate": smallest_dnmt3l_only,
                "interpretation": interpretation,
                "tree_file": str((TREE_DIR / filename).relative_to(PROJECT)),
            }
        )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUTPUT, sep="\t", index=False)
    print(OUTPUT)


if __name__ == "__main__":
    main()
