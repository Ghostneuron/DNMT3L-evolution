#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd
from Bio import Phylo


PROJECT = Path(__file__).resolve().parents[1]
WORK = PROJECT / "analysis/gene_family_tree_337"
TREE = WORK / "DNMT3_family_337_LG_F_R4_fast.treefile"
MATRIX = PROJECT / "results/DNMT3L_337_preliminary_orthology_matrix.tsv"
RESULTS = PROJECT / "results"
PROCESSED = PROJECT / "data/processed/ncbi_337_representatives"


def exact_split(tree, target: set[str]) -> tuple[bool, float | None]:
    all_tips = {tip.name for tip in tree.get_terminals()}
    for clade in tree.find_clades(order="level"):
        descendants = {tip.name for tip in clade.get_terminals()}
        if descendants == target or all_tips - descendants == target:
            return True, clade.confidence
    return False, None


def write_fasta(rows: pd.DataFrame, column: str, path: Path) -> None:
    with path.open("w") as handle:
        for row in rows.itertuples(index=False):
            sequence = str(getattr(row, column))
            handle.write(f">{row.safe_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    tree = Phylo.read(TREE, "newick")
    tips = [tip.name for tip in tree.get_terminals()]
    groups = {
        group: {tip for tip in tips if tip.startswith(f"{group}__")}
        for group in ["DNMT3L", "DNMT3A", "DNMT3B"]
    }
    catalytic = groups["DNMT3A"] | groups["DNMT3B"]
    summary_rows = []
    split_status: dict[str, bool] = {}
    for label, members in [
        ("DNMT3L", groups["DNMT3L"]),
        ("DNMT3A", groups["DNMT3A"]),
        ("DNMT3B", groups["DNMT3B"]),
        ("DNMT3A_plus_DNMT3B", catalytic),
    ]:
        separated, support = exact_split(tree, members)
        split_status[label] = separated
        summary_rows.append(
            {
                "group": label,
                "n_tips": len(members),
                "exact_root_independent_split": separated,
                "SH_aLRT_support": support,
            }
        )

    distance_rows = []
    for candidate in sorted(groups["DNMT3L"]):
        other_l = groups["DNMT3L"] - {candidate}
        nearest_l = min(other_l, key=lambda tip: tree.distance(candidate, tip))
        nearest_l_distance = tree.distance(candidate, nearest_l)
        nearest_paralog = min(catalytic, key=lambda tip: tree.distance(candidate, tip))
        nearest_paralog_distance = tree.distance(candidate, nearest_paralog)
        distance_rows.append(
            {
                "candidate_tip": candidate,
                "nearest_DNMT3L_tip": nearest_l,
                "nearest_DNMT3L_distance": nearest_l_distance,
                "nearest_DNMT3A_or_DNMT3B_tip": nearest_paralog,
                "nearest_DNMT3A_or_DNMT3B_distance": nearest_paralog_distance,
                "nearest_group_is_DNMT3L": nearest_l_distance < nearest_paralog_distance,
                "paralog_to_DNMT3L_distance_ratio": (
                    nearest_paralog_distance / nearest_l_distance
                    if nearest_l_distance > 0
                    else float("inf")
                ),
            }
        )
    summary = pd.DataFrame(summary_rows)
    distances = pd.DataFrame(distance_rows).sort_values("paralog_to_DNMT3L_distance_ratio")
    summary.to_csv(RESULTS / "DNMT3_family_337_group_separation_summary.tsv", sep="\t", index=False)
    distances.to_csv(RESULTS / "DNMT3_family_337_candidate_distance_audit.tsv", sep="\t", index=False)

    matrix = pd.read_csv(MATRIX, sep="\t", dtype={"gene_id": str}).fillna("")
    nearest_pass = {
        row.candidate_tip.removeprefix("DNMT3L__"): bool(row.nearest_group_is_DNMT3L)
        for row in distances.itertuples(index=False)
    }
    statuses = []
    for row in matrix.itertuples(index=False):
        if not bool(row.Preliminary_Alignment_Ready):
            statuses.append("NOT_TESTED_NOT_ALIGNMENT_READY")
        elif not split_status["DNMT3L"]:
            statuses.append("FAIL_GROUP_NOT_EXACT_SPLIT")
        elif nearest_pass.get(row.safe_id, False):
            statuses.append("PASS")
        else:
            statuses.append("FAIL_NEAREST_PARALOG")
    matrix["Family_Tree_Orthology_Status"] = statuses
    matrix["Selection_Ready"] = (
        matrix["Preliminary_Alignment_Ready"].eq(True)
        & matrix["Family_Tree_Orthology_Status"].eq("PASS")
    )
    matrix.to_csv(RESULTS / "DNMT3L_337_orthology_evidence_matrix.tsv", sep="\t", index=False)

    ready = matrix.loc[matrix["Selection_Ready"]].copy()
    write_fasta(ready, "protein_seq", PROCESSED / "DNMT3L_337_selection_ready_all_protein.fasta")
    write_fasta(ready, "cds_seq", PROCESSED / "DNMT3L_337_selection_ready_all_cds.fasta")
    subset_definitions = {
        "Mammalia": {"Afrotheria", "Euarchontoglires", "Laurasiatheria", "Xenarthra", "Marsupialia"},
        "Placentalia": {"Afrotheria", "Euarchontoglires", "Laurasiatheria", "Xenarthra"},
        "Euarchontoglires": {"Euarchontoglires"},
        "Laurasiatheria": {"Laurasiatheria"},
        "Sauropsida": {"Sauropsida"},
        "Amphibia": {"Amphibia"},
    }
    subset_rows = []
    for subset, clades in subset_definitions.items():
        selected = ready.loc[ready["clade"].isin(clades)].copy()
        write_fasta(selected, "protein_seq", PROCESSED / f"DNMT3L_337_selection_ready_{subset}_protein.fasta")
        write_fasta(selected, "cds_seq", PROCESSED / f"DNMT3L_337_selection_ready_{subset}_cds.fasta")
        subset_rows.append({"subset": subset, "included_clades": ";".join(sorted(clades)), "n": len(selected)})
    pd.DataFrame(subset_rows).to_csv(
        RESULTS / "ncbi_337_selection_ready_subsets.tsv", sep="\t", index=False
    )

    print(summary.to_string(index=False))
    print("\nCandidates closer to a paralog reference")
    print(distances.loc[~distances["nearest_group_is_DNMT3L"]].to_string(index=False))
    print(f"\nFinal selection-ready records: {len(ready)}")


if __name__ == "__main__":
    main()
