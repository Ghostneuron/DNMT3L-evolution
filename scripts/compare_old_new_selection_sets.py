#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
OLD = (
    PROJECT
    / "../../../source_package/analysis/orthology_rebuild_20260826/orthology_audit/"
    "DNMT3L_orthology_evidence_matrix.csv"
).resolve()
NEW = PROJECT / "results/DNMT3L_337_orthology_evidence_matrix.tsv"
RESULTS = PROJECT / "results"


def main() -> None:
    old = pd.read_csv(OLD, dtype={"gene_id": str}).fillna("")
    new = pd.read_csv(NEW, sep="\t", dtype={"gene_id": str}).fillna("")
    old_by_gene = old.set_index("gene_id")
    new_by_gene = new.set_index("gene_id")
    rows = []
    for gene_id in sorted(set(old_by_gene.index) | set(new_by_gene.index), key=int):
        old_row = old_by_gene.loc[gene_id] if gene_id in old_by_gene.index else None
        new_row = new_by_gene.loc[gene_id] if gene_id in new_by_gene.index else None
        old_ready = bool(old_row["Selection_Ready"]) if old_row is not None else False
        new_ready = bool(new_row["Selection_Ready"]) if new_row is not None else False
        if old_row is None:
            membership = "NEW_GENEID"
        elif new_row is None:
            membership = "REMOVED_GENEID"
        else:
            membership = "SHARED_GENEID"
        if not old_ready and new_ready:
            readiness_change = "GAINED_READY"
        elif old_ready and not new_ready:
            readiness_change = "LOST_READY"
        elif old_ready and new_ready:
            readiness_change = "READY_BOTH"
        else:
            readiness_change = "NOT_READY_BOTH"
        old_accession = str(old_row["protein_accession"]) if old_row is not None else ""
        new_accession = str(new_row["protein_accession"]) if new_row is not None else ""
        rows.append(
            {
                "gene_id": gene_id,
                "membership": membership,
                "old_species": str(old_row["species"]) if old_row is not None else "",
                "new_species": str(new_row["species"]) if new_row is not None else "",
                "old_protein_accession": old_accession,
                "new_protein_accession": new_accession,
                "protein_accession_changed": bool(old_accession and new_accession and old_accession != new_accession),
                "old_selection_ready": old_ready,
                "new_selection_ready": new_ready,
                "readiness_change": readiness_change,
            }
        )
    comparison = pd.DataFrame(rows)
    comparison.to_csv(RESULTS / "old_321_vs_new_337_selection_set_comparison.tsv", sep="\t", index=False)
    summary = pd.DataFrame(
        [
            ("old_candidate_GeneIDs", old["gene_id"].nunique()),
            ("new_candidate_GeneIDs", new["gene_id"].nunique()),
            ("old_selection_ready", int(old["Selection_Ready"].astype(bool).sum())),
            ("new_selection_ready", int(new["Selection_Ready"].astype(bool).sum())),
            ("shared_GeneIDs", int(comparison["membership"].eq("SHARED_GENEID").sum())),
            ("new_GeneIDs", int(comparison["membership"].eq("NEW_GENEID").sum())),
            ("removed_GeneIDs", int(comparison["membership"].eq("REMOVED_GENEID").sum())),
            ("gained_ready", int(comparison["readiness_change"].eq("GAINED_READY").sum())),
            ("lost_ready", int(comparison["readiness_change"].eq("LOST_READY").sum())),
            ("ready_in_both", int(comparison["readiness_change"].eq("READY_BOTH").sum())),
            (
                "accession_changed_among_ready_in_both",
                int(
                    (
                        comparison["readiness_change"].eq("READY_BOTH")
                        & comparison["protein_accession_changed"]
                    ).sum()
                ),
            ),
        ],
        columns=["metric", "value"],
    )
    summary.to_csv(RESULTS / "old_321_vs_new_337_selection_set_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print("\nGained readiness")
    print(comparison.loc[comparison["readiness_change"].eq("GAINED_READY")].to_string(index=False))
    print("\nLost readiness")
    print(comparison.loc[comparison["readiness_change"].eq("LOST_READY")].to_string(index=False))


if __name__ == "__main__":
    main()
