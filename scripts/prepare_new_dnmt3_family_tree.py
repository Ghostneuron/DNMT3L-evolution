#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd
from Bio import SeqIO


PROJECT = Path(__file__).resolve().parents[1]
MATRIX = PROJECT / "results/DNMT3L_337_preliminary_orthology_matrix.tsv"
OLD_FAMILY_FASTA = (
    PROJECT
    / "../../../source_package/analysis/orthology_rebuild_20260826/gene_family_tree/"
    "DNMT3_family_candidates_plus_paralog_references.fasta"
).resolve()
OUT = PROJECT / "analysis/gene_family_tree_337"


def write_record(handle, name: str, sequence: str) -> None:
    handle.write(f">{name}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_csv(MATRIX, sep="\t", dtype={"gene_id": str}).fillna("")
    candidates = matrix.loc[matrix["Preliminary_Alignment_Ready"].eq(True)].copy()
    paralogs = [
        record
        for record in SeqIO.parse(OLD_FAMILY_FASTA, "fasta")
        if record.id.startswith(("DNMT3A__", "DNMT3B__"))
    ]
    if len(paralogs) != 21:
        raise ValueError(f"Expected 21 frozen DNMT3A/B references, found {len(paralogs)}")

    output = OUT / "DNMT3_family_337_candidates_plus_paralog_references.fasta"
    with output.open("w") as handle:
        for row in candidates.itertuples(index=False):
            write_record(handle, f"DNMT3L__{row.safe_id}", row.protein_seq)
        for record in paralogs:
            write_record(handle, record.id, str(record.seq).upper())

    reference_rows = []
    for record in paralogs:
        parts = record.id.split("__")
        reference_rows.append(
            {
                "tip": record.id,
                "group": parts[0],
                "species": parts[1].replace("_", " "),
                "protein_length": len(record.seq),
                "source": "frozen_20260826_reference_panel",
            }
        )
    pd.DataFrame(reference_rows).to_csv(
        OUT / "DNMT3A_DNMT3B_reference_panel.tsv", sep="\t", index=False
    )
    print(f"DNMT3L preliminary-ready candidates: {len(candidates)}")
    print(f"DNMT3A/B references: {len(paralogs)}")
    print(f"Combined FASTA: {output}")


if __name__ == "__main__":
    main()
