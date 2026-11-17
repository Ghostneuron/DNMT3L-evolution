#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

from Bio import SeqIO


PROJECT = Path(__file__).resolve().parents[1]
MATRIX = PROJECT / "results/DNMT3L_337_orthology_evidence_matrix.tsv"
REFERENCE_PANEL = PROJECT / "data/processed/dnmt3_representatives.fasta"
OUTPUT = PROJECT / "analysis/locus_history_337/targeted_partial_candidate_panel.fasta"
SPECIES = (
    "Ambystoma mexicanum",
    "Geotrypetes seraphini",
    "Alligator mississippiensis",
    "Crocodylus porosus",
    "Anolis carolinensis",
    "Pogona vitticeps",
    "Pelodiscus sinensis",
    "Monodelphis domestica",
    "Mus musculus",
    "Homo sapiens",
)


def write_record(handle, name: str, sequence: str) -> None:
    handle.write(f">{name}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    with MATRIX.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_species = {row["species"]: row for row in rows}
    missing = set(SPECIES) - set(by_species)
    if missing:
        raise ValueError(f"Missing DNMT3L panel species: {sorted(missing)}")
    paralogs = [
        record
        for record in SeqIO.parse(REFERENCE_PANEL, "fasta")
        if record.id.startswith(("DNMT3A__", "DNMT3B__", "DNMT3C__"))
    ]
    with OUTPUT.open("w") as handle:
        for species in SPECIES:
            row = by_species[species]
            write_record(handle, f"DNMT3L__{row['safe_id']}", row["protein_seq"])
        for record in paralogs:
            write_record(handle, record.id, str(record.seq).upper())
    print(f"Wrote {OUTPUT.relative_to(PROJECT)}: {len(SPECIES)} DNMT3L candidates + {len(paralogs)} paralogs")


if __name__ == "__main__":
    main()
