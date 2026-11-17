#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from Bio import SeqIO


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = (PROJECT / "../../../analysis_v2/orthology_sensitivity/DNMT3_family_plus_DNMT3C.fasta").resolve()
OUTPUT = PROJECT / "data" / "processed" / "dnmt3_representatives.fasta"
TARGETS = {
    "DNMT3L__Homo_sapiens",
    "DNMT3L__Anolis_carolinensis",
    "DNMT3L__Alligator_mississippiensis",
    "DNMT3L__Ambystoma_mexicanum",
    "DNMT3A__Homo_sapiens__NP_072046_2",
    "DNMT3B__Homo_sapiens__XP_047295902_1",
    "DNMT3C__Mus_musculus__P0DOY1",
}


def main() -> None:
    records = [record for record in SeqIO.parse(SOURCE, "fasta") if record.id in TARGETS]
    recovered = {record.id for record in records}
    missing = TARGETS - recovered
    if missing:
        raise ValueError(f"Missing representative sequences: {sorted(missing)}")
    for record in records:
        record.description = ""
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT.with_suffix(OUTPUT.suffix + ".tmp")
    SeqIO.write(records, temporary, "fasta")
    temporary.replace(OUTPUT)
    print(f"Wrote {OUTPUT.relative_to(PROJECT)} ({len(records)} DNMT3-family representatives)")


if __name__ == "__main__":
    main()

