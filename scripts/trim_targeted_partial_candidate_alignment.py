#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "analysis/locus_history_337/targeted_partial_candidate_family_mafft.fasta"
OUTPUT = (
    PROJECT
    / "analysis/locus_history_337/targeted_partial_candidate_family_human_dnmt3l_columns.fasta"
)
REFERENCE = "DNMT3L__Homo_sapiens"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--reference", default=REFERENCE)
    args = parser.parse_args()
    alignment = AlignIO.read(args.input, "fasta")
    reference = next(record for record in alignment if record.id == args.reference)
    retained = [index for index, residue in enumerate(str(reference.seq)) if residue != "-"]
    trimmed_records = [
        record.__class__(
            record.seq.__class__("".join(str(record.seq)[index] for index in retained)),
            id=record.id,
            description="",
        )
        for record in alignment
    ]
    trimmed = MultipleSeqAlignment(trimmed_records)
    if trimmed.get_alignment_length() != 387:
        raise ValueError(
            f"Expected 387 human DNMT3L columns, got {trimmed.get_alignment_length()}"
        )
    AlignIO.write(trimmed, args.output, "fasta")
    output_path = args.output.resolve()
    try:
        display_path = output_path.relative_to(PROJECT)
    except ValueError:
        display_path = output_path
    print(f"Wrote {display_path} with {len(trimmed)} sequences and 387 columns")


if __name__ == "__main__":
    main()
