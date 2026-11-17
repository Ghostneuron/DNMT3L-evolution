#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd
from Bio import SeqIO


PROJECT = Path(__file__).resolve().parents[1]
WORK = PROJECT / "analysis/gene_family_tree_337"
SOURCE = WORK / "DNMT3_family_337_MAFFT.fasta"
OUTPUT = WORK / "DNMT3_family_337_MAFFT_human_DNMT3L_columns.fasta"
REFERENCE = "DNMT3L__Homo_sapiens"


def main() -> None:
    records = list(SeqIO.parse(SOURCE, "fasta"))
    lengths = {len(record.seq) for record in records}
    if len(lengths) != 1:
        raise ValueError(f"Unequal alignment lengths: {lengths}")
    references = [str(record.seq).upper() for record in records if record.id == REFERENCE]
    if len(references) != 1:
        raise ValueError(f"Expected one {REFERENCE}, found {len(references)}")
    retained = [index for index, residue in enumerate(references[0]) if residue != "-"]
    if len(retained) != 387:
        raise ValueError(f"Expected 387 human DNMT3L columns, found {len(retained)}")

    qc_rows = []
    with OUTPUT.open("w") as handle:
        for record in records:
            sequence = str(record.seq).upper()
            trimmed = "".join(sequence[index] for index in retained)
            nongap = sum(residue not in {"-", "X", "?"} for residue in trimmed)
            qc_rows.append(
                {
                    "tip": record.id,
                    "group": record.id.split("__", 1)[0],
                    "reference_columns": len(trimmed),
                    "nongap_unambiguous": nongap,
                    "reference_column_coverage": nongap / len(trimmed),
                }
            )
            handle.write(f">{record.id}\n")
            for start in range(0, len(trimmed), 80):
                handle.write(trimmed[start : start + 80] + "\n")

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(WORK / "DNMT3_family_337_alignment_coverage.tsv", sep="\t", index=False)
    print(qc.groupby("group")["reference_column_coverage"].describe().to_string())


if __name__ == "__main__":
    main()
