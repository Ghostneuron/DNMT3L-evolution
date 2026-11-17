#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "analysis/codon_alignments_337"
RESULTS = PROJECT / "results"
CLADES = ["Euarchontoglires", "Laurasiatheria", "Sauropsida", "Amphibia"]


def read_fasta(path: Path) -> dict[str, str]:
    return {record.id: str(record.seq).upper().replace("U", "T") for record in SeqIO.parse(path, "fasta")}


def allowed(sequence: str) -> set[str]:
    values = {sequence}
    if str(Seq(sequence[-3:]).translate()) == "*":
        values.add(sequence[:-3])
    return values


def main() -> None:
    rows = []
    for clade in CLADES:
        source = read_fasta(ROOT / f"inputs/{clade}_cds_with_coordinate_anchor.fasta")
        alignment = read_fasta(ROOT / f"{clade}/MACSE/{clade}_MACSE_NT.fasta")
        if set(source) != set(alignment):
            raise ValueError(f"{clade} identifier mismatch")
        lengths = {len(sequence) for sequence in alignment.values()}
        if len(lengths) != 1 or next(iter(lengths)) % 3:
            raise ValueError(f"{clade} alignment is not rectangular codon data")
        exact = 0
        frameshifts = 0
        records_with_frameshifts = 0
        for name, sequence in alignment.items():
            ungapped = sequence.replace("-", "")
            record_frameshifts = ungapped.count("!")
            frameshifts += record_frameshifts
            records_with_frameshifts += int(record_frameshifts > 0)
            if ungapped.replace("!", "") in allowed(source[name]):
                exact += 1
        reference_id = "Homo_sapiens" if "Homo_sapiens" in alignment else "Homo_sapiens_anchor"
        human = alignment[reference_id]
        mapped = 0
        for start in range(0, len(human), 3):
            codon = human[start : start + 3]
            if codon == "---":
                continue
            if set(codon) <= set("ACGT") and str(Seq(codon).translate()) == "*":
                continue
            if set(codon) - set("ACGT"):
                raise ValueError(f"{clade} has unresolved human anchor codon")
            mapped += 1
        if mapped != 387 or exact != len(source):
            raise ValueError(f"{clade} failed human mapping or source-content gate")
        rows.append(
            {
                "clade": clade,
                "biological_records": len(source) - int(reference_id.endswith("_anchor")),
                "alignment_records": len(source),
                "human_coordinate_anchor": reference_id,
                "alignment_codon_columns": next(iter(lengths)) // 3,
                "human_mapped_residues": mapped,
                "exact_source_content_records": exact,
                "frameshift_symbols": frameshifts,
                "records_with_frameshift_symbols": records_with_frameshifts,
                "status": "PASS_RAW_REQUIRES_MASKING" if frameshifts else "PASS",
            }
        )
    output = pd.DataFrame(rows)
    output.to_csv(RESULTS / "clade_MACSE_337_alignment_audit.tsv", sep="\t", index=False)
    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
