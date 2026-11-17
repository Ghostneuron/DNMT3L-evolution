#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


PROJECT = Path(__file__).resolve().parents[1]
MATRIX = PROJECT / "results/DNMT3L_337_orthology_evidence_matrix.tsv"
FULL_ALIGNMENT = PROJECT / "analysis/gene_family_tree_337/DNMT3_family_337_MAFFT.fasta"
STRICT_PROTEIN = PROJECT / "analysis/switching_helix_asr_337/inputs/DNMT3L_289_human_coordinate_strict_focal_masked.fasta"
OUT = PROJECT / "analysis/switching_helix_asr_337/inputs/DNMT3L_289_human_coordinate_strict_focal_masked_codon.fasta"
AUDIT = PROJECT / "results/switching_helix_337_codon_projection_audit.tsv"
FOCAL = {349, 351, 352, 355}


def write_fasta(records: dict[str, str], path: Path) -> None:
    with path.open("w") as handle:
        for name in sorted(records):
            handle.write(f">{name}\n")
            sequence = records[name]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    with MATRIX.open(newline="") as handle:
        selected = {
            row["safe_id"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["Selection_Ready"] == "True"
        }
    strict = {record.id: str(record.seq).upper() for record in SeqIO.parse(STRICT_PROTEIN, "fasta")}
    full = {
        record.id.removeprefix("DNMT3L__"): str(record.seq).upper()
        for record in SeqIO.parse(FULL_ALIGNMENT, "fasta")
        if record.id.startswith("DNMT3L__")
    }
    if set(strict) != set(selected) or set(full) != set(selected):
        raise ValueError("Selected-record mismatch among matrix and protein alignments")
    human_full = full["Homo_sapiens"]
    human_columns = [index for index, state in enumerate(human_full) if state != "-"]
    if len(human_columns) != 387:
        raise ValueError(f"Expected 387 human columns, observed {len(human_columns)}")

    output: dict[str, str] = {}
    audit_rows: list[dict[str, object]] = []
    for safe_id, row in selected.items():
        protein = row["protein_seq"].upper()
        cds = row["cds_seq"].upper().replace("U", "T")
        aligned = full[safe_id]
        if aligned.replace("-", "") != protein:
            raise ValueError(f"Full alignment does not reproduce selected protein for {safe_id}")
        if len(cds) < len(protein) * 3:
            raise ValueError(f"CDS is shorter than protein for {safe_id}")

        codon_by_column: dict[int, str] = {}
        amino_index = 0
        for column, state in enumerate(aligned):
            if state == "-":
                continue
            codon = cds[amino_index * 3 : amino_index * 3 + 3]
            translated = str(Seq(codon).translate())
            if translated != state:
                raise ValueError(
                    f"Translation mismatch for {safe_id}, residue {amino_index + 1}: "
                    f"{codon}->{translated}, expected {state}"
                )
            codon_by_column[column] = codon
            amino_index += 1
        if amino_index != len(protein):
            raise ValueError(f"Protein traversal mismatch for {safe_id}")

        projected: list[str] = []
        for human_position, column in enumerate(human_columns, start=1):
            protein_state = aligned[column]
            strict_state = strict[safe_id][human_position - 1]
            if human_position in FOCAL and strict_state == "X":
                projected.append("NNN")
            elif protein_state == "-":
                projected.append("---")
            else:
                codon = codon_by_column[column]
                if str(Seq(codon).translate()) != strict_state:
                    raise ValueError(
                        f"Strict projection mismatch for {safe_id} NP{human_position}: "
                        f"{codon} versus {strict_state}"
                    )
                projected.append(codon)
        sequence = "".join(projected)
        if len(sequence) != 1161:
            raise ValueError(f"Projected codon length mismatch for {safe_id}")
        output[safe_id] = sequence

    write_fasta(output, OUT)
    for position in range(1, 388):
        codons = [sequence[(position - 1) * 3 : position * 3] for sequence in output.values()]
        resolved = [codon for codon in codons if len(codon) == 3 and set(codon) <= set("ACGT")]
        stops = sum(str(Seq(codon).translate()) == "*" for codon in resolved)
        audit_rows.append(
            {
                "human_np_position": position,
                "taxa": len(codons),
                "resolved_codons": len(resolved),
                "occupancy": len(resolved) / len(codons),
                "stop_codons": stops,
                "focal_site": position in FOCAL,
            }
        )
    with AUDIT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)
    if any(row["stop_codons"] for row in audit_rows):
        raise ValueError("Projected alignment contains resolved stop codons")
    print(f"Wrote {len(output)} taxa x 387 codons to {OUT.relative_to(PROJECT)}")
    for row in audit_rows:
        if row["focal_site"]:
            print(row)


if __name__ == "__main__":
    main()
