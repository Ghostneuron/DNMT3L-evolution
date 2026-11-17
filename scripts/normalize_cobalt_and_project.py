#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def codons_for_cds(sequence: str) -> list[str]:
    sequence = sequence.upper().replace("U", "T")
    if len(sequence) % 3:
        raise ValueError("CDS length is not divisible by three")
    codons = [sequence[start : start + 3] for start in range(0, len(sequence), 3)]
    if codons and str(Seq(codons[-1]).translate()) == "*":
        codons.pop()
    return codons


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument("--cds", type=Path, required=True)
    parser.add_argument("--raw-alignment", type=Path, required=True)
    parser.add_argument("--named-alignment", type=Path, required=True)
    parser.add_argument("--codon-alignment", type=Path, required=True)
    parser.add_argument("--qc-json", type=Path, required=True)
    args = parser.parse_args()

    proteins = [(record.id, str(record.seq).upper()) for record in SeqIO.parse(args.protein, "fasta")]
    cdss = {record.id: str(record.seq).upper() for record in SeqIO.parse(args.cds, "fasta")}
    aligned = [(record.id, str(record.seq).upper()) for record in SeqIO.parse(args.raw_alignment, "fasta")]
    if len(proteins) != len(cdss) or set(dict(proteins)) != set(cdss):
        raise ValueError("Protein and CDS identifiers are not one-to-one")
    if len(aligned) != len(proteins):
        raise ValueError(f"COBALT emitted {len(aligned)} of {len(proteins)} records")
    if len({len(sequence) for _, sequence in aligned}) != 1:
        raise ValueError("COBALT output is not rectangular")

    named: list[tuple[str, str]] = []
    projected: list[tuple[str, str]] = []
    for index, ((raw_id, aligned_sequence), (expected_id, source_protein)) in enumerate(
        zip(aligned, proteins), start=1
    ):
        if raw_id not in {f"lcl|{index}", str(index), expected_id}:
            raise ValueError(f"Unexpected COBALT identifier {raw_id}; expected local id {index}")
        if aligned_sequence.replace("-", "") != source_protein:
            raise ValueError(f"COBALT changed residue content for {expected_id}")
        named.append((expected_id, aligned_sequence))
        codons = codons_for_cds(cdss[expected_id])
        translated = "".join(str(Seq(codon).translate()) for codon in codons)
        if translated != source_protein:
            raise ValueError(f"Protein-CDS translation mismatch for {expected_id}")
        codon_index = 0
        output: list[str] = []
        for residue in aligned_sequence:
            if residue == "-":
                output.append("---")
            else:
                output.append(codons[codon_index])
                codon_index += 1
        if codon_index != len(codons):
            raise ValueError(f"Unused codons remain for {expected_id}")
        projected.append((expected_id, "".join(output)))

    write_fasta(named, args.named_alignment)
    write_fasta(projected, args.codon_alignment)
    aligned_columns = len(named[0][1])
    qc = {
        "sequences": len(named),
        "unique_identifiers": len({name for name, _ in named}),
        "aligned_amino_acid_columns": aligned_columns,
        "aligned_codon_columns": aligned_columns,
        "ungapped_protein_content_preserved_for_all_sequences": True,
        "exact_cds_projection_for_all_sequences": True,
        "raw_alignment_sha256": sha256(args.raw_alignment),
        "named_alignment_sha256": sha256(args.named_alignment),
        "codon_alignment_sha256": sha256(args.codon_alignment),
        "cobalt_mode": "local_similarity_only_-norps_T",
    }
    args.qc_json.write_text(json.dumps(qc, indent=2) + "\n")
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
