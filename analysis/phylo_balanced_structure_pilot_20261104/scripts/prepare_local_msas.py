#!/usr/bin/env python3
"""Build frozen, fully local MSAs for the DNMT3L structure pilot.

The MSA policy is identical across targets.  DNMT3L rows come from the
validated 289-orthologue human-coordinate alignment.  DNMT3A complex rows
come from the validated MACSE alignment, projected onto human residues
623--912.  Matched species share integer pairing keys; unmatched rows are
unpaired (key -1).  No sequence is sent to an external MSA service.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import SeqIO


AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_text(path: Path, content: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    partial.write_text(content)
    partial.replace(path)


def csv_msa(query: str, rows: list[tuple[int, str]]) -> str:
    if set(query) - AA:
        raise ValueError("Query contains non-standard residues")
    kept = [(-1, query)]
    seen = {query}
    for key, sequence in rows:
        if len(sequence) != len(query):
            raise ValueError("MSA/query length mismatch")
        if set(sequence) - (AA | {"-"}):
            raise ValueError("MSA contains unsupported characters")
        if set(sequence) == {"-"} or sequence in seen:
            continue
        kept.append((key, sequence))
        seen.add(sequence)
    return "key,sequence\n" + "\n".join(f"{key},{sequence}" for key, sequence in kept) + "\n"


def yaml_monomer(sequence: str, msa_relpath: str) -> str:
    return (
        "version: 1\nsequences:\n  - protein:\n      id: A\n"
        f"      sequence: {sequence}\n      msa: {msa_relpath}\n"
    )


def yaml_complex(dnmt3a: str, dnmt3l: str, msa_a: str, msa_b: str) -> str:
    return (
        "version: 1\nsequences:\n"
        "  - protein:\n      id: A\n"
        f"      sequence: {dnmt3a}\n      msa: {msa_a}\n"
        "  - protein:\n      id: B\n"
        f"      sequence: {dnmt3l}\n      msa: {msa_b}\n"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-root", type=Path, required=True)
    parser.add_argument("--dnmt3l-alignment", type=Path, required=True)
    parser.add_argument("--dnmt3a-alignment", type=Path, required=True)
    parser.add_argument("--matched-inventory", type=Path, required=True)
    args = parser.parse_args()

    panel = read_tsv(args.panel_root / "metadata/panel_manifest.tsv")
    crosswalk_rows = read_tsv(args.panel_root / "metadata/human_coordinate_crosswalk.tsv")
    crosswalk: dict[str, list[dict[str, str]]] = {}
    for row in crosswalk_rows:
        crosswalk.setdefault(row["safe_id"], []).append(row)

    dnmt3l_records = list(SeqIO.parse(args.dnmt3l_alignment, "fasta"))
    if len(dnmt3l_records) != 289 or any(len(record.seq) != 387 for record in dnmt3l_records):
        raise ValueError("Expected the frozen 289 x 387 DNMT3L alignment")

    matched_ids = [row["safe_id"] for row in read_tsv(args.matched_inventory)]
    pairing_key = {safe_id: index for index, safe_id in enumerate(sorted(set(matched_ids)))}

    dnmt3a_records = list(SeqIO.parse(args.dnmt3a_alignment, "fasta"))
    human_a = next(record for record in dnmt3a_records if record.id == "Homo_sapiens")
    human_aligned = str(human_a.seq).upper()
    residue_number = 0
    selected_columns = []
    for column, residue in enumerate(human_aligned):
        if residue in AA:
            residue_number += 1
            if 623 <= residue_number <= 912:
                selected_columns.append(column)
    if len(selected_columns) != 290:
        raise ValueError(f"Expected 290 human DNMT3A columns, found {len(selected_columns)}")

    dnmt3a_query = "".join(human_aligned[column] for column in selected_columns)
    dnmt3a_rows = []
    for record in dnmt3a_records:
        aligned = str(record.seq).upper()
        row = "".join(aligned[column] if aligned[column] in AA else "-" for column in selected_columns)
        dnmt3a_rows.append((pairing_key.get(record.id, -1), row))
    msa_dir = args.panel_root / "inputs_local_msa/msa"
    atomic_text(msa_dir / "DNMT3A_human_623_912.csv", csv_msa(dnmt3a_query, dnmt3a_rows))

    summary = []
    for item in panel:
        safe_id = item["safe_id"]
        query = str(next(SeqIO.parse(args.panel_root / f"inputs/monomer_fasta/{safe_id}.fasta", "fasta")).seq)
        rows = sorted(crosswalk[safe_id], key=lambda row: int(row["human_position"]))
        human_to_input = {int(row["human_position"]): int(row["monomer_input_position"]) for row in rows}
        ctd_rows = sorted(
            (row for row in rows if row["complex_chain_B_position"]),
            key=lambda row: int(row["complex_chain_B_position"]),
        )
        ctd_positions = [int(row["human_position"]) for row in ctd_rows]
        ctd_query = "".join(row["monomer_residue"] for row in ctd_rows)

        full_homologs = []
        ctd_homologs = []
        for record in dnmt3l_records:
            aligned = str(record.seq).upper()
            full = ["-"] * len(query)
            for human_position, input_position in human_to_input.items():
                residue = aligned[human_position - 1]
                if residue in AA:
                    full[input_position - 1] = residue
            full_homologs.append((pairing_key.get(record.id, -1), "".join(full)))
            ctd_homologs.append((
                pairing_key.get(record.id, -1),
                "".join(aligned[position - 1] if aligned[position - 1] in AA else "-" for position in ctd_positions),
            ))

        full_msa = f"{safe_id}_DNMT3L_full.csv"
        ctd_msa = f"{safe_id}_DNMT3L_CTD.csv"
        full_csv = csv_msa(query, full_homologs)
        ctd_csv = csv_msa(ctd_query, ctd_homologs)
        atomic_text(msa_dir / full_msa, full_csv)
        atomic_text(msa_dir / ctd_msa, ctd_csv)
        atomic_text(
            args.panel_root / f"inputs_local_msa/monomer_yaml/{safe_id}.yaml",
            yaml_monomer(query, f"../msa/{full_msa}"),
        )
        atomic_text(
            args.panel_root / f"inputs_local_msa/complex_yaml/{safe_id}.yaml",
            yaml_complex(dnmt3a_query, ctd_query, "../msa/DNMT3A_human_623_912.csv", f"../msa/{ctd_msa}"),
        )
        summary.append({
            "safe_id": safe_id,
            "monomer_length": len(query),
            "ctd_length": len(ctd_query),
            "monomer_msa_rows": len(full_csv.splitlines()) - 1,
            "ctd_msa_rows": len(ctd_csv.splitlines()) - 1,
        })

    fields = list(summary[0])
    out = args.panel_root / "metadata/local_msa_summary.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.with_suffix(".tsv.partial").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    out.with_suffix(".tsv.partial").replace(out)
    print(f"Validated local MSAs and YAML inputs for {len(summary)} targets")


if __name__ == "__main__":
    main()
