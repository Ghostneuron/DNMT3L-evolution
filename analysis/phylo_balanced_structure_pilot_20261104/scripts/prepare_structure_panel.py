#!/usr/bin/env python3
"""Prepare the frozen DNMT3L monomer and fixed-partner complex pilot.

Extant monomers use the curated full protein sequence. Complexes use the
human-coordinate DNMT3L C-terminal interval (NP 178-380) and a constant human
DNMT3A1 catalytic-domain construct (residues 623-912). Amniota and Theria use
maximum-posterior IQ-TREE sequences while preserving per-site uncertainty.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

from Bio import SeqIO


AA = set("ACDEFGHIKLMNPQRSTVWY")
ANCESTORS = {"ANCESTRAL_AMNIOTA": "NCBI_32524", "ANCESTRAL_THERIA": "NCBI_32525"}
DNMT3L_CTD = (178, 380)
DNMT3A_CONSTRUCT = (623, 912)


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_text(path: Path, content: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    partial.write_text(content)
    partial.replace(path)


def write_tsv(path: Path, rows, fields=None):
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def validate_sequence(name: str, sequence: str, low: int, high: int):
    invalid = sorted(set(sequence) - AA)
    if invalid:
        raise ValueError(f"{name} contains invalid residues: {invalid}")
    if not low <= len(sequence) <= high:
        raise ValueError(f"{name} length {len(sequence)} outside {low}-{high}")


def fasta(name: str, sequence: str) -> str:
    wrapped = "\n".join(sequence[i : i + 80] for i in range(0, len(sequence), 80))
    return f">{name}\n{wrapped}\n"


def monomer_yaml(name: str, sequence: str) -> str:
    return (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        f"      sequence: {sequence}\n"
    )


def complex_yaml(name: str, dnmt3a: str, dnmt3l: str) -> str:
    return (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        f"      sequence: {dnmt3a}\n"
        "  - protein:\n"
        "      id: B\n"
        f"      sequence: {dnmt3l}\n"
    )


def asr_rows(path: Path, wanted_nodes: set[str]):
    with path.open() as handle:
        data_lines = (line for line in handle if not line.startswith("#"))
        rows = [row for row in csv.DictReader(data_lines, delimiter="\t") if row["Node"] in wanted_nodes]
    by_node = {node: [] for node in wanted_nodes}
    posterior_rows = []
    for row in rows:
        states = sorted(
            ((key.removeprefix("p_"), float(value)) for key, value in row.items() if key.startswith("p_")),
            key=lambda item: item[1],
            reverse=True,
        )
        top_state, top_pp = states[0]
        second_state, second_pp = states[1]
        if top_state != row["State"]:
            raise ValueError(f"IQ-TREE state mismatch at {row['Node']} site {row['Site']}")
        by_node[row["Node"]].append((int(row["Site"]), top_state, top_pp))
        posterior_rows.append({
            "node": row["Node"],
            "human_position": row["Site"],
            "maximum_posterior_state": top_state,
            "maximum_posterior_probability": f"{top_pp:.5f}",
            "second_state": second_state,
            "second_probability": f"{second_pp:.5f}",
            "resolved_at_0.80": "yes" if top_pp >= 0.80 else "no",
        })
    sequences = {}
    confidence = {}
    for node, calls in by_node.items():
        calls.sort()
        if [position for position, _, _ in calls] != list(range(1, 388)):
            raise ValueError(f"{node} does not contain exactly sites 1-387")
        sequences[node] = "".join(state for _, state, _ in calls)
        confidence[node] = {position: pp for position, _, pp in calls}
    return sequences, confidence, posterior_rows


def subsequence_positions(projected: str, full: str):
    """Map an ungapped human-coordinate projection into its curated full sequence."""
    positions = []
    cursor = 0
    for residue in projected:
        found = full.find(residue, cursor)
        if found < 0:
            raise ValueError("Projected sequence is not an ordered subsequence of curated full sequence")
        positions.append(found + 1)
        cursor = found + 1
    return positions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--dnmt3l-matrix", type=Path, required=True)
    parser.add_argument("--dnmt3l-projection", type=Path, required=True)
    parser.add_argument("--matched-inventory", type=Path, required=True)
    parser.add_argument("--asr-state", type=Path, required=True)
    parser.add_argument("--dnmt3a-proteins", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()

    panel = read_tsv(args.panel)
    matrix = {row["safe_id"]: row for row in read_tsv(args.dnmt3l_matrix)}
    projection = {record.id: str(record.seq).upper() for record in SeqIO.parse(args.dnmt3l_projection, "fasta")}
    matched = {row["safe_id"] for row in read_tsv(args.matched_inventory)}
    dnmt3a = {record.id: str(record.seq).upper() for record in SeqIO.parse(args.dnmt3a_proteins, "fasta")}

    human_dnmt3a = dnmt3a["Homo_sapiens"]
    dnmt3a_construct = human_dnmt3a[DNMT3A_CONSTRUCT[0] - 1 : DNMT3A_CONSTRUCT[1]]
    validate_sequence("human DNMT3A1 623-912", dnmt3a_construct, 290, 290)
    expected = {725: "E", 726: "G", 729: "R", 732: "F", 772: "F"}
    for position, residue in expected.items():
        if human_dnmt3a[position - 1] != residue:
            raise ValueError(f"DNMT3A reference mismatch at {position}")

    ancestor_sequences, ancestor_confidence, posterior_rows = asr_rows(args.asr_state, set(ANCESTORS.values()))
    write_tsv(args.output_root / "metadata/ancestral_site_posteriors.tsv", posterior_rows)

    manifest = []
    crosswalk = []
    hashes = []
    seen = set()
    for choice in panel:
        safe_id = choice["safe_id"]
        if safe_id in seen:
            raise ValueError(f"Duplicate panel entry: {safe_id}")
        seen.add(safe_id)

        if safe_id in ANCESTORS:
            node = ANCESTORS[safe_id]
            full = ancestor_sequences[node]
            ctd = full[DNMT3L_CTD[0] - 1 : DNMT3L_CTD[1]]
            accession = f"IQTREE:{node}"
            source = "maximum-posterior ancestral reconstruction"
            resolved_full = sum(pp >= 0.80 for pp in ancestor_confidence[node].values())
            ctd_pp = [ancestor_confidence[node][p] for p in range(DNMT3L_CTD[0], DNMT3L_CTD[1] + 1)]
            resolved_ctd = sum(pp >= 0.80 for pp in ctd_pp)
            for position in range(1, 388):
                crosswalk.append({
                    "safe_id": safe_id,
                    "human_position": position,
                    "monomer_input_position": position,
                    "monomer_residue": full[position - 1],
                    "complex_chain_B_position": position - DNMT3L_CTD[0] + 1 if DNMT3L_CTD[0] <= position <= DNMT3L_CTD[1] else "",
                    "ancestral_posterior_probability": f"{ancestor_confidence[node][position]:.5f}",
                })
        else:
            if safe_id not in matrix or matrix[safe_id].get("Selection_Ready") != "True":
                raise ValueError(f"Panel entry is not selection-ready DNMT3L: {safe_id}")
            if choice["complex_eligible"] == "yes" and safe_id not in matched:
                raise ValueError(f"Complex-eligible entry lacks matched DNMT3A inventory record: {safe_id}")
            full = matrix[safe_id]["protein_seq"].upper()
            projected = projection[safe_id]
            if len(projected) != 387:
                raise ValueError(f"Projection length mismatch for {safe_id}: {len(projected)}")
            ungapped = "".join(residue for residue in projected if residue in AA)
            input_positions = subsequence_positions(ungapped, full)
            human_to_input = {}
            index = 0
            for human_position, residue in enumerate(projected, start=1):
                if residue in AA:
                    human_to_input[human_position] = input_positions[index]
                    index += 1
            ctd_positions = [p for p in range(DNMT3L_CTD[0], DNMT3L_CTD[1] + 1) if p in human_to_input]
            ctd = "".join(projected[p - 1] for p in ctd_positions)
            accession = matrix[safe_id]["protein_accession"]
            source = "curated extant orthologue"
            resolved_full = len(human_to_input)
            resolved_ctd = len(ctd_positions)
            ctd_input_index = {human_position: i for i, human_position in enumerate(ctd_positions, start=1)}
            for human_position in range(1, 388):
                if human_position not in human_to_input:
                    continue
                crosswalk.append({
                    "safe_id": safe_id,
                    "human_position": human_position,
                    "monomer_input_position": human_to_input[human_position],
                    "monomer_residue": projected[human_position - 1],
                    "complex_chain_B_position": ctd_input_index.get(human_position, ""),
                    "ancestral_posterior_probability": "",
                })

        # Curated orthologues include lineage-specific N-terminal extensions;
        # Geotrypetes and dog are deliberately retained rather than truncated.
        validate_sequence(f"{safe_id} full DNMT3L", full, 300, 700)
        validate_sequence(f"{safe_id} DNMT3L CTD", ctd, 175, 210)

        atomic_text(args.output_root / f"inputs/monomer_fasta/{safe_id}.fasta", fasta(safe_id, full))
        atomic_text(args.output_root / f"inputs/monomer_yaml/{safe_id}.yaml", monomer_yaml(safe_id, full))
        atomic_text(args.output_root / f"inputs/complex_yaml/{safe_id}.yaml", complex_yaml(safe_id, dnmt3a_construct, ctd))

        full_hash = sha256_text(full)
        ctd_hash = sha256_text(ctd)
        manifest.append({
            **choice,
            "source": source,
            "protein_accession_or_node": accession,
            "monomer_length": len(full),
            "homologous_human_positions_resolved": resolved_full,
            "ctd_length": len(ctd),
            "ctd_human_positions_resolved": resolved_ctd,
            "ancestral_sites_pp_ge_0.80": resolved_full if safe_id in ANCESTORS else "",
            "ancestral_ctd_sites_pp_ge_0.80": resolved_ctd if safe_id in ANCESTORS else "",
            "monomer_sequence_sha256": full_hash,
            "ctd_sequence_sha256": ctd_hash,
            "monomer_yaml": f"inputs/monomer_yaml/{safe_id}.yaml",
            "complex_yaml": f"inputs/complex_yaml/{safe_id}.yaml",
        })
        hashes.extend([
            {"safe_id": safe_id, "artifact": "DNMT3L_full", "length": len(full), "sha256": full_hash},
            {"safe_id": safe_id, "artifact": "DNMT3L_CTD_human178_380", "length": len(ctd), "sha256": ctd_hash},
        ])

    atomic_text(args.output_root / "inputs/dnmt3a_human_623_912.fasta", fasta("Homo_sapiens_DNMT3A1_623_912", dnmt3a_construct))
    hashes.append({
        "safe_id": "Homo_sapiens",
        "artifact": "DNMT3A1_623_912_fixed_partner",
        "length": len(dnmt3a_construct),
        "sha256": sha256_text(dnmt3a_construct),
    })
    write_tsv(args.output_root / "metadata/panel_manifest.tsv", manifest)
    write_tsv(args.output_root / "metadata/human_coordinate_crosswalk.tsv", crosswalk)
    write_tsv(args.output_root / "metadata/sequence_hashes.tsv", hashes)

    if len(manifest) != 24:
        raise ValueError(f"Expected 24 panel inputs, observed {len(manifest)}")
    print(f"Validated panel inputs: {len(manifest)}")
    print(f"Extant: {len(manifest) - len(ANCESTORS)}; ancestral: {len(ANCESTORS)}")
    print(f"Fixed DNMT3A1 partner length: {len(dnmt3a_construct)}")
    for ancestral_id, node in ANCESTORS.items():
        total = sum(pp >= 0.80 for pp in ancestor_confidence[node].values())
        ctd = sum(ancestor_confidence[node][p] >= 0.80 for p in range(178, 381))
        print(f"{ancestral_id}: PP>=0.80 full {total}/387; CTD {ctd}/203")


if __name__ == "__main__":
    main()
