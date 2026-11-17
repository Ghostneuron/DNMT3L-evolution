#!/usr/bin/env python3
"""Test rodent conservation at predicted G9a residues contacting DNMT3L 139-141.

Contacts are defined from rank 1 of ModelArchive ma-dm-mmk9-101 using a
5-A heavy-atom cutoff. Orthologue residues are mapped with a local alignment
of the 25-residue mouse G9a context surrounding the contact set.
"""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1


ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "data" / "raw" / "modelarchive_ma-dm-mmk9-101" / "ma-dm-mmk9-101.cif"
PROTEINS = ROOT / "data" / "raw" / "ncbi_ehmt2_orthologs" / "protein.faa"
RODENTS = ROOT / "results" / "g9a_rodent_tree_states.tsv"
DETAIL = ROOT / "results" / "g9a_partner_contact_states.tsv"
SUMMARY = ROOT / "results" / "g9a_partner_contact_summary.tsv"


def read_fasta(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    header = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                match = re.search(r"\[organism=([^]]+)\]", header)
                if match:
                    rows.append((header.split()[0], match.group(1), "".join(chunks)))
            header = line[1:]
            chunks = []
        else:
            chunks.append(line.strip())
    if header is not None:
        match = re.search(r"\[organism=([^]]+)\]", header)
        if match:
            rows.append((header.split()[0], match.group(1), "".join(chunks)))
    return rows


def model_contact_definition() -> tuple[str, list[int], list[int], str]:
    structure = MMCIFParser(QUIET=True).get_structure("ma-dm-mmk9-101", MODEL)
    model = structure[0]
    g9a = model["A"]
    dnmt3l_atoms = [
        atom
        for residue in model["B"]
        if residue.id[1] in (139, 140, 141)
        for atom in residue
        if atom.element != "H"
    ]
    contacts: list[int] = []
    for residue in g9a:
        distance = min((atom - other for atom in residue for other in dnmt3l_atoms), default=999.0)
        if distance <= 5.0:
            contacts.append(residue.id[1])
    if contacts != [699, 700, 701, 702, 703]:
        raise ValueError(f"Unexpected rank-1 contact set: {contacts}")
    context_positions = list(range(689, 714))
    context = "".join(seq1(g9a[(" ", pos, " ")].resname) for pos in context_positions)
    contact_state = "".join(seq1(g9a[(" ", pos, " ")].resname) for pos in contacts)
    return context, context_positions, contacts, contact_state


def mapped_state(
    reference: str,
    reference_positions: list[int],
    contact_positions: list[int],
    query: str,
    aligner: Align.PairwiseAligner,
) -> tuple[float, float, str, str]:
    alignment = aligner.align(reference, query)[0]
    indices = alignment.indices
    mapping: dict[int, int] = {}
    for column in range(indices.shape[1]):
        ref_idx = int(indices[0, column])
        qry_idx = int(indices[1, column])
        if ref_idx >= 0 and qry_idx >= 0:
            mapping[reference_positions[ref_idx]] = qry_idx
    if not all(pos in mapping for pos in contact_positions):
        return alignment.score, 0.0, "", ""
    state = "".join(query[mapping[pos]] for pos in contact_positions)
    aligned_ref_positions = [pos for pos in reference_positions if pos in mapping]
    identity = sum(
        reference[reference_positions.index(pos)] == query[mapping[pos]]
        for pos in aligned_ref_positions
    ) / len(reference_positions)
    local = "".join(query[mapping[pos]] if pos in mapping else "-" for pos in reference_positions)
    return alignment.score, identity, state, local


def main() -> None:
    reference, reference_positions, contacts, reference_state = model_contact_definition()
    with RODENTS.open(newline="") as handle:
        rodent_rows = list(csv.DictReader(handle, delimiter="\t"))
    target_species = {r["species"] for r in rodent_rows}

    candidates: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for accession, species, sequence in read_fasta(PROTEINS):
        if species in target_species:
            candidates[species].append((accession, sequence))

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    detail_rows: list[dict[str, str]] = []
    for rodent in rodent_rows:
        species = rodent["species"]
        ranked: list[tuple[float, float, int, str, str, str]] = []
        for accession, sequence in candidates.get(species, []):
            score, identity, state, local = mapped_state(
                reference, reference_positions, contacts, sequence, aligner
            )
            ranked.append((score, identity, len(sequence), accession, state, local))
        if not ranked:
            detail_rows.append(
                {
                    "safe_id": rodent["safe_id"],
                    "species": species,
                    "muroid": rodent["muroid"],
                    "dnmt3l_triplet": rodent["triplet_state"],
                    "ehmt2_accession": "",
                    "ehmt2_length": "",
                    "local_identity_25aa": "",
                    "g9a_contact_state_699_703": "not available",
                    "g9a_local_context_689_713": "",
                }
            )
            continue
        score, identity, length, accession, state, local = max(ranked)
        detail_rows.append(
            {
                "safe_id": rodent["safe_id"],
                "species": species,
                "muroid": rodent["muroid"],
                "dnmt3l_triplet": rodent["triplet_state"],
                "ehmt2_accession": accession,
                "ehmt2_length": str(length),
                "local_identity_25aa": f"{identity:.3f}",
                "g9a_contact_state_699_703": state,
                "g9a_local_context_689_713": local,
            }
        )

    with DETAIL.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(detail_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(detail_rows)

    available = [r for r in detail_rows if r["g9a_contact_state_699_703"] != "not available"]
    muroid = [r for r in available if r["muroid"] == "True"]
    nonmuroid = [r for r in available if r["muroid"] == "False"]
    summary_rows = [
        {"metric": "modelarchive_accession", "value": "ma-dm-mmk9-101"},
        {"metric": "model_rank", "value": "1"},
        {"metric": "heavy_atom_cutoff_angstrom", "value": "5.0"},
        {"metric": "dnmt3l_residues", "value": "139-141"},
        {"metric": "g9a_contact_positions", "value": ",".join(map(str, contacts))},
        {"metric": "mouse_g9a_contact_state", "value": reference_state},
        {"metric": "rodents_with_ehmt2_mapping", "value": str(len(available))},
        {"metric": "muroids_with_ehmt2_mapping", "value": str(len(muroid))},
        {"metric": "muroid_contact_state_counts", "value": ";".join(f"{k}:{v}" for k, v in sorted(Counter(r["g9a_contact_state_699_703"] for r in muroid).items()))},
        {"metric": "nonmuroid_rodents_with_ehmt2_mapping", "value": str(len(nonmuroid))},
        {"metric": "nonmuroid_contact_state_counts", "value": ";".join(f"{k}:{v}" for k, v in sorted(Counter(r["g9a_contact_state_699_703"] for r in nonmuroid).items()))},
    ]
    with SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote {DETAIL}")
    print(f"Wrote {SUMMARY}")
    for row in summary_rows:
        print(f"{row['metric']}: {row['value']}")


if __name__ == "__main__":
    main()
