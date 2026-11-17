#!/usr/bin/env python3
"""Test conservation of predicted G9a contact residues across eutherians.

The target set is the 219 eutherian DNMT3L records with accepted local-context
mappings for the NP105/NP106/NP107 triplet. G9a contacts are defined from rank
1 of ModelArchive ma-dm-mmk9-101 using a 5-A heavy-atom cutoff. EHMT2
orthologues are mapped with a local alignment of the 25-residue mouse context.
"""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1


ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "data" / "raw" / "modelarchive_ma-dm-mmk9-101" / "ma-dm-mmk9-101.cif"
PROTEINS = ROOT / "data" / "raw" / "ncbi_ehmt2_orthologs" / "protein.faa"
TRIPLETS = ROOT / "results" / "g9a_triplet_codon_states.tsv"
DETAIL = ROOT / "results" / "g9a_partner_contact_states_eutherian.tsv"
SUMMARY = ROOT / "results" / "g9a_partner_contact_summary_eutherian.tsv"


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
    aligned = [pos for pos in reference_positions if pos in mapping]
    identity = sum(
        reference[reference_positions.index(pos)] == query[mapping[pos]] for pos in aligned
    ) / len(reference_positions)
    local = "".join(query[mapping[pos]] if pos in mapping else "-" for pos in reference_positions)
    return alignment.score, identity, state, local


def main() -> None:
    reference, reference_positions, contacts, reference_state = model_contact_definition()
    with TRIPLETS.open(newline="") as handle:
        target_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(target_rows) != 219:
        raise ValueError(f"Expected 219 accepted eutherian triplets, found {len(target_rows)}")

    target_species = {row["species"] for row in target_rows}
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
    for target in target_rows:
        species = target["species"]
        ranked: list[tuple[float, float, int, str, str, str]] = []
        for accession, sequence in candidates.get(species, []):
            score, identity, state, local = mapped_state(
                reference, reference_positions, contacts, sequence, aligner
            )
            ranked.append((score, identity, len(sequence), accession, state, local))
        if ranked:
            score, identity, length, accession, state, local = max(ranked)
        else:
            identity, length, accession, state, local = 0.0, 0, "", "not available", ""
        detail_rows.append(
            {
                "safe_id": target["safe_id"],
                "species": species,
                "dnmt3l_triplet": target["triplet_state"],
                "muroid": target["muroid"],
                "ehmt2_accession": accession,
                "ehmt2_length": str(length) if accession else "",
                "local_identity_25aa": f"{identity:.3f}" if accession else "",
                "g9a_contact_state_699_703": state,
                "g9a_local_context_689_713": local,
            }
        )

    with DETAIL.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(detail_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(detail_rows)

    available = [row for row in detail_rows if row["g9a_contact_state_699_703"] != "not available"]
    missing = [row for row in detail_rows if row["g9a_contact_state_699_703"] == "not available"]
    states = Counter(row["g9a_contact_state_699_703"] for row in available)
    muroids = [row for row in available if row["muroid"] == "True"]
    other = [row for row in available if row["muroid"] != "True"]
    summary_rows = [
        {"metric": "modelarchive_accession", "value": "ma-dm-mmk9-101"},
        {"metric": "model_rank", "value": "1"},
        {"metric": "heavy_atom_cutoff_angstrom", "value": "5.0"},
        {"metric": "dnmt3l_residues", "value": "139-141"},
        {"metric": "g9a_contact_positions", "value": ",".join(map(str, contacts))},
        {"metric": "mouse_g9a_contact_state", "value": reference_state},
        {"metric": "accepted_eutherian_dnmt3l_triplets", "value": str(len(target_rows))},
        {"metric": "eutherians_with_ehmt2_mapping", "value": str(len(available))},
        {"metric": "eutherians_without_ehmt2_mapping", "value": str(len(missing))},
        {"metric": "eutherian_contact_state_counts", "value": ";".join(f"{k}:{v}" for k, v in sorted(states.items()))},
        {"metric": "muroids_with_ehmt2_mapping", "value": str(len(muroids))},
        {"metric": "muroid_contact_state_counts", "value": ";".join(f"{k}:{v}" for k, v in sorted(Counter(row["g9a_contact_state_699_703"] for row in muroids).items()))},
        {"metric": "other_eutherians_with_ehmt2_mapping", "value": str(len(other))},
        {"metric": "other_eutherian_contact_state_counts", "value": ";".join(f"{k}:{v}" for k, v in sorted(Counter(row["g9a_contact_state_699_703"] for row in other).items()))},
        {"metric": "missing_species", "value": ";".join(row["species"] for row in missing)},
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
