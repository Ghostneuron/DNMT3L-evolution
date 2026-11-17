#!/usr/bin/env python3
"""Map the experimentally tested mouse DNMT3L G9a-contact motif.

The L139/F140/I141 motif was nominated by AlphaFold-Multimer and tested by
triple-alanine mutagenesis and co-immunoprecipitation (PMID 40972525). This
script requires MAFFT/MACSE agreement for the mouse-to-human coordinate map,
then reports observed residue triplets in the human-coordinate orthologue
atlas. Conservation supports retention of the tested local chemistry; it does
not demonstrate G9a binding in another species.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO


AA = set("ACDEFGHIKLMNPQRSTVWY")
MOUSE_SITES = ((139, "L"), (140, "F"), (141, "I"))


def load_alignment(path: Path):
    return {record.id: str(record.seq).upper() for record in SeqIO.parse(path, "fasta")}


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def coordinate_map(mouse_aligned: str, human_aligned: str, mouse_position: int):
    mouse_coordinate = 0
    human_coordinate = 0
    for column, (mouse_residue, human_residue) in enumerate(zip(mouse_aligned, human_aligned), start=1):
        if human_residue in AA:
            human_coordinate += 1
        if mouse_residue in AA:
            mouse_coordinate += 1
        if mouse_coordinate == mouse_position and mouse_residue in AA:
            return column, mouse_residue, human_coordinate if human_residue in AA else None, human_residue
    raise ValueError(f"Mouse position {mouse_position} was not found")


def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mafft", type=Path, required=True)
    parser.add_argument("--macse", type=Path, required=True)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--orthology", type=Path, required=True)
    parser.add_argument("--out-crosswalk", type=Path, required=True)
    parser.add_argument("--out-species", type=Path, required=True)
    parser.add_argument("--out-clades", type=Path, required=True)
    args = parser.parse_args()

    alignments = {"MAFFT": load_alignment(args.mafft), "MACSE": load_alignment(args.macse)}
    mapped = []
    for mouse_position, expected in MOUSE_SITES:
        calls = {}
        for method, aln in alignments.items():
            calls[method] = coordinate_map(aln["Mus_musculus"], aln["Homo_sapiens"], mouse_position)
            if calls[method][1] != expected:
                raise ValueError(f"Mouse {mouse_position}: {calls[method][1]} != {expected} in {method}")
        if calls["MAFFT"][2:] != calls["MACSE"][2:]:
            raise ValueError(f"Aligners disagree at mouse {expected}{mouse_position}: {calls}")
        if calls["MAFFT"][2] is None:
            raise ValueError(f"Mouse {expected}{mouse_position} maps to a human gap")
        mapped.append((mouse_position, expected, calls))

    crosswalk_rows = []
    human_positions = []
    for mouse_position, expected, calls in mapped:
        human_position = calls["MAFFT"][2]
        human_positions.append(human_position)
        crosswalk_rows.append({
            "mouse_position": mouse_position,
            "mouse_residue": expected,
            "mafft_column": calls["MAFFT"][0],
            "macse_column": calls["MACSE"][0],
            "human_position": human_position,
            "human_residue": calls["MAFFT"][3],
            "mapping_status": "CONCORDANT_MAPPED",
            "functional_evidence": "triple-alanine mutant reduced G9a co-IP and chromatin association",
            "claim_boundary": "residue conservation does not establish G9a binding",
        })
    write_tsv(args.out_crosswalk, crosswalk_rows, list(crosswalk_rows[0]))

    atlas = load_alignment(args.atlas)
    metadata = {row["safe_id"]: row for row in read_tsv(args.orthology)}
    if set(atlas) - set(metadata):
        raise ValueError("Atlas contains records without orthology metadata")

    source_triplet = "".join(expected for _, expected in MOUSE_SITES)
    species_rows = []
    clade_triplets = defaultdict(list)
    for safe_id, sequence in sorted(atlas.items()):
        states = [sequence[position - 1] for position in human_positions]
        resolved = all(state in AA for state in states)
        triplet = "".join(states) if resolved else "UNRESOLVED"
        clade = metadata[safe_id]["clade"]
        if resolved:
            clade_triplets[clade].append(triplet)
        species_rows.append({
            "safe_id": safe_id,
            "species": metadata[safe_id].get("species", safe_id.replace("_", " ")),
            "clade": clade,
            "human_positions": ";".join(map(str, human_positions)),
            "states": ";".join(states),
            "triplet": triplet,
            "mouse_LFI_retained": "yes" if triplet == source_triplet else "no" if resolved else "unknown",
        })
    write_tsv(args.out_species, species_rows, list(species_rows[0]))

    clade_rows = []
    for clade in sorted(clade_triplets):
        counts = Counter(clade_triplets[clade])
        n = sum(counts.values())
        clade_rows.append({
            "clade": clade,
            "resolved_n": n,
            "mouse_LFI_n": counts[source_triplet],
            "mouse_LFI_fraction": f"{counts[source_triplet] / n:.6f}" if n else "",
            "distinct_triplets": len(counts),
            "triplet_counts": ";".join(f"{triplet}:{counts[triplet]}" for triplet in sorted(counts)),
        })
    write_tsv(args.out_clades, clade_rows, list(clade_rows[0]))

    overall = Counter(row["triplet"] for row in species_rows if row["triplet"] != "UNRESOLVED")
    print("Mouse-to-human positions:", human_positions)
    print("Resolved orthologues:", sum(overall.values()))
    print("Mouse LFI retained:", overall[source_triplet])
    print("Triplets:", dict(overall))


if __name__ == "__main__":
    main()
