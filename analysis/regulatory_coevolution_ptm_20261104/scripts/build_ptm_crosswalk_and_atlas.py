#!/usr/bin/env python3
"""Map published mouse DNMT3L PTMs to human coordinates and summarize conservation.

High-confidence coordinate mappings require exact agreement between independent
MAFFT and MACSE placental alignments. The downstream 289-ortholog atlas uses the
human-coordinate global alignment and reports residue preservation, not PTM
occupancy or modification in another species.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO


AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_alignment(path: Path):
    return {record.id: str(record.seq).upper() for record in SeqIO.parse(path, "fasta")}


def coordinate_map(mouse_aligned: str, human_aligned: str, mouse_position: int):
    mouse_coordinate = 0
    human_coordinate = 0
    for column, (mouse_residue, human_residue) in enumerate(zip(mouse_aligned, human_aligned), start=1):
        if human_residue in AA:
            human_coordinate += 1
        if mouse_residue in AA:
            mouse_coordinate += 1
        if mouse_coordinate == mouse_position and mouse_residue in AA:
            return {
                "column": column,
                "mouse_residue": mouse_residue,
                "human_position": human_coordinate if human_residue in AA else None,
                "human_residue": human_residue if human_residue in AA else "-",
            }
    raise ValueError(f"Mouse position {mouse_position} was not found")


def region(position):
    if position is None:
        return "unmapped"
    if 41 <= position <= 173:
        return "ADD"
    if 178 <= position <= 380:
        return "C-terminal"
    return "other"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--mafft", type=Path, required=True)
    parser.add_argument("--macse", type=Path, required=True)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--orthology", type=Path, required=True)
    parser.add_argument("--out-crosswalk", type=Path, required=True)
    parser.add_argument("--out-sites", type=Path, required=True)
    parser.add_argument("--out-clades", type=Path, required=True)
    args = parser.parse_args()

    sites = read_tsv(args.sites)
    alignments = {"MAFFT": load_alignment(args.mafft), "MACSE": load_alignment(args.macse)}
    atlas = load_alignment(args.atlas)
    orthology = {row["safe_id"]: row for row in read_tsv(args.orthology)}

    crosswalk = []
    for site in sites:
        mappings = {}
        for method, aln in alignments.items():
            mappings[method] = coordinate_map(
                aln["Mus_musculus"], aln["Homo_sapiens"], int(site["mouse_position"])
            )
            if mappings[method]["mouse_residue"] != site["mouse_residue"]:
                raise ValueError(
                    f"{site['published_site']} is {mappings[method]['mouse_residue']} in {method}, "
                    f"expected {site['mouse_residue']}"
                )
        same_target = (
            mappings["MAFFT"]["human_position"] == mappings["MACSE"]["human_position"]
            and mappings["MAFFT"]["human_residue"] == mappings["MACSE"]["human_residue"]
        )
        mapped = mappings["MAFFT"]["human_position"] is not None and same_target
        status = "CONCORDANT_MAPPED" if mapped else (
            "CONCORDANT_HUMAN_GAP" if same_target else "ALIGNER_DISCORDANT"
        )
        human_position = mappings["MAFFT"]["human_position"] if mapped else None
        human_residue = mappings["MAFFT"]["human_residue"] if mapped else ""
        mapped_region = region(human_position)
        if status == "CONCORDANT_MAPPED" and mapped_region in {"ADD", "C-terminal"}:
            confidence = "HIGH_CROSS_ALIGNER_STRUCTURED"
        elif status == "CONCORDANT_MAPPED":
            confidence = "LOW_DISORDERED_OR_TERMINAL"
        else:
            confidence = "UNRESOLVED"
        crosswalk.append({
            **site,
            "mafft_column": mappings["MAFFT"]["column"],
            "mafft_human_position": mappings["MAFFT"]["human_position"] or "",
            "mafft_human_residue": mappings["MAFFT"]["human_residue"],
            "macse_column": mappings["MACSE"]["column"],
            "macse_human_position": mappings["MACSE"]["human_position"] or "",
            "macse_human_residue": mappings["MACSE"]["human_residue"],
            "mapping_status": status,
            "mapping_confidence": confidence,
            "human_position": human_position or "",
            "human_residue": human_residue,
            "human_region": mapped_region,
            "same_modifiable_residue_in_human": (
                "yes" if mapped and human_residue == site["mouse_residue"] else "no" if mapped else "unknown"
            ),
        })

    crosswalk_fields = list(crosswalk[0])
    write_tsv(args.out_crosswalk, crosswalk, crosswalk_fields)

    clade_by_id = {key: row["clade"] for key, row in orthology.items()}
    missing_metadata = sorted(set(atlas) - set(clade_by_id))
    if missing_metadata:
        raise ValueError(f"Atlas records lack clade metadata: {missing_metadata[:10]}")

    site_rows = []
    clade_rows = []
    for mapped_site in crosswalk:
        if mapped_site["mapping_status"] != "CONCORDANT_MAPPED":
            continue
        pos = int(mapped_site["human_position"])
        reference_residue = mapped_site["human_residue"]
        source_residue = mapped_site["mouse_residue"]
        if "Homo_sapiens" not in atlas or atlas["Homo_sapiens"][pos - 1] != reference_residue:
            observed = atlas.get("Homo_sapiens", "")[pos - 1:pos]
            raise ValueError(
                f"Human-coordinate atlas mismatch at NP{pos}: {observed!r} != {reference_residue!r}"
            )
        per_clade = defaultdict(list)
        all_states = []
        for taxon, sequence in atlas.items():
            state = sequence[pos - 1]
            if state in AA:
                all_states.append(state)
                per_clade[clade_by_id[taxon]].append(state)

        counts = Counter(all_states)
        n = len(all_states)
        site_rows.append({
            "ptm_type": mapped_site["ptm_type"],
            "published_site": mapped_site["published_site"],
            "mouse_position": mapped_site["mouse_position"],
            "mouse_residue": source_residue,
            "human_position": pos,
            "human_residue": reference_residue,
            "human_region": mapped_site["human_region"],
            "functional_followup": mapped_site["functional_followup"],
            "resolved_n": n,
            "mouse_ptm_residue_n": counts[source_residue],
            "mouse_ptm_residue_fraction": f"{counts[source_residue] / n:.6f}" if n else "",
            "human_residue_n": counts[reference_residue],
            "human_residue_fraction": f"{counts[reference_residue] / n:.6f}" if n else "",
            "distinct_residue_states": len(counts),
            "state_counts": ";".join(f"{aa}:{counts[aa]}" for aa in sorted(counts)),
        })
        for clade in sorted(per_clade):
            clade_counts = Counter(per_clade[clade])
            clade_n = sum(clade_counts.values())
            clade_rows.append({
                "ptm_type": mapped_site["ptm_type"],
                "published_site": mapped_site["published_site"],
                "mouse_position": mapped_site["mouse_position"],
                "mouse_residue": source_residue,
                "human_position": pos,
                "human_residue": reference_residue,
                "human_region": mapped_site["human_region"],
                "clade": clade,
                "resolved_n": clade_n,
                "mouse_ptm_residue_n": clade_counts[source_residue],
                "mouse_ptm_residue_fraction": f"{clade_counts[source_residue] / clade_n:.6f}",
                "human_residue_n": clade_counts[reference_residue],
                "human_residue_fraction": f"{clade_counts[reference_residue] / clade_n:.6f}",
                "distinct_residue_states": len(clade_counts),
                "state_counts": ";".join(f"{aa}:{clade_counts[aa]}" for aa in sorted(clade_counts)),
            })

    write_tsv(args.out_sites, site_rows, list(site_rows[0]))
    write_tsv(args.out_clades, clade_rows, list(clade_rows[0]))

    status_counts = Counter(row["mapping_status"] for row in crosswalk)
    print(f"PTM source sites: {len(crosswalk)}")
    print("Mapping status:", dict(status_counts))
    print(f"Mapped sites summarized in 289-record atlas: {len(site_rows)}")


if __name__ == "__main__":
    main()
