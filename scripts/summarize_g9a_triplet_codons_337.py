#!/usr/bin/env python3
"""Recover codons underlying the DNMT3L G9a-contact triplet.

The protein-coordinate mapping is taken from the accepted local-context atlas.
Each mapped residue is cross-checked against the paired representative CDS.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from Bio import Phylo


ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "results" / "functional_site_states_by_record.tsv"
RODENTS = ROOT / "results" / "g9a_rodent_tree_states.tsv"
CDS = (
    ROOT
    / "data"
    / "processed"
    / "ncbi_337_representatives"
    / "DNMT3L_337_selection_ready_Placentalia_cds.fasta"
)
DETAIL = ROOT / "results" / "g9a_triplet_codon_states.tsv"
SUMMARY = ROOT / "results" / "g9a_muroid_codon_summary.tsv"
TREE = (
    ROOT
    / "analysis"
    / "selection_reanalysis_337"
    / "trees"
    / "Placentalia_CONSENSUS_selection_ncbi_taxonomy_topology.nwk"
)
MUROIDEA_NODE = "NCBI_337687"


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            if name is not None:
                records[name] = "".join(chunks).upper()
            name = raw[1:].split()[0]
            chunks = []
        else:
            chunks.append(raw.strip())
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


def main() -> None:
    cds_by_id = read_fasta(CDS)
    with RODENTS.open(newline="") as handle:
        rodent_rows = {r["safe_id"]: r for r in csv.DictReader(handle, delimiter="\t")}

    atlas_rows: dict[str, dict[str, dict[str, str]]] = {}
    with ATLAS.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (
                row["feature_id"] in {
                    "G9A_INTERFACE_L139",
                    "G9A_INTERFACE_F140",
                    "G9A_INTERFACE_I141",
                }
                and row["broad_lineage"] == "Placentalia"
                and row["mapping_status"] == "PASS"
            ):
                atlas_rows.setdefault(row["safe_id"], {})[row["feature_id"]] = row

    detail_rows: list[dict[str, str]] = []
    order = ["G9A_INTERFACE_L139", "G9A_INTERFACE_F140", "G9A_INTERFACE_I141"]
    for safe_id, sites in sorted(atlas_rows.items()):
        if not all(site in sites for site in order):
            continue
        if safe_id not in cds_by_id:
            raise KeyError(f"No paired CDS for {safe_id}")
        seq = cds_by_id[safe_id]
        codons: list[str] = []
        residues: list[str] = []
        for site in order:
            row = sites[site]
            pos = int(float(row["query_position"]))
            codon = seq[(pos - 1) * 3 : pos * 3]
            residue = CODON_TABLE.get(codon, "?")
            if residue != row["query_residue"]:
                raise ValueError(
                    f"CDS/protein mismatch for {safe_id} {site}: "
                    f"{codon}->{residue}, expected {row['query_residue']}"
                )
            codons.append(codon)
            residues.append(residue)
        rodent = rodent_rows.get(safe_id, {})
        detail_rows.append(
            {
                "safe_id": safe_id,
                "species": sites[order[0]]["species"],
                "triplet_state": "".join(residues),
                "triplet_codons": "-".join(codons),
                "middle_residue": residues[1],
                "middle_codon": codons[1],
                "rodent_family_group": rodent.get("family", ""),
                "muroid": rodent.get("muroid", ""),
            }
        )

    fields = list(detail_rows[0])
    with DETAIL.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(detail_rows)

    muroids = [r for r in detail_rows if r["muroid"] == "True"]
    nonmuroid_rodents = [
        r for r in detail_rows if r["muroid"] == "False" and r["rodent_family_group"]
    ]
    jaculus = [r for r in nonmuroid_rodents if r["safe_id"] == "Jaculus_jaculus"]
    if len(muroids) != 27 or len(jaculus) != 1:
        raise ValueError(f"Expected 27 muroids and one Jaculus; got {len(muroids)}, {len(jaculus)}")

    muroid_codons = Counter(r["middle_codon"] for r in muroids)
    other_rodent_codons = Counter(r["middle_codon"] for r in nonmuroid_rodents)
    jaculus_codon = jaculus[0]["middle_codon"]
    dominant_muroid_codon = muroid_codons.most_common(1)[0][0]
    dominant_change_count = sum(
        a != b for a, b in zip(jaculus_codon, dominant_muroid_codon)
    )

    muroid_state = {row["safe_id"]: row["middle_codon"] for row in muroids}
    tree = Phylo.read(TREE, "newick")
    muroid_clade = next(clade for clade in tree.find_clades() if clade.name == MUROIDEA_NODE)
    ttt_ids = {safe_id for safe_id, codon in muroid_state.items() if codon == "TTT"}
    ttt_mrca_descendants = {
        terminal.name for terminal in tree.common_ancestor(ttt_ids).get_terminals()
    }

    def fitch(clade) -> tuple[set[str], int]:
        if clade.is_terminal():
            return {muroid_state[clade.name]}, 0
        children = [fitch(child) for child in clade.clades]
        states = [child[0] for child in children]
        score = sum(child[1] for child in children)
        intersection = set.intersection(*states)
        if intersection:
            return intersection, score
        return set.union(*states), score + 1

    _, minimum_muroid_codon_changes = fitch(muroid_clade)
    summary_rows = [
        {"metric": "muroid_records", "value": str(len(muroids))},
        {"metric": "muroid_middle_residue", "value": "F"},
        {"metric": "muroid_middle_codon_counts", "value": ";".join(f"{k}:{v}" for k, v in sorted(muroid_codons.items()))},
        {"metric": "nonmuroid_rodent_records", "value": str(len(nonmuroid_rodents))},
        {"metric": "nonmuroid_rodent_middle_codon_counts", "value": ";".join(f"{k}:{v}" for k, v in sorted(other_rodent_codons.items()))},
        {"metric": "jaculus_middle_codon", "value": jaculus_codon},
        {"metric": "same_middle_codon_in_all_muroids", "value": str(len(muroid_codons) == 1)},
        {"metric": "dominant_muroid_middle_codon", "value": dominant_muroid_codon},
        {"metric": "nucleotide_differences_jaculus_to_dominant_muroid_codon", "value": str(dominant_change_count)},
        {"metric": "muroid_F_codons_differ_synonymously", "value": str(set(muroid_codons) == {"TTC", "TTT"})},
        {"metric": "ttt_muroid_species", "value": ";".join(sorted(ttt_ids))},
        {"metric": "ttt_muroids_form_single_clade", "value": str(ttt_mrca_descendants == ttt_ids)},
        {"metric": "minimum_synonymous_codon_changes_within_muroids", "value": str(minimum_muroid_codon_changes)},
    ]
    with SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote {DETAIL} ({len(detail_rows)} accepted eutherian mappings)")
    print(f"Wrote {SUMMARY}")
    for row in summary_rows:
        print(f"{row['metric']}: {row['value']}")


if __name__ == "__main__":
    main()
