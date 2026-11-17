#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "analysis/switching_helix_asr_337"
INPUT = ROOT / "inputs/DNMT3L_289_human_coordinate_strict_focal_masked_codon.fasta"
OUTPUT = ROOT / "inputs/DNMT3L_286_human_coordinate_strict_focal_masked_codon_deduplicated.fasta"
MANIFEST = PROJECT / "results/switching_helix_337_selection_exact_duplicate_manifest.tsv"
TREE_NAMES = [
    "DNMT3L_289_taxonomy_therian_stem_foreground.nwk",
    "DNMT3L_289_gene_tree_therian_stem_foreground.nwk",
    "DNMT3L_289_taxonomy_theria_test_sauropsida_reference.nwk",
    "DNMT3L_289_gene_tree_theria_test_sauropsida_reference.nwk",
]


def main() -> None:
    records = list(SeqIO.parse(INPUT, "fasta"))
    by_sequence: dict[str, list[str]] = defaultdict(list)
    record_by_name = {record.id: record for record in records}
    for record in records:
        by_sequence[str(record.seq).upper()].append(record.id)

    keep: set[str] = set()
    rows: list[dict[str, object]] = []
    for group_index, names in enumerate(sorted(by_sequence.values(), key=lambda values: sorted(values)), start=1):
        names = sorted(names)
        representative = "Homo_sapiens" if "Homo_sapiens" in names else names[0]
        keep.add(representative)
        for name in names:
            rows.append(
                {
                    "duplicate_group": group_index,
                    "group_size": len(names),
                    "safe_id": name,
                    "representative": representative,
                    "retained": name == representative,
                    "rule": "retain human when present, otherwise lexicographically first safe_id",
                }
            )
    dropped = set(record_by_name) - keep
    if len(dropped) != 3 or len(keep) != 286:
        raise ValueError(f"Expected 3 dropped and 286 retained, observed {len(dropped)} and {len(keep)}")

    output_records = []
    for name in sorted(keep):
        source = record_by_name[name]
        output_records.append(SeqRecord(Seq(str(source.seq)), id=name, description=""))
    with OUTPUT.open("w") as handle:
        SeqIO.write(output_records, handle, "fasta")
    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    for name in TREE_NAMES:
        source = ROOT / "trees" / name
        target = ROOT / "trees" / name.replace("289", "286").replace(".nwk", "_deduplicated.nwk")
        tree = Phylo.read(source, "newick")
        for tip in list(tree.get_terminals()):
            clean_name = (tip.name or "").replace("{Foreground}", "").replace("{Test}", "").replace("{Reference}", "")
            if clean_name in dropped:
                tree.prune(tip)
        observed = {
            (tip.name or "").replace("{Foreground}", "").replace("{Test}", "").replace("{Reference}", "")
            for tip in tree.get_terminals()
        }
        if observed != keep:
            raise ValueError(f"Pruned tree mismatch for {name}")
        Phylo.write(tree, target, "newick")
        print(target.relative_to(PROJECT))
    print(f"Dropped exact duplicates: {sorted(dropped)}")
    print(OUTPUT.relative_to(PROJECT))


if __name__ == "__main__":
    main()
