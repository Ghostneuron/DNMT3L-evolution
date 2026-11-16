#!/usr/bin/env python3
"""Leave-one-amphibian sensitivity for the NP349 ancestral reconstruction."""

from __future__ import annotations

import copy
import csv
import os
import subprocess
from pathlib import Path

from Bio import Phylo, SeqIO


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "analysis/switching_helix_asr_337"
ALIGNMENT = ROOT / "inputs/DNMT3L_289_human_coordinate_strict_focal_masked.fasta"
ATLAS = PROJECT / "results/functional_site_states_by_record.tsv"
OUT = ROOT / "runs/NP349_LEAVE_ONE_AMPHIBIAN"
RESULT = PROJECT / "results/switching_helix_337_np349_leave_one_amphibian.tsv"
AMPHIBIANS = (
    "Geotrypetes_seraphini",
    "Lissotriton_helveticus",
    "Lissotriton_vulgaris",
    "Triturus_cristatus",
)
TOPOLOGIES = {
    "taxonomy": ROOT / "trees/DNMT3L_289_ncbi_taxonomy_topology.nwk",
    "gene_tree": ROOT / "trees/DNMT3L_289_gene_tree_topology.nwk",
}
IQTREE2_BIN = os.environ.get("IQTREE2_BIN", "iqtree2")


def read_lineages() -> dict[str, set[str]]:
    by_tip: dict[str, str] = {}
    with ATLAS.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_tip.setdefault(row["safe_id"], row["broad_lineage"])
    lineages: dict[str, set[str]] = {}
    for tip, lineage in by_tip.items():
        lineages.setdefault(lineage, set()).add(tip)
    lineages["Theria"] = lineages["Marsupialia"] | lineages["Placentalia"]
    lineages["Amniota"] = lineages["Theria"] | lineages["Sauropsida"]
    return lineages


def read_np349_state(path: Path) -> dict[str, tuple[str, float]]:
    with path.open() as handle:
        reader = csv.DictReader(
            (line for line in handle if not line.startswith("#")), delimiter="\t"
        )
        output = {}
        for row in reader:
            if int(row["Site"]) == 349:
                state = row["State"]
                output[row["Node"]] = (state, float(row[f"p_{state}"]))
    return output


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(ALIGNMENT, "fasta"))
    record_ids = {record.id for record in records}
    if not set(AMPHIBIANS) <= record_ids:
        raise ValueError("One or more amphibian tips are missing from the ASR alignment")
    lineages = read_lineages()
    rows: list[dict[str, object]] = []

    for omitted in AMPHIBIANS:
        retained = [record for record in records if record.id != omitted]
        retained_ids = {record.id for record in retained}
        for topology, source_tree in TOPOLOGIES.items():
            run_dir = OUT / f"{topology}_minus_{omitted}"
            run_dir.mkdir(parents=True, exist_ok=True)
            alignment_path = run_dir / "alignment.fasta"
            tree_path = run_dir / "topology.nwk"
            SeqIO.write(retained, alignment_path, "fasta")

            tree = Phylo.read(source_tree, "newick")
            tree.prune(next(tip for tip in tree.get_terminals() if tip.name == omitted))
            Phylo.write(tree, tree_path, "newick")
            prefix = run_dir / "asr"
            subprocess.run(
                [
                    IQTREE2_BIN,
                    "-s",
                    str(alignment_path),
                    "-st",
                    "AA",
                    "-m",
                    "LG+F+R4",
                    "-te",
                    str(tree_path),
                    "-asr",
                    "-seed",
                    "20261102",
                    "-nt",
                    "4",
                    "-pre",
                    str(prefix),
                    "-redo",
                ],
                check=True,
            )
            states = read_np349_state(prefix.with_suffix(".state"))
            output_tree = Phylo.read(prefix.with_suffix(".treefile"), "newick")
            tip_names = {tip.name for tip in output_tree.get_terminals()}
            if tip_names != retained_ids:
                raise ValueError(f"Unexpected tips in {topology} minus {omitted}")

            for lineage in ("Amphibia", "Amniota", "Sauropsida", "Theria"):
                requested = lineages[lineage] & retained_ids
                query = copy.deepcopy(output_tree)
                query.root_with_outgroup(sorted(tip_names - requested)[0])
                mrca = query.common_ancestor(sorted(requested))
                descendants = {tip.name for tip in mrca.get_terminals()}
                state, posterior = states[mrca.name]
                rows.append(
                    {
                        "omitted_tip": omitted,
                        "topology": topology,
                        "lineage": lineage,
                        "mrca_node": mrca.name,
                        "requested_tips": len(requested),
                        "mrca_descendant_tips": len(descendants),
                        "monophyletic_exact": descendants == requested,
                        "np349_state": state,
                        "np349_posterior": posterior,
                    }
                )

    with RESULT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {RESULT}")


if __name__ == "__main__":
    main()
