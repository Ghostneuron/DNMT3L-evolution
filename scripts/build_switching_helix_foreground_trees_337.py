#!/usr/bin/env python3
from __future__ import annotations

import csv
import copy
from pathlib import Path

from Bio import Phylo, SeqIO


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "analysis/switching_helix_asr_337"
ALIGNMENT = ROOT / "inputs/DNMT3L_289_human_coordinate_strict_focal_masked_codon.fasta"
ATLAS = PROJECT / "results/functional_site_states_by_record.tsv"
TAX_TREE = ROOT / "trees/DNMT3L_289_ncbi_taxonomy_topology.nwk"
GENE_TREE = ROOT / "runs/GENE_TREE_STRICT_LG_F_R4/DNMT3L_289_gene_tree_strict.treefile"
OUT = ROOT / "trees"


def replace_once(source: Path, output: Path, node: str, label: str) -> None:
    text = source.read_text()
    old = f"{node}:"
    new = f"{node}{{{label}}}:"
    if text.count(old) != 1:
        raise ValueError(f"Expected one {old} in {source}, observed {text.count(old)}")
    output.write_text(text.replace(old, new))


def append_label(clade, label: str) -> None:
    clade.name = f"{clade.name or ''}{{{label}}}"


def lineage_sets() -> dict[str, set[str]]:
    by_tip: dict[str, str] = {}
    with ATLAS.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_tip.setdefault(row["safe_id"], row["broad_lineage"])
    output: dict[str, set[str]] = {}
    for tip, lineage in by_tip.items():
        output.setdefault(lineage, set()).add(tip)
    output["Theria"] = output["Marsupialia"] | output["Placentalia"]
    return output


def build_relax_tree(source: Path, output: Path, lineages: dict[str, set[str]], reroot: bool) -> None:
    tree = Phylo.read(source, "newick")
    if reroot:
        outgroup = sorted(lineages["Amphibia"])[0]
        tree.root_with_outgroup(outgroup)
    test = tree.common_ancestor(list(lineages["Theria"]))
    reference = tree.common_ancestor(list(lineages["Sauropsida"]))
    if {tip.name for tip in test.get_terminals()} != lineages["Theria"]:
        raise ValueError("Theria is not an exact clade in RELAX tree")
    if {tip.name for tip in reference.get_terminals()} != lineages["Sauropsida"]:
        raise ValueError("Sauropsida is not an exact clade in RELAX tree")
    for clade in test.find_clades(order="preorder"):
        append_label(clade, "Test")
    for clade in reference.find_clades(order="preorder"):
        append_label(clade, "Reference")
    Phylo.write(tree, output, "newick")


def validate_tips(path: Path, expected: set[str]) -> None:
    observed = {
        (tip.name or "").replace("{Foreground}", "").replace("{Test}", "").replace("{Reference}", "")
        for tip in Phylo.read(path, "newick").get_terminals()
    }
    if observed != expected:
        raise ValueError(f"Tip-set mismatch in {path}")


def main() -> None:
    expected = {record.id for record in SeqIO.parse(ALIGNMENT, "fasta")}
    lineages = lineage_sets()
    tax_stem = OUT / "DNMT3L_289_taxonomy_therian_stem_foreground.nwk"
    gene_stem = OUT / "DNMT3L_289_gene_tree_therian_stem_foreground.nwk"
    tax_relax = OUT / "DNMT3L_289_taxonomy_theria_test_sauropsida_reference.nwk"
    gene_relax = OUT / "DNMT3L_289_gene_tree_theria_test_sauropsida_reference.nwk"
    replace_once(TAX_TREE, tax_stem, "NCBI_32525", "Foreground")
    replace_once(GENE_TREE, gene_stem, "Node92", "Foreground")
    build_relax_tree(TAX_TREE, tax_relax, lineages, reroot=False)
    build_relax_tree(GENE_TREE, gene_relax, lineages, reroot=True)
    for path in (tax_stem, gene_stem, tax_relax, gene_relax):
        validate_tips(path, expected)
        print(path.relative_to(PROJECT))


if __name__ == "__main__":
    main()
