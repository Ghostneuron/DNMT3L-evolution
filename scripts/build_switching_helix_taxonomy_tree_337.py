#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

from Bio import Phylo


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "analysis/switching_helix_asr_337"
TAXONOMY = ROOT / "taxonomy/ncbi_taxonomy_efetch.xml"
MAP = ROOT / "taxonomy/DNMT3L_289_taxid_map.tsv"
OUTPUT = ROOT / "trees/DNMT3L_289_ncbi_taxonomy_topology.nwk"
SUMMARY = PROJECT / "results/switching_helix_337_taxonomy_tree_summary.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    with MAP.open(newline="") as handle:
        map_rows = list(csv.DictReader(handle, delimiter="\t"))
    safe_by_taxid = {row["taxid"]: row["safe_id"] for row in map_rows}

    xml_root = ET.parse(TAXONOMY).getroot()
    paths: dict[str, list[str]] = {}
    name_by_taxid: dict[str, str] = {}
    for taxon in xml_root.findall("./Taxon"):
        tip_taxid = taxon.findtext("TaxId")
        if tip_taxid is None or tip_taxid not in safe_by_taxid:
            continue
        path: list[str] = []
        for ancestor in taxon.findall("./LineageEx/Taxon"):
            taxid = ancestor.findtext("TaxId")
            name = ancestor.findtext("ScientificName")
            if taxid and name:
                path.append(taxid)
                name_by_taxid[taxid] = name
        path.append(tip_taxid)
        name_by_taxid[tip_taxid] = taxon.findtext("ScientificName", default=tip_taxid)
        paths[tip_taxid] = path

    expected = set(safe_by_taxid)
    if set(paths) != expected:
        raise ValueError(f"Taxonomy mismatch: missing={sorted(expected-set(paths))}")

    ordered_paths = list(paths.values())
    common: list[str] = []
    for values in zip(*ordered_paths):
        if len(set(values)) != 1:
            break
        common.append(values[0])
    if not common:
        raise ValueError("No common NCBI taxonomy ancestor")
    root_taxid = common[-1]

    children: dict[str, set[str]] = defaultdict(set)
    for path in ordered_paths:
        trimmed = path[path.index(root_taxid) :]
        for parent, child in zip(trimmed, trimmed[1:]):
            children[parent].add(child)
    selected_internal = set(safe_by_taxid) & set(children)
    if selected_internal:
        raise ValueError(f"Selected tips are internal taxonomy nodes: {selected_internal}")

    collapsed_unary = 0

    def render(node: str, length: int | None) -> str:
        nonlocal collapsed_unary
        node_children = sorted(children.get(node, set()), key=lambda child: (name_by_taxid[child], child))
        if node in safe_by_taxid:
            body = safe_by_taxid[node]
        elif len(node_children) == 1 and node != root_taxid:
            collapsed_unary += 1
            return render(node_children[0], (length or 0) + 1)
        else:
            body = "(" + ",".join(render(child, 1) for child in node_children) + f")NCBI_{node}"
        return body if length is None else f"{body}:{length}"

    OUTPUT.write_text(render(root_taxid, None) + ";\n")
    tree = Phylo.read(OUTPUT, "newick")
    observed = {tip.name for tip in tree.get_terminals()}
    if observed != set(safe_by_taxid.values()):
        raise ValueError("Taxonomy tree tip set failed exact validation")
    polytomies = sum(len(clade.clades) > 2 for clade in tree.get_nonterminals())
    row = {
        "source": str(TAXONOMY.relative_to(PROJECT)),
        "source_sha256": sha256(TAXONOMY),
        "tip_count": len(observed),
        "tip_set_verified": True,
        "root_taxid": root_taxid,
        "root_name": name_by_taxid[root_taxid],
        "collapsed_unary_nodes": collapsed_unary,
        "polytomous_nodes": polytomies,
        "tree": str(OUTPUT.relative_to(PROJECT)),
        "tree_sha256": sha256(OUTPUT),
        "branch_length_rule": "one unit per NCBI taxonomy edge; IQ-TREE will re-estimate lengths",
    }
    with SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    print(row)


if __name__ == "__main__":
    main()
