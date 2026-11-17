#!/usr/bin/env python3
from __future__ import annotations

import csv
import copy
import xml.etree.ElementTree as ET
from pathlib import Path

from Bio import Phylo, SeqIO


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "analysis/switching_helix_asr_337"
RESULTS = PROJECT / "results"
ALIGNMENT = ROOT / "inputs/DNMT3L_289_human_coordinate_strict_focal_masked.fasta"
ATLAS = RESULTS / "functional_site_states_by_record.tsv"
TAXONOMY_XML = ROOT / "taxonomy/ncbi_taxonomy_efetch.xml"
TAXONOMY_INPUT_TREE = ROOT / "trees/DNMT3L_289_ncbi_taxonomy_topology.nwk"
RUNS = {
    "taxonomy": ROOT / "runs/TAXONOMY_STRICT_LG_F_R4/DNMT3L_289_taxonomy_strict",
    "gene_tree": ROOT / "runs/GENE_TREE_STRICT_LG_F_R4/DNMT3L_289_gene_tree_strict",
}
SITES = (349, 352)
THRESHOLD = 0.80


def read_state(path: Path) -> dict[tuple[str, int], tuple[str, float]]:
    with path.open() as handle:
        lines = [line for line in handle if not line.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    output: dict[tuple[str, int], tuple[str, float]] = {}
    for row in reader:
        site = int(row["Site"])
        if site not in SITES:
            continue
        state = row["State"]
        output[(row["Node"], site)] = (state, float(row[f"p_{state}"]))
    return output


def observed_states() -> dict[tuple[str, int], tuple[str, float]]:
    output: dict[tuple[str, int], tuple[str, float]] = {}
    for record in SeqIO.parse(ALIGNMENT, "fasta"):
        sequence = str(record.seq).upper()
        for site in SITES:
            state = sequence[site - 1]
            if state not in {"-", "X", "?"}:
                output[(record.id, site)] = (state, 1.0)
    return output


def node_state(
    lookup: dict[tuple[str, int], tuple[str, float]], node: str | None, site: int
) -> tuple[str, float, str]:
    if node is None or (node, site) not in lookup:
        return "UNAVAILABLE", 0.0, "unavailable"
    state, posterior = lookup[(node, site)]
    call = state if posterior >= THRESHOLD else "UNCERTAIN"
    return call, posterior, state


def lineage_sets() -> dict[str, set[str]]:
    by_tip: dict[str, str] = {}
    with ATLAS.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_tip.setdefault(row["safe_id"], row["broad_lineage"])
    output: dict[str, set[str]] = {}
    for tip, lineage in by_tip.items():
        output.setdefault(lineage, set()).add(tip)
    return output


def taxonomy_names() -> dict[str, str]:
    output: dict[str, str] = {}
    root = ET.parse(TAXONOMY_XML).getroot()
    for taxon in root.findall("./Taxon"):
        taxid = taxon.findtext("TaxId")
        name = taxon.findtext("ScientificName")
        if taxid and name:
            output[taxid] = name
        for ancestor in taxon.findall("./LineageEx/Taxon"):
            taxid = ancestor.findtext("TaxId")
            name = ancestor.findtext("ScientificName")
            if taxid and name:
                output[taxid] = name
    return output


def main() -> None:
    observed = observed_states()
    expected_tips = {record.id for record in SeqIO.parse(ALIGNMENT, "fasta")}
    lineages = lineage_sets()
    lineages["Theria"] = lineages["Marsupialia"] | lineages["Placentalia"]
    lineages["Amniota"] = lineages["Theria"] | lineages["Sauropsida"]
    names = taxonomy_names()
    major_rows: list[dict[str, object]] = []
    transition_rows: list[dict[str, object]] = []

    for topology, prefix in RUNS.items():
        state_lookup = read_state(prefix.with_suffix(".state"))
        lookup = {**state_lookup, **observed}
        tree = Phylo.read(prefix.with_suffix(".treefile"), "newick")
        tip_names = {tip.name for tip in tree.get_terminals()}
        if tip_names != expected_tips:
            raise ValueError(f"{topology} output tree has incorrect tips")
        for lineage in ("Amphibia", "Amniota", "Sauropsida", "Theria", "Marsupialia", "Placentalia"):
            requested = lineages[lineage]
            # IQ-TREE writes reversible-model trees with an arbitrary display
            # root. Root each query on one tip outside the requested lineage so
            # a clade containing the display root is not falsely nonmonophyletic.
            query_tree = copy.deepcopy(tree)
            outgroup_tip = sorted(tip_names - requested)[0]
            query_tree.root_with_outgroup(outgroup_tip)
            mrca = query_tree.common_ancestor(list(requested))
            descendants = {tip.name for tip in mrca.get_terminals()}
            for site in SITES:
                call, posterior, raw_state = node_state(lookup, mrca.name, site)
                major_rows.append(
                    {
                        "topology": topology,
                        "lineage": lineage,
                        "mrca_node": mrca.name,
                        "requested_tips": len(requested),
                        "mrca_descendant_tips": len(descendants),
                        "monophyletic_exact": descendants == requested,
                        "site": site,
                        "called_state_pp80": call,
                        "maximum_posterior_state": raw_state,
                        "maximum_posterior": posterior,
                    }
                )

        if topology == "taxonomy":
            input_tree = Phylo.read(TAXONOMY_INPUT_TREE, "newick")

            def walk(parent) -> None:
                for child in parent.clades:
                    for site in SITES:
                        parent_call, parent_pp, parent_raw = node_state(lookup, parent.name, site)
                        child_call, child_pp, child_raw = node_state(lookup, child.name, site)
                        if (
                            parent_call not in {"UNAVAILABLE", "UNCERTAIN"}
                            and child_call not in {"UNAVAILABLE", "UNCERTAIN"}
                            and parent_call != child_call
                        ):
                            transition_rows.append(
                                {
                                    "topology": topology,
                                    "site": site,
                                    "parent_node": parent.name,
                                    "parent_taxon": names.get((parent.name or "").removeprefix("NCBI_"), ""),
                                    "parent_state": parent_call,
                                    "parent_posterior": parent_pp,
                                    "child_node": child.name,
                                    "child_taxon": names.get((child.name or "").removeprefix("NCBI_"), child.name or ""),
                                    "child_state": child_call,
                                    "child_posterior": child_pp,
                                    "child_tip_count": len(child.get_terminals()),
                                }
                            )
                    walk(child)

            walk(input_tree.root)

    key_nodes = {
        "NCBI_8292": "Amphibia",
        "NCBI_32524": "Amniota",
        "NCBI_32561": "Sauria",
        "NCBI_32525": "Theria",
        "NCBI_9263": "Marsupialia",
        "NCBI_9347": "Placentalia",
    }
    tax_lookup = {**read_state(RUNS["taxonomy"].with_suffix(".state")), **observed}
    for node, lineage in key_nodes.items():
        for site in SITES:
            call, posterior, raw_state = node_state(tax_lookup, node, site)
            major_rows.append(
                {
                    "topology": "taxonomy_named_node",
                    "lineage": lineage,
                    "mrca_node": node,
                    "requested_tips": "",
                    "mrca_descendant_tips": "",
                    "monophyletic_exact": True,
                    "site": site,
                    "called_state_pp80": call,
                    "maximum_posterior_state": raw_state,
                    "maximum_posterior": posterior,
                }
            )

    major_path = RESULTS / "switching_helix_337_asr_major_nodes.tsv"
    with major_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(major_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(major_rows)
    transition_path = RESULTS / "switching_helix_337_asr_high_confidence_transitions.tsv"
    with transition_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(transition_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(transition_rows)

    print(f"Major-node rows: {len(major_rows)}")
    print(f"High-confidence taxonomy transitions: {len(transition_rows)}")
    print("\nNamed taxonomy nodes")
    for row in major_rows:
        if row["topology"] == "taxonomy_named_node":
            print(row)
    print("\nCross-topology clade MRCAs")
    for row in major_rows:
        if row["topology"] in {"taxonomy", "gene_tree"}:
            print(row)


if __name__ == "__main__":
    main()
