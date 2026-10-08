from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "source_package" / "analysis" / "orthology_rebuild_20260826"
INPUTS = PROJECT / "analysis_v2" / "selection_inputs"
OUT = PROJECT / "analysis_v2" / "species_trees"
SCOPE_DATASETS = tuple(f"{scope}_MACSE" for scope in ("Placentalia", "Euarchontoglires", "Laurasiatheria", "Sauropsida"))


def read_names(path: Path) -> list[str]:
    return [line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")]


def clean_gene_id(value: object) -> str:
    text = str(value)
    return text[:-2] if text.endswith(".0") else text


def fetch_taxonomy(tax_ids: list[int]) -> dict[int, dict]:
    OUT.mkdir(parents=True, exist_ok=True)
    records: dict[int, dict] = {}
    for batch_index, start in enumerate(range(0, len(tax_ids), 80), start=1):
        batch = tax_ids[start : start + 80]
        joined = ",".join(str(value) for value in batch)
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{urllib.parse.quote(joined, safe=',')}"
        destination = OUT / f"ncbi_taxonomy_batch_{batch_index:03d}.json"
        if not destination.exists():
            request = urllib.request.Request(url, headers={"User-Agent": "DNMT3L-reanalysis/20261006"})
            with urllib.request.urlopen(request, timeout=120) as response:
                destination.write_bytes(response.read())
            time.sleep(0.4)
        payload = json.loads(destination.read_text())
        for node in payload.get("taxonomy_nodes", []):
            taxonomy = node.get("taxonomy", {})
            if "tax_id" in taxonomy:
                records[int(taxonomy["tax_id"])] = taxonomy
    missing = sorted(set(tax_ids) - set(records))
    if missing:
        raise ValueError(f"Missing NCBI taxonomy records: {missing}")
    return records


def make_tree(name_to_taxid: dict[str, int], taxonomy: dict[int, dict]) -> Tree:
    children: dict[int, set[int]] = defaultdict(set)
    tips_by_taxid: dict[int, list[str]] = defaultdict(list)
    nodes: set[int] = set()
    for name, tax_id in name_to_taxid.items():
        lineage = [int(value) for value in taxonomy[tax_id].get("lineage", [])]
        path = lineage + ([] if lineage and lineage[-1] == tax_id else [tax_id])
        if not path:
            raise ValueError(f"Empty lineage for {name} ({tax_id})")
        nodes.update(path)
        for parent, child in zip(path, path[1:]):
            children[parent].add(child)
        tips_by_taxid[tax_id].append(name)

    child_nodes = {child for values in children.values() for child in values}
    roots = sorted(nodes - child_nodes)
    if len(roots) != 1:
        raise ValueError(f"Expected one taxonomy root, found {roots}")

    def build(tax_id: int) -> Clade:
        clade = Clade(branch_length=1.0, name=f"tax_{tax_id}")
        for child in sorted(children.get(tax_id, set())):
            clade.clades.append(build(child))
        for tip_name in sorted(tips_by_taxid.get(tax_id, [])):
            clade.clades.append(Clade(branch_length=1.0, name=tip_name))
        return clade

    root = build(roots[0])

    def compress(clade: Clade, is_root: bool = False) -> Clade:
        clade.clades = [compress(child) for child in clade.clades]
        while len(clade.clades) == 1 and not clade.is_terminal():
            child = clade.clades[0]
            child.branch_length = (clade.branch_length or 0.0) + (child.branch_length or 0.0)
            clade = child
        return clade

    root = compress(root, is_root=True)
    root.branch_length = 0.0
    tree = Tree(root=root, rooted=True)
    observed = {tip.name for tip in tree.get_terminals()}
    expected = set(name_to_taxid)
    if observed != expected:
        raise ValueError(f"Species tree tip mismatch: missing={sorted(expected-observed)}, extra={sorted(observed-expected)}")
    return tree


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_csv(SOURCE / "orthology_audit" / "DNMT3L_orthology_evidence_matrix.csv", usecols=["safe_id", "gene_id", "Selection_Ready"])
    matrix = matrix.loc[matrix["Selection_Ready"].astype(bool)].copy()
    matrix["gene_id_text"] = matrix["gene_id"].map(clean_gene_id)

    report = SOURCE / "ncbi" / "DNMT3L_orthologs_ncbi" / "ncbi_dataset" / "data" / "data_report.jsonl"
    gene_to_taxid: dict[str, int] = {}
    for line in report.read_text().splitlines():
        record = json.loads(line)
        gene_to_taxid[str(record["geneId"])] = int(record["taxId"])
    matrix["tax_id"] = matrix["gene_id_text"].map(gene_to_taxid)
    missing_taxid = matrix.loc[matrix["tax_id"].isna(), ["safe_id", "gene_id_text"]]
    if not missing_taxid.empty:
        raise ValueError(f"No tax ID for records:\n{missing_taxid.to_string(index=False)}")
    safe_to_taxid = dict(zip(matrix["safe_id"].astype(str), matrix["tax_id"].astype(int)))

    required_names: set[str] = set()
    dataset_names: dict[str, list[str]] = {}
    for dataset in SCOPE_DATASETS:
        path = INPUTS / f"{dataset}_selection_input.fasta"
        names = read_names(path)
        dataset_names[dataset] = names
        required_names.update(names)
    missing_names = sorted(required_names - set(safe_to_taxid))
    if missing_names:
        raise ValueError(f"Selection names absent from orthology matrix: {missing_names}")

    tax_ids = sorted({safe_to_taxid[name] for name in required_names})
    taxonomy = fetch_taxonomy(tax_ids)
    mapping_rows: list[dict[str, object]] = []
    for dataset, names in dataset_names.items():
        mapping = {name: safe_to_taxid[name] for name in names}
        tree = make_tree(mapping, taxonomy)
        output = OUT / f"{dataset}_NCBI_taxonomy_tree.nwk"
        Phylo.write(tree, output, "newick")
        reread = Phylo.read(output, "newick")
        if {tip.name for tip in reread.get_terminals()} != set(names):
            raise ValueError(f"Round-trip tree mismatch for {dataset}")
        for name in names:
            mapping_rows.append(
                {
                    "Dataset": dataset,
                    "Sequence_Name": name,
                    "Tax_ID": mapping[name],
                    "Organism_Name": taxonomy[mapping[name]].get("organism_name", ""),
                    "NCBI_Lineage_Tax_IDs": ";".join(str(value) for value in taxonomy[mapping[name]].get("lineage", [])),
                }
            )
        print(f"{dataset}: {len(names)} tips -> {output.name}")
    pd.DataFrame(mapping_rows).to_csv(OUT / "species_tree_taxonomy_manifest.csv", index=False)


if __name__ == "__main__":
    main()
