from __future__ import annotations

import argparse
import copy
import json
import math
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp


THERIAN_CLADES = {
    "Afrotheria",
    "Euarchontoglires",
    "Laurasiatheria",
    "Marsupialia",
    "Xenarthra",
}


def clean_gene_id(value: object) -> str:
    text = str(value)
    return text[:-2] if text.endswith(".0") else text


def load_gene_taxids(report_path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    for line in report_path.read_text().splitlines():
        record = json.loads(line)
        result[str(record["geneId"])] = int(record["taxId"])
    return result


def load_cached_taxonomy(cache_paths: list[Path]) -> dict[int, dict]:
    records: dict[int, dict] = {}
    for path in cache_paths:
        payload = json.loads(path.read_text())
        for node in payload.get("taxonomy_nodes", []):
            taxonomy = node.get("taxonomy", {})
            if "tax_id" in taxonomy:
                records[int(taxonomy["tax_id"])] = taxonomy
    return records


def fetch_missing_taxonomy(tax_ids: list[int], outdir: Path) -> dict[int, dict]:
    records: dict[int, dict] = {}
    for batch_index, start in enumerate(range(0, len(tax_ids), 80), start=1):
        batch = tax_ids[start : start + 80]
        joined = ",".join(str(value) for value in batch)
        url = (
            "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/"
            + urllib.parse.quote(joined, safe=",")
        )
        destination = outdir / f"ncbi_taxonomy_added_{batch_index:03d}.json"
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "DNMT3L-revision/20261025"},
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            payload_bytes = response.read()
        destination.write_bytes(payload_bytes)
        payload = json.loads(payload_bytes)
        for node in payload.get("taxonomy_nodes", []):
            taxonomy = node.get("taxonomy", {})
            if "tax_id" in taxonomy:
                records[int(taxonomy["tax_id"])] = taxonomy
        time.sleep(0.4)
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

    def compress(clade: Clade) -> Clade:
        clade.clades = [compress(child) for child in clade.clades]
        while len(clade.clades) == 1 and not clade.is_terminal():
            child = clade.clades[0]
            child.branch_length = (clade.branch_length or 0.0) + (child.branch_length or 0.0)
            clade = child
        return clade

    root = compress(build(roots[0]))
    root.branch_length = 0.0
    tree = Tree(root=root, rooted=True)
    observed = {tip.name for tip in tree.get_terminals()}
    expected = set(name_to_taxid)
    if observed != expected:
        raise ValueError(
            f"Species tree tip mismatch: missing={sorted(expected-observed)}, extra={sorted(observed-expected)}"
        )
    return tree


def prune_to_names(tree: Tree, keep_names: set[str]) -> Tree:
    result = copy.deepcopy(tree)
    for tip in list(result.get_terminals()):
        if tip.name not in keep_names:
            result.prune(tip)
    result.root.branch_length = 0.0
    return result


def sankoff_costs(tree: Tree, tip_states: dict[str, str], states: list[str]) -> tuple[dict[int, np.ndarray], int]:
    index = {state: i for i, state in enumerate(states)}
    costs: dict[int, np.ndarray] = {}
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            values = np.full(len(states), np.inf)
            values[index[tip_states[node.name]]] = 0.0
        else:
            values = np.zeros(len(states))
            for child in node.clades:
                child_cost = costs[id(child)]
                transition = np.ones((len(states), len(states)))
                np.fill_diagonal(transition, 0.0)
                values += np.min(transition + child_cost[None, :], axis=1)
        costs[id(node)] = values
    minimum = int(np.min(costs[id(tree.root)]))
    return costs, minimum


def er_log_transition(k: int, rate: float, branch_length: float) -> np.ndarray:
    branch_length = max(float(branch_length), 1e-9)
    exponent = math.exp(-(k * rate / (k - 1)) * branch_length) if k > 1 else 1.0
    off = (1.0 - exponent) / k
    matrix = np.full((k, k), off)
    np.fill_diagonal(matrix, (1.0 / k) + ((k - 1.0) / k) * exponent)
    return np.log(np.maximum(matrix, 1e-300))


def er_messages(
    tree: Tree,
    tip_states: dict[str, str],
    states: list[str],
    rate: float,
) -> tuple[float, dict[int, np.ndarray]]:
    index = {state: i for i, state in enumerate(states)}
    upward: dict[int, np.ndarray] = {}
    child_message: dict[tuple[int, int], np.ndarray] = {}
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            values = np.full(len(states), -np.inf)
            values[index[tip_states[node.name]]] = 0.0
        else:
            values = np.zeros(len(states))
            for child in node.clades:
                log_p = er_log_transition(len(states), rate, child.branch_length or 1.0)
                message = logsumexp(log_p + upward[id(child)][None, :], axis=1)
                child_message[(id(node), id(child))] = message
                values += message
        upward[id(node)] = values

    log_prior = np.full(len(states), -math.log(len(states)))
    log_likelihood = float(logsumexp(log_prior + upward[id(tree.root)]))
    downward: dict[int, np.ndarray] = {id(tree.root): log_prior}
    for node in tree.find_clades(order="preorder"):
        if node.is_terminal():
            continue
        for child in node.clades:
            outside_parent = downward[id(node)].copy()
            for sibling in node.clades:
                if sibling is not child:
                    outside_parent += child_message[(id(node), id(sibling))]
            log_p = er_log_transition(len(states), rate, child.branch_length or 1.0)
            downward[id(child)] = logsumexp(outside_parent[:, None] + log_p, axis=0)

    posterior: dict[int, np.ndarray] = {}
    for node in tree.find_clades():
        values = downward[id(node)] + upward[id(node)]
        posterior[id(node)] = np.exp(values - logsumexp(values))
    return log_likelihood, posterior


def fit_er(tree: Tree, tip_states: dict[str, str], states: list[str]) -> tuple[float, float, dict[int, np.ndarray]]:
    if len(states) == 1:
        posterior = {id(node): np.array([1.0]) for node in tree.find_clades()}
        return 0.0, 0.0, posterior

    def objective(log_rate: float) -> float:
        rate = math.exp(log_rate)
        log_likelihood, _ = er_messages(tree, tip_states, states, rate)
        return -log_likelihood

    fitted = minimize_scalar(objective, bounds=(-12.0, 5.0), method="bounded")
    rate = math.exp(float(fitted.x))
    log_likelihood, posterior = er_messages(tree, tip_states, states, rate)
    return rate, log_likelihood, posterior


def top_posterior(states: list[str], probabilities: np.ndarray, n: int = 5) -> str:
    order = np.argsort(probabilities)[::-1][:n]
    return ";".join(f"{states[i]}:{probabilities[i]:.4f}" for i in order)


def analyze_position(
    full_tree: Tree,
    metadata: pd.DataFrame,
    residue_table: pd.DataFrame,
    position: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows = residue_table.loc[
        (residue_table["Human_Position"] == position)
        & residue_table["Query_Residue"].notna()
        & (residue_table["Query_Residue"].astype(str).str.len() == 1)
    ].copy()
    allowed = set(metadata["safe_id"])
    rows = rows.loc[rows["safe_id"].isin(allowed)]
    tip_states = dict(zip(rows["safe_id"].astype(str), rows["Query_Residue"].astype(str)))
    tree = prune_to_names(full_tree, set(tip_states))
    states = sorted(set(tip_states.values()))
    meta_by_id = metadata.set_index("safe_id")
    therian_names = [name for name in tip_states if meta_by_id.loc[name, "Clade"] in THERIAN_CLADES]
    sauropsid_names = [name for name in tip_states if meta_by_id.loc[name, "Clade"] == "Sauropsida"]
    therian_mrca = tree.common_ancestor(therian_names)
    sauropsid_mrca = tree.common_ancestor(sauropsid_names)
    sankoff, minimum_changes = sankoff_costs(tree, tip_states, states)

    summary_rows: list[dict[str, object]] = []
    for scheme in ("taxonomy_step_lengths", "unit_branch_lengths"):
        model_tree = copy.deepcopy(tree)
        if scheme == "unit_branch_lengths":
            for node in model_tree.find_clades():
                if node is not model_tree.root:
                    node.branch_length = 1.0
        model_therian_mrca = model_tree.common_ancestor(therian_names)
        model_sauropsid_mrca = model_tree.common_ancestor(sauropsid_names)
        rate, log_likelihood, posterior = fit_er(model_tree, tip_states, states)
        therian_prob = posterior[id(model_therian_mrca)]
        sauropsid_prob = posterior[id(model_sauropsid_mrca)]
        state_index = {state: i for i, state in enumerate(states)}
        summary_rows.append(
            {
                "Human_Position": position,
                "Branch_Length_Scheme": scheme,
                "Resolved_Tips": len(tip_states),
                "Therian_Tips": len(therian_names),
                "Sauropsid_Tips": len(sauropsid_names),
                "Observed_States": ";".join(states),
                "Parsimony_Minimum_Changes": minimum_changes,
                "Therian_MRCA_Descendant_Parsimony_States": ";".join(
                    state for state, value in zip(states, sankoff[id(therian_mrca)]) if value == np.min(sankoff[id(therian_mrca)])
                ),
                "Sauropsid_MRCA_Descendant_Parsimony_States": ";".join(
                    state for state, value in zip(states, sankoff[id(sauropsid_mrca)]) if value == np.min(sankoff[id(sauropsid_mrca)])
                ),
                "ER_Rate": rate,
                "ER_LogLikelihood": log_likelihood,
                "Therian_MRCA_Posterior": top_posterior(states, therian_prob),
                "Sauropsid_MRCA_Posterior": top_posterior(states, sauropsid_prob),
                "Therian_MRCA_Q_Probability": therian_prob[state_index["Q"]] if "Q" in state_index else 0.0,
                "Therian_MRCA_K_Probability": therian_prob[state_index["K"]] if "K" in state_index else 0.0,
                "Sauropsid_MRCA_Q_Probability": sauropsid_prob[state_index["Q"]] if "Q" in state_index else 0.0,
                "Sauropsid_MRCA_K_Probability": sauropsid_prob[state_index["K"]] if "K" in state_index else 0.0,
            }
        )

    count_rows: list[dict[str, object]] = []
    for group, names in (("Theria", therian_names), ("Sauropsida", sauropsid_names)):
        counts = Counter(tip_states[name] for name in names)
        for residue, count in sorted(counts.items()):
            count_rows.append(
                {
                    "Human_Position": position,
                    "Group": group,
                    "Residue": residue,
                    "Count": count,
                    "Resolved_N": len(names),
                    "Frequency": count / len(names),
                }
            )
    return summary_rows, count_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    source = project / "source_package" / "analysis" / "orthology_rebuild_20260826"
    matrix_path = source / "orthology_audit" / "DNMT3L_orthology_evidence_matrix.csv"
    report_path = source / "ncbi" / "DNMT3L_orthologs_ncbi" / "ncbi_dataset" / "data" / "data_report.jsonl"
    residue_path = project / "analysis_v2" / "summaries" / "structure_informed_residue_states_by_record.csv"
    cache_paths = sorted((project / "analysis_v2" / "species_trees").glob("ncbi_taxonomy_batch_*.json"))

    matrix = pd.read_csv(matrix_path)
    metadata = matrix.loc[
        matrix["Selection_Ready"].astype(bool)
        & ~matrix["Hybrid_Taxon"].astype(bool)
        & matrix["Clade"].isin(THERIAN_CLADES | {"Sauropsida"}),
        ["safe_id", "species", "gene_id", "Clade"],
    ].copy()
    gene_to_taxid = load_gene_taxids(report_path)
    metadata["gene_id_text"] = metadata["gene_id"].map(clean_gene_id)
    metadata["tax_id"] = metadata["gene_id_text"].map(gene_to_taxid)
    if metadata["tax_id"].isna().any():
        raise ValueError("Missing tax IDs for one or more selected records")
    metadata["tax_id"] = metadata["tax_id"].astype(int)

    taxonomy = load_cached_taxonomy(cache_paths)
    missing_taxids = sorted(set(metadata["tax_id"]) - set(taxonomy))
    if missing_taxids:
        taxonomy.update(fetch_missing_taxonomy(missing_taxids, outdir))
    still_missing = sorted(set(metadata["tax_id"]) - set(taxonomy))
    if still_missing:
        raise ValueError(f"Missing NCBI taxonomy records after fetch: {still_missing}")

    name_to_taxid = dict(zip(metadata["safe_id"].astype(str), metadata["tax_id"].astype(int)))
    tree = make_tree(name_to_taxid, taxonomy)
    tree_path = outdir / "theria_sauropsida_NCBI_taxonomy_tree.nwk"
    Phylo.write(tree, tree_path, "newick")

    taxonomy_rows = []
    for row in metadata.itertuples(index=False):
        record = taxonomy[int(row.tax_id)]
        taxonomy_rows.append(
            {
                "safe_id": row.safe_id,
                "species": row.species,
                "Clade": row.Clade,
                "Gene_ID": row.gene_id,
                "Tax_ID": int(row.tax_id),
                "NCBI_Organism_Name": record.get("organism_name", ""),
                "NCBI_Lineage_Tax_IDs": ";".join(str(value) for value in record.get("lineage", [])),
            }
        )
    pd.DataFrame(taxonomy_rows).to_csv(outdir / "ancestral_state_taxonomy_manifest.csv", index=False)

    residue_table = pd.read_csv(residue_path)
    all_summary: list[dict[str, object]] = []
    all_counts: list[dict[str, object]] = []
    for position in (349, 352):
        summary_rows, count_rows = analyze_position(tree, metadata, residue_table, position)
        all_summary.extend(summary_rows)
        all_counts.extend(count_rows)
    pd.DataFrame(all_summary).to_csv(outdir / "ancestral_state_summary.csv", index=False)
    pd.DataFrame(all_counts).to_csv(outdir / "ancestral_state_tip_counts.csv", index=False)
    print(pd.DataFrame(all_summary).to_string(index=False))


if __name__ == "__main__":
    main()
