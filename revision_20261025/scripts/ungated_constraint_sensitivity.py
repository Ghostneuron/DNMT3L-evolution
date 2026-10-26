from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ancestral_state_analysis import (  # noqa: E402
    clean_gene_id,
    fetch_missing_taxonomy,
    load_cached_taxonomy,
    load_gene_taxids,
    make_tree,
    prune_to_names,
    sankoff_costs,
)


PLACENTAL_CLADES = {"Afrotheria", "Euarchontoglires", "Laurasiatheria", "Xenarthra"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_fasta(frame: pd.DataFrame, path: Path) -> None:
    with path.open("w") as handle:
        for row in frame.itertuples(index=False):
            handle.write(f">{row.safe_id}\n")
            sequence = str(row.protein_seq).upper()
            for start in range(0, len(sequence), 90):
                handle.write(sequence[start : start + 90] + "\n")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def human_position_map(aligned_reference: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    position = 0
    for column, residue in enumerate(aligned_reference):
        if residue not in {"-", "."}:
            position += 1
            mapping[position] = column
    return mapping


def domain_for(position: int) -> str:
    if position <= 40:
        return "N_terminal"
    if position <= 173:
        return "ADD"
    if position <= 177:
        return "linker"
    if position <= 380:
        return "C_terminal_MTase_like"
    return "C_terminal_tail"


def site_table(records: list[tuple[str, str]], label: str) -> pd.DataFrame:
    by_name = dict(records)
    reference = by_name["Homo_sapiens"]
    mapping = human_position_map(reference)
    rows = []
    for human_position, column in mapping.items():
        counts = Counter(
            sequence[column]
            for name, sequence in records
            if sequence[column] not in {"-", ".", "X", "?", "*"}
        )
        resolved = sum(counts.values())
        occupancy = resolved / len(records)
        if occupancy < 0.80:
            continue
        entropy = 0.0
        for count in counts.values():
            frequency = count / resolved
            entropy -= frequency * math.log2(frequency)
        rows.append(
            {
                "Cohort": label,
                "Human_Position": human_position,
                "Domain": domain_for(human_position),
                "Resolved_N": resolved,
                "Occupancy": occupancy,
                "Entropy": entropy,
                "Major_Residue": counts.most_common(1)[0][0],
                "Major_Residue_Frequency": counts.most_common(1)[0][1] / resolved,
                "Amino_Acid_Counts": ";".join(f"{aa}:{count}" for aa, count in sorted(counts.items())),
            }
        )
    return pd.DataFrame(rows)


def parsimony_table(
    records: list[tuple[str, str]],
    sites: pd.DataFrame,
    full_tree,
    label: str,
) -> pd.DataFrame:
    by_name = dict(records)
    mapping = human_position_map(by_name["Homo_sapiens"])
    rows = []
    for site in sites.itertuples(index=False):
        column = mapping[int(site.Human_Position)]
        tip_states = {
            name: sequence[column]
            for name, sequence in records
            if sequence[column] not in {"-", ".", "X", "?", "*"}
        }
        tree = prune_to_names(full_tree, set(tip_states))
        states = sorted(set(tip_states.values()))
        _, minimum_changes = sankoff_costs(tree, tip_states, states)
        rows.append(
            {
                "Cohort": label,
                "Human_Position": int(site.Human_Position),
                "Domain": site.Domain,
                "Resolved_N": len(tip_states),
                "Minimum_Changes": minimum_changes,
                "Minimum_Changes_Per_100_Tips": 100.0 * minimum_changes / len(tip_states),
            }
        )
    return pd.DataFrame(rows)


def summarize_sites(sites: pd.DataFrame, parsimony: pd.DataFrame, label: str) -> dict[str, object]:
    add = sites.loc[sites["Domain"] == "ADD", "Entropy"]
    ctd = sites.loc[sites["Domain"] == "C_terminal_MTase_like", "Entropy"]
    add_tree = parsimony.loc[parsimony["Domain"] == "ADD", "Minimum_Changes_Per_100_Tips"]
    ctd_tree = parsimony.loc[
        parsimony["Domain"] == "C_terminal_MTase_like", "Minimum_Changes_Per_100_Tips"
    ]
    entropy_test = mannwhitneyu(add, ctd, alternative="two-sided")
    tree_test = mannwhitneyu(add_tree, ctd_tree, alternative="two-sided")
    result: dict[str, object] = {
        "Cohort": label,
        "Sequences": int(sites["Resolved_N"].max()),
        "Stable_Sites": len(sites),
        "ADD_Sites": len(add),
        "CTD_Sites": len(ctd),
        "ADD_Median_Entropy": float(add.median()),
        "CTD_Median_Entropy": float(ctd.median()),
        "ADD_vs_CTD_Entropy_P": float(entropy_test.pvalue),
        "ADD_Median_Minimum_Changes_Per_100_Tips": float(add_tree.median()),
        "CTD_Median_Minimum_Changes_Per_100_Tips": float(ctd_tree.median()),
        "ADD_vs_CTD_Tree_P": float(tree_test.pvalue),
    }
    for position in (349, 352):
        row = sites.loc[sites["Human_Position"] == position]
        result[f"NP{position}_Counts"] = row.iloc[0]["Amino_Acid_Counts"] if not row.empty else "unresolved"
        result[f"NP{position}_Resolved_N"] = int(row.iloc[0]["Resolved_N"]) if not row.empty else 0
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--mafft", default="mafft")
    args = parser.parse_args()
    project = args.project_root.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    source = project / "source_package" / "analysis" / "orthology_rebuild_20260826"
    matrix_path = source / "orthology_audit" / "DNMT3L_orthology_evidence_matrix.csv"
    report_path = source / "ncbi" / "DNMT3L_orthologs_ncbi" / "ncbi_dataset" / "data" / "data_report.jsonl"
    matrix = pd.read_csv(matrix_path)
    placental = matrix.loc[matrix["Clade"].isin(PLACENTAL_CLADES)].copy()

    gated = placental.loc[placental["Selection_Ready"].astype(bool)].copy()
    ungated_mask = (
        placental["NCBI_DNMT3L_Ortholog_Group"].astype(bool)
        & ~placental["NCBI_Paralog_Group_Conflict"].astype(bool)
        & placental["Sequence_QC_Pass"].astype(bool)
        & (placental["Full_Reference_Coverage"] >= 0.80)
        & (placental["CTD_Reference_Coverage"] >= 0.90)
        & ~placental["Hybrid_Taxon"].astype(bool)
    )
    ungated = placental.loc[ungated_mask].copy()
    added = ungated.loc[~ungated["safe_id"].isin(gated["safe_id"])].copy()
    if len(added) != 5:
        raise ValueError(f"Expected five h358-gate additions in Placentalia, found {len(added)}")
    added[
        [
            "safe_id",
            "species",
            "gene_id",
            "protein_accession",
            "Full_Reference_Coverage",
            "CTD_Reference_Coverage",
        ]
    ].to_csv(outdir / "ungated_added_placental_records.csv", index=False)

    gene_to_taxid = load_gene_taxids(report_path)
    all_metadata = ungated[["safe_id", "species", "gene_id", "Clade"]].copy()
    all_metadata["gene_id_text"] = all_metadata["gene_id"].map(clean_gene_id)
    all_metadata["tax_id"] = all_metadata["gene_id_text"].map(gene_to_taxid)
    if all_metadata["tax_id"].isna().any():
        raise ValueError("Missing tax ID in ungated cohort")
    all_metadata["tax_id"] = all_metadata["tax_id"].astype(int)
    cache_paths = sorted((project / "analysis_v2" / "species_trees").glob("ncbi_taxonomy_batch_*.json"))
    cache_paths += sorted(outdir.parent.glob("ancestral_state/ncbi_taxonomy_added_*.json"))
    taxonomy = load_cached_taxonomy(cache_paths)
    missing_taxids = sorted(set(all_metadata["tax_id"]) - set(taxonomy))
    if missing_taxids:
        taxonomy.update(fetch_missing_taxonomy(missing_taxids, outdir))
    if set(all_metadata["tax_id"]) - set(taxonomy):
        raise ValueError("Missing taxonomy after fetch")

    summary_rows = []
    site_frames = []
    parsimony_frames = []
    command_rows = []
    for label, frame in (("h358_gated", gated), ("h358_ungated", ungated)):
        frame = frame.sort_values("safe_id").copy()
        input_fasta = outdir / f"Placentalia_{label}_protein.fasta"
        alignment_fasta = outdir / f"Placentalia_{label}_protein_MAFFT.fasta"
        log_path = outdir / f"Placentalia_{label}_MAFFT.log"
        write_fasta(frame, input_fasta)
        command = [args.mafft, "--auto", "--thread", "8", str(input_fasta)]
        with alignment_fasta.open("w") as output, log_path.open("w") as log:
            completed = subprocess.run(command, stdout=output, stderr=log, text=True, check=False)
        if completed.returncode:
            raise subprocess.CalledProcessError(completed.returncode, command)
        records = read_fasta(alignment_fasta)
        if len(records) != len(frame):
            raise ValueError(f"Alignment record mismatch for {label}")
        sites = site_table(records, label)

        metadata = all_metadata.loc[all_metadata["safe_id"].isin(frame["safe_id"])].copy()
        name_to_taxid = dict(zip(metadata["safe_id"].astype(str), metadata["tax_id"].astype(int)))
        tree = make_tree(name_to_taxid, taxonomy)
        tree_path = outdir / f"Placentalia_{label}_NCBI_taxonomy_tree.nwk"
        Phylo.write(tree, tree_path, "newick")
        parsimony = parsimony_table(records, sites, tree, label)

        site_frames.append(sites)
        parsimony_frames.append(parsimony)
        summary_rows.append(summarize_sites(sites, parsimony, label))
        command_rows.append(
            {
                "Cohort": label,
                "Input_FASTA": str(input_fasta),
                "Input_SHA256": sha256(input_fasta),
                "Alignment_FASTA": str(alignment_fasta),
                "Alignment_SHA256": sha256(alignment_fasta),
                "Command": " ".join(command),
                "Sequence_Count": len(frame),
            }
        )

    summary = pd.DataFrame(summary_rows)
    gated_row = summary.loc[summary["Cohort"] == "h358_gated"].iloc[0]
    ungated_row = summary.loc[summary["Cohort"] == "h358_ungated"].iloc[0]
    summary["Delta_ADD_Median_Entropy_vs_Gated"] = summary["ADD_Median_Entropy"] - gated_row["ADD_Median_Entropy"]
    summary["Delta_CTD_Median_Entropy_vs_Gated"] = summary["CTD_Median_Entropy"] - gated_row["CTD_Median_Entropy"]
    summary["Delta_ADD_Median_Minimum_Changes_vs_Gated"] = (
        summary["ADD_Median_Minimum_Changes_Per_100_Tips"]
        - gated_row["ADD_Median_Minimum_Changes_Per_100_Tips"]
    )
    summary["Delta_CTD_Median_Minimum_Changes_vs_Gated"] = (
        summary["CTD_Median_Minimum_Changes_Per_100_Tips"]
        - gated_row["CTD_Median_Minimum_Changes_Per_100_Tips"]
    )
    summary.to_csv(outdir / "ungated_constraint_summary.csv", index=False)
    pd.concat(site_frames, ignore_index=True).to_csv(outdir / "ungated_constraint_by_site.csv", index=False)
    pd.concat(parsimony_frames, ignore_index=True).to_csv(outdir / "ungated_tree_parsimony_by_site.csv", index=False)
    pd.DataFrame(command_rows).to_csv(outdir / "ungated_alignment_manifest.csv", index=False)
    print("Added records:")
    print(added[["safe_id", "species", "Clade"]].to_string(index=False))
    print("\nSummary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
