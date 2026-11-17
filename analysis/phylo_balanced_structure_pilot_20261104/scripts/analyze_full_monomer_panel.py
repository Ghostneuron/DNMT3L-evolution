#!/usr/bin/env python3
"""Confidence-gated analysis of 24 DNMT3L monomers across three prediction seeds."""

from __future__ import annotations

import csv
import itertools
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.SVDSuperimposer import SVDSuperimposer


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/full_monomer_analysis"
DOMAINS = {"ADD": (49, 181), "CTD": (178, 380)}
ANCESTORS = {"ANCESTRAL_AMNIOTA", "ANCESTRAL_THERIA"}


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    if not rows:
        raise ValueError(f"Empty table: {path}")
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def one(path, pattern):
    found = [item for item in path.glob(pattern) if not item.name.startswith("._")]
    if len(found) != 1:
        raise ValueError((path, pattern, len(found)))
    return found[0]


def rmsd(reference, moving):
    if len(reference) < 3 or len(moving) < 3:
        return float("nan")
    sup = SVDSuperimposer()
    sup.set(np.asarray(reference, float), np.asarray(moving, float))
    sup.run()
    return float(sup.get_rms())


def plot_heatmap(matrix, labels, groups, title, path):
    masked = np.ma.masked_invalid(matrix)
    fig, ax = plt.subplots(figsize=(10.5, 9.2))
    image = ax.imshow(masked, cmap="viridis", vmin=0, vmax=float(np.nanpercentile(matrix, 95)))
    ax.set_xticks(range(len(labels)), labels, rotation=90, fontsize=6)
    ax.set_yticks(range(len(labels)), labels, fontsize=6)
    ax.set_title(title)
    previous = groups[0]
    for index, group in enumerate(groups[1:], start=1):
        if group != previous:
            ax.axhline(index - 0.5, color="white", linewidth=0.8)
            ax.axvline(index - 0.5, color="white", linewidth=0.8)
        previous = group
    fig.colorbar(image, ax=ax, label="Matched-seed mean C-alpha RMSD (A)", shrink=0.75)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=240)
    plt.close(fig)


def main():
    status = read_tsv(ROOT / "metadata/full_monomer_run_status.tsv")
    if len(status) != 72 or any(row["status"] not in {"VALIDATED", "VALIDATED_REUSED"} for row in status):
        raise ValueError("Full monomer completion gate has not passed")
    panel_rows = read_tsv(ROOT / "metadata/panel_manifest.tsv")
    panel = {row["safe_id"]: row for row in panel_rows}
    targets = [row["safe_id"] for row in panel_rows]
    seeds = [4101, 4102, 4103]
    crosswalk = defaultdict(dict)
    for row in read_tsv(ROOT / "metadata/human_coordinate_crosswalk.tsv"):
        crosswalk[row["safe_id"]][int(row["human_position"])] = {
            "index": int(row["monomer_input_position"]) - 1,
            "pp": float(row["ancestral_posterior_probability"]) if row["ancestral_posterior_probability"] else None,
        }

    predictions = {}
    prediction_rows = []
    for row in status:
        target, seed = row["target"], int(row["seed"])
        output = ROOT / row["output"]
        plddt = np.load(one(output, "plddt_*.npz"))["plddt"]
        residues = list(MMCIFParser(QUIET=True).get_structure("prediction", one(output, "*.cif"))[0]["A"].get_residues())
        predictions[(target, seed)] = {"plddt": plddt, "residues": residues}
        result = {"target": target, "lineage_group": panel[target]["lineage_group"], "seed": seed}
        for domain, (start, end) in DOMAINS.items():
            indices = [crosswalk[target][human]["index"] for human in range(start, end + 1) if human in crosswalk[target]]
            fraction = float((plddt[indices] >= .70).mean())
            result[f"{domain}_mapped_n"] = len(indices)
            result[f"{domain}_plddt_ge_0.70_fraction"] = f"{fraction:.6f}"
            result[f"{domain}_confidence_pass"] = str(fraction >= .80)
        prediction_rows.append(result)
    write_tsv(OUT / "per_prediction_domain_confidence.tsv", prediction_rows)

    within_rows = []
    within_max = defaultdict(dict)
    target_gate_rows = []
    for target in targets:
        gate = {"target": target, "lineage_group": panel[target]["lineage_group"]}
        for domain, (start, end) in DOMAINS.items():
            pass_n = sum(row[f"{domain}_confidence_pass"] == "True" for row in prediction_rows if row["target"] == target)
            values = []
            for seed_a, seed_b in itertools.combinations(seeds, 2):
                left, right = predictions[(target, seed_a)], predictions[(target, seed_b)]
                positions = []
                for human in range(start, end + 1):
                    mapping = crosswalk[target].get(human)
                    if not mapping:
                        continue
                    index = mapping["index"]
                    if left["plddt"][index] < .70 or right["plddt"][index] < .70:
                        continue
                    if target in ANCESTORS and mapping["pp"] < .80:
                        continue
                    positions.append(human)
                value = rmsd(
                    [left["residues"][crosswalk[target][h]["index"]]["CA"].coord for h in positions],
                    [right["residues"][crosswalk[target][h]["index"]]["CA"].coord for h in positions],
                )
                values.append(value)
                within_rows.append({
                    "target": target, "lineage_group": panel[target]["lineage_group"], "domain": domain,
                    "seed_a": seed_a, "seed_b": seed_b, "mapped_high_confidence_n": len(positions), "rmsd_A": f"{value:.6f}",
                })
            finite_values = [value for value in values if np.isfinite(value)]
            within_value = max(finite_values) if finite_values else float("inf")
            within_max[target][domain] = within_value
            gate[f"{domain}_confidence_pass_n"] = pass_n
            gate[f"{domain}_within_seed_rmsd_max_A"] = f"{within_value:.6f}"
            gate[f"{domain}_target_gate"] = str(pass_n >= 2 and within_value <= 2.0)
        target_gate_rows.append(gate)
    write_tsv(OUT / "within_target_seed_dispersion.tsv", within_rows)
    write_tsv(OUT / "per_target_domain_gate.tsv", target_gate_rows)
    target_gate = {(row["target"], domain): row[f"{domain}_target_gate"] == "True" for row in target_gate_rows for domain in DOMAINS}

    contrast_rows = []
    matrices = {domain: np.full((len(targets), len(targets)), np.nan) for domain in DOMAINS}
    for domain in DOMAINS:
        for index, target in enumerate(targets):
            if target_gate[(target, domain)]:
                matrices[domain][index, index] = 0.0
    for left_index, right_index in itertools.combinations(range(len(targets)), 2):
        target_a, target_b = targets[left_index], targets[right_index]
        for domain, (start, end) in DOMAINS.items():
            values, counts = [], []
            for seed in seeds:
                left, right = predictions[(target_a, seed)], predictions[(target_b, seed)]
                positions = []
                for human in range(start, end + 1):
                    map_a, map_b = crosswalk[target_a].get(human), crosswalk[target_b].get(human)
                    if not map_a or not map_b:
                        continue
                    if left["plddt"][map_a["index"]] < .70 or right["plddt"][map_b["index"]] < .70:
                        continue
                    if target_a in ANCESTORS and map_a["pp"] < .80:
                        continue
                    if target_b in ANCESTORS and map_b["pp"] < .80:
                        continue
                    positions.append(human)
                value = rmsd(
                    [left["residues"][crosswalk[target_a][h]["index"]]["CA"].coord for h in positions],
                    [right["residues"][crosswalk[target_b][h]["index"]]["CA"].coord for h in positions],
                )
                values.append(value)
                counts.append(len(positions))
            threshold = 2 * max(within_max[target_a][domain], within_max[target_b][domain])
            both_gated = target_gate[(target_a, domain)] and target_gate[(target_b, domain)]
            finite_values = [value for value in values if np.isfinite(value)]
            robust = both_gated and len(finite_values) == 3 and all(value > threshold for value in finite_values)
            mean_value = float(np.mean(finite_values)) if finite_values else float("nan")
            if both_gated:
                matrices[domain][left_index, right_index] = matrices[domain][right_index, left_index] = mean_value
            contrast_rows.append({
                "target_a": target_a, "lineage_a": panel[target_a]["lineage_group"],
                "target_b": target_b, "lineage_b": panel[target_b]["lineage_group"], "domain": domain,
                "both_target_gates_pass": str(both_gated), "matched_seed_rmsd_mean_A": f"{mean_value:.6f}",
                "matched_seed_rmsd_min_A": f"{min(finite_values):.6f}" if finite_values else "nan",
                "matched_seed_rmsd_max_A": f"{max(finite_values):.6f}" if finite_values else "nan",
                "mapped_high_confidence_n_min": min(counts), "two_x_larger_within_seed_max_A": f"{threshold:.6f}",
                "all_three_between_exceed_threshold": str(robust),
            })
    write_tsv(OUT / "between_target_domain_contrasts.tsv", contrast_rows)

    lineage_rows = []
    groups = []
    for row in panel_rows:
        if row["lineage_group"] not in groups:
            groups.append(row["lineage_group"])
    for domain in DOMAINS:
        for group_a_index, group_a in enumerate(groups):
            for group_b in groups[group_a_index:]:
                selected = []
                robust_n = 0
                for row in contrast_rows:
                    if row["domain"] != domain or row["both_target_gates_pass"] != "True":
                        continue
                    pair = {row["lineage_a"], row["lineage_b"]}
                    wanted = {group_a, group_b}
                    if pair != wanted:
                        continue
                    if group_a == group_b and not (row["lineage_a"] == group_a and row["lineage_b"] == group_b):
                        continue
                    selected.append(float(row["matched_seed_rmsd_mean_A"]))
                    robust_n += row["all_three_between_exceed_threshold"] == "True"
                if selected:
                    lineage_rows.append({
                        "domain": domain, "lineage_a": group_a, "lineage_b": group_b, "target_pair_n": len(selected),
                        "mean_pair_rmsd_A": f"{np.mean(selected):.6f}", "median_pair_rmsd_A": f"{np.median(selected):.6f}",
                        "minimum_pair_rmsd_A": f"{min(selected):.6f}", "maximum_pair_rmsd_A": f"{max(selected):.6f}",
                        "robust_pair_n": robust_n,
                    })
    write_tsv(OUT / "lineage_pair_summary.tsv", lineage_rows)

    OUT.mkdir(parents=True, exist_ok=True)
    for domain, matrix in matrices.items():
        np.savetxt(OUT / f"{domain}_matched_seed_mean_rmsd_matrix.tsv", matrix, delimiter="\t", fmt="%.6f")
        plot_heatmap(matrix, targets, [panel[target]["lineage_group"] for target in targets], f"DNMT3L {domain} structural distances", OUT / f"{domain}_structural_distance_heatmap.png")
    print(f"Wrote full monomer analysis to {OUT}")


if __name__ == "__main__":
    main()
