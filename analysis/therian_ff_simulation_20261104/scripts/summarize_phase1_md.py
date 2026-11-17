#!/usr/bin/env python3
"""Validate and summarize the six prespecified phase-1 OpenMM replicas."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPLICAS = [
    ("amniote_YEE", "Y/E/E", 1, 1101),
    ("amniote_YEE", "Y/E/E", 2, 1102),
    ("amniote_YEE", "Y/E/E", 3, 1103),
    ("therian_FTD", "F/T/D", 1, 2101),
    ("therian_FTD", "F/T/D", 2, 2102),
    ("therian_FTD", "F/T/D", 3, 2103),
]

METRICS = {
    "min_B261_A728_nm": ("B261-A728 minimum distance", "nm"),
    "min_B261_A768_nm": ("B261-A768 minimum distance", "nm"),
    "min_B291_A725_R729ref_nm": ("B291-A725 (DNMT3A1 R729) minimum distance", "nm"),
    "min_B291_A721_E725ref_nm": ("B291-A721 (DNMT3A1 E725) minimum distance", "nm"),
    "min_B293_A725_R729ref_nm": ("B293-A725 (DNMT3A1 R729) minimum distance", "nm"),
    "min_B293_A721_E725ref_nm": ("B293-A721 (DNMT3A1 E725) minimum distance", "nm"),
    "min_B293_A722_G726ref_nm": ("B293-A722 (DNMT3A1 G726) minimum distance", "nm"),
    "protein_backbone_rmsd_nm": ("Protein-backbone RMSD", "nm"),
    "contact05_B261_A728_nm": ("B261-A728 contact occupancy (<0.5 nm)", "fraction"),
    "contact05_B261_A768_nm": ("B261-A768 contact occupancy (<0.5 nm)", "fraction"),
    "contact05_B291_A725_R729ref_nm": ("B291-A725 contact occupancy (<0.5 nm)", "fraction"),
    "contact05_B291_A721_E725ref_nm": ("B291-A721 contact occupancy (<0.5 nm)", "fraction"),
    "contact05_B293_A725_R729ref_nm": ("B293-A725 contact occupancy (<0.5 nm)", "fraction"),
    "contact05_B293_A721_E725ref_nm": ("B293-A721 contact occupancy (<0.5 nm)", "fraction"),
    "contact05_B293_A722_G726ref_nm": ("B293-A722 contact occupancy (<0.5 nm)", "fraction"),
    "active_loop_mean_ca_rmsf_nm": ("DNMT3A active-loop mean CA RMSF", "nm"),
}

TABLE_METRICS = [
    "min_B261_A728_nm",
    "contact05_B261_A728_nm",
    "min_B261_A768_nm",
    "min_B291_A725_R729ref_nm",
    "min_B293_A725_R729ref_nm",
    "min_B293_A722_G726ref_nm",
    "protein_backbone_rmsd_nm",
    "active_loop_mean_ca_rmsf_nm",
]


def load_replica(root: Path, variant: str, state: str, pair: int, seed: int) -> dict:
    directory = root / f"{variant}_seed{seed}"
    completion = json.loads((directory / "COMPLETE.json").read_text())
    analysis = directory / "analysis"
    metrics = pd.read_csv(analysis / "frame_metrics.tsv", sep="\t")
    loop = pd.read_csv(analysis / "active_loop_rmsf.tsv", sep="\t")
    summary = json.loads((analysis / "summary.json").read_text())

    expected_atoms = 142702 if variant == "amniote_YEE" else 142689
    checks = {
        "final_step": completion.get("final_step") == 1_275_000,
        "production_time": completion.get("production_time_ns") == 5.0,
        "frames": len(metrics) == 500 and summary.get("frames") == 500,
        "frame_index": metrics["frame"].tolist() == list(range(500)),
        "atoms": summary.get("atoms") == expected_atoms,
        "loop_residues": loop["pdb_residue"].tolist() == list(range(700, 722)),
        "finite": np.isfinite(metrics.select_dtypes(include=[np.number])).all().all()
        and np.isfinite(loop["ca_rmsf_nm"]).all(),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise RuntimeError(f"Validation failed for {directory.name}: {failed}")

    row = {
        "variant": variant,
        "state": state,
        "replicate_pair": pair,
        "seed": seed,
        "frames": len(metrics),
        "production_ns": float(completion["production_time_ns"]),
        "first_time_ps": float(metrics["time_ps"].iloc[0]),
        "last_time_ps": float(metrics["time_ps"].iloc[-1]),
    }
    for column in metrics.columns:
        if column not in {"frame", "time_ps"}:
            row[column] = float(metrics[column].mean())
    row["active_loop_mean_ca_rmsf_nm"] = float(loop["ca_rmsf_nm"].mean())
    return row


def format_value(mean: float, sd: float, unit_name: str) -> str:
    if unit_name == "fraction":
        return f"{100 * mean:.1f} +/- {100 * sd:.1f}%"
    return f"{mean:.3f} +/- {sd:.3f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--production-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--figure-prefix", required=True, type=Path)
    args = parser.parse_args()

    rows = [load_replica(args.production_root, *spec) for spec in REPLICAS]
    replica = pd.DataFrame(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    args.figure_prefix.parent.mkdir(parents=True, exist_ok=True)
    replica.to_csv(args.output / "replica_metrics.tsv", sep="\t", index=False)

    group_rows = []
    paired_rows = []
    for metric, (label, units) in METRICS.items():
        for variant, state in [("amniote_YEE", "Y/E/E"), ("therian_FTD", "F/T/D")]:
            values = replica.loc[replica["variant"].eq(variant), metric].to_numpy()
            group_rows.append(
                {
                    "metric": metric,
                    "label": label,
                    "units": units,
                    "variant": variant,
                    "state": state,
                    "n_replicas": len(values),
                    "mean_of_replica_means": values.mean(),
                    "sd_across_replica_means": values.std(ddof=1),
                    "median_of_replica_means": np.median(values),
                    "minimum_replica_mean": values.min(),
                    "maximum_replica_mean": values.max(),
                }
            )
        amniote = (
            replica.loc[replica["variant"].eq("amniote_YEE")]
            .set_index("replicate_pair")[metric]
            .sort_index()
        )
        therian = (
            replica.loc[replica["variant"].eq("therian_FTD")]
            .set_index("replicate_pair")[metric]
            .sort_index()
        )
        difference = therian - amniote
        paired_rows.append(
            {
                "metric": metric,
                "label": label,
                "units": units,
                "therian_minus_amniote_mean": difference.mean(),
                "pair1_difference": difference.loc[1],
                "pair2_difference": difference.loc[2],
                "pair3_difference": difference.loc[3],
                "therian_lower_pairs": int((difference < 0).sum()),
                "therian_higher_pairs": int((difference > 0).sum()),
            }
        )

    groups = pd.DataFrame(group_rows)
    paired = pd.DataFrame(paired_rows)
    groups.to_csv(args.output / "group_metrics.tsv", sep="\t", index=False)
    paired.to_csv(args.output / "paired_contrasts.tsv", sep="\t", index=False)

    table_rows = []
    for metric in TABLE_METRICS:
        label, units = METRICS[metric]
        subset = groups.loc[groups["metric"].eq(metric)].set_index("variant")
        contrast = paired.loc[paired["metric"].eq(metric)].iloc[0]
        if contrast.therian_lower_pairs == 3:
            direction = "Lower in F/T/D in 3/3 pairs"
        elif contrast.therian_higher_pairs == 3:
            direction = "Higher in F/T/D in 3/3 pairs"
        else:
            direction = (
                f"Mixed ({int(contrast.therian_lower_pairs)}/3 lower; "
                f"{int(contrast.therian_higher_pairs)}/3 higher in F/T/D)"
            )
        table_rows.append(
            {
                "metric": label,
                "units": "% frames" if units == "fraction" else units,
                "amniote_YEE_mean_sd": format_value(
                    subset.loc["amniote_YEE", "mean_of_replica_means"],
                    subset.loc["amniote_YEE", "sd_across_replica_means"],
                    units,
                ),
                "therian_FTD_mean_sd": format_value(
                    subset.loc["therian_FTD", "mean_of_replica_means"],
                    subset.loc["therian_FTD", "sd_across_replica_means"],
                    units,
                ),
                "paired_direction": direction,
            }
        )
    manuscript_table = pd.DataFrame(table_rows)
    manuscript_table.to_csv(args.output / "manuscript_table.tsv", sep="\t", index=False)

    plot_metrics = [
        ("min_B261_A728_nm", "B261-A728 distance (nm)"),
        ("min_B291_A725_R729ref_nm", "B291-R729(ref.) distance (nm)"),
        ("min_B293_A722_G726ref_nm", "B293-G726(ref.) distance (nm)"),
        ("active_loop_mean_ca_rmsf_nm", "Active-loop mean CA RMSF (nm)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2), constrained_layout=True)
    colors = ["#0072B2", "#D55E00", "#009E73"]
    for panel, (axis, (metric, ylabel)) in enumerate(zip(axes.flat, plot_metrics)):
        for pair, color in zip((1, 2, 3), colors):
            y0 = replica.loc[
                replica["variant"].eq("amniote_YEE")
                & replica["replicate_pair"].eq(pair), metric
            ].iloc[0]
            y1 = replica.loc[
                replica["variant"].eq("therian_FTD")
                & replica["replicate_pair"].eq(pair), metric
            ].iloc[0]
            axis.plot([0, 1], [y0, y1], color=color, alpha=0.75, linewidth=1.2)
            axis.scatter([0, 1], [y0, y1], color=color, s=34, zorder=3, label=f"Pair {pair}")
        axis.set_xticks([0, 1], ["Y/E/E", "F/T/D"])
        axis.set_ylabel(ylabel)
        axis.set_title(chr(97 + panel), loc="left", fontweight="bold")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.6)
    axes[0, 0].legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle("Replicated 5-ns DNMT3A-DNMT3L interface simulations", fontsize=11)
    fig.savefig(args.figure_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(args.figure_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    lookup = groups.set_index(["metric", "variant"])
    pair_lookup = paired.set_index("metric")
    summary_text = f"""# Phase-1 molecular-dynamics summary

## Validation

- Six canonical replicas passed the completion and analysis gates (3 Y/E/E and 3 F/T/D).
- Each replica contributed 500 production frames from a 5.0-ns trajectory; the aggregate sampled production time was 30.0 ns.
- Each Y/E/E trajectory contained 142,702 atoms and each F/T/D trajectory contained 142,689 atoms.
- The preserved duplicate-corrupted directory was not admitted.

## Main findings

- The B261-A728 aromatic contact was present in 100% of frames in every replica. Replica-mean distances were {lookup.loc[("min_B261_A728_nm", "amniote_YEE"), "mean_of_replica_means"]:.3f} +/- {lookup.loc[("min_B261_A728_nm", "amniote_YEE"), "sd_across_replica_means"]:.3f} nm for Y/E/E and {lookup.loc[("min_B261_A728_nm", "therian_FTD"), "mean_of_replica_means"]:.3f} +/- {lookup.loc[("min_B261_A728_nm", "therian_FTD"), "sd_across_replica_means"]:.3f} nm for F/T/D (mean +/- SD across three replica means). The paired direction was mixed ({int(pair_lookup.loc["min_B261_A728_nm", "therian_lower_pairs"])}/3 lower in F/T/D).
- B291-A725 and B293-A725 distances varied among replicas, without a common direction across all three pairs. Sub-0.5-nm contacts at these sites were infrequent and replica-specific.
- B293-A722 was shorter in F/T/D in all three matched replica indices, but the mean difference was small ({pair_lookup.loc["min_B293_A722_G726ref_nm", "therian_minus_amniote_mean"]:.3f} nm) relative to between-replica variation, and neither endpoint formed a persistent sub-0.5-nm contact.
- Protein-backbone RMSD averaged {lookup.loc[("protein_backbone_rmsd_nm", "amniote_YEE"), "mean_of_replica_means"]:.3f} +/- {lookup.loc[("protein_backbone_rmsd_nm", "amniote_YEE"), "sd_across_replica_means"]:.3f} nm for Y/E/E and {lookup.loc[("protein_backbone_rmsd_nm", "therian_FTD"), "mean_of_replica_means"]:.3f} +/- {lookup.loc[("protein_backbone_rmsd_nm", "therian_FTD"), "sd_across_replica_means"]:.3f} nm for F/T/D.
- Mean DNMT3A active-loop CA RMSF overlapped ({lookup.loc[("active_loop_mean_ca_rmsf_nm", "amniote_YEE"), "mean_of_replica_means"]:.3f} +/- {lookup.loc[("active_loop_mean_ca_rmsf_nm", "amniote_YEE"), "sd_across_replica_means"]:.3f} nm versus {lookup.loc[("active_loop_mean_ca_rmsf_nm", "therian_FTD"), "mean_of_replica_means"]:.3f} +/- {lookup.loc[("active_loop_mean_ca_rmsf_nm", "therian_FTD"), "sd_across_replica_means"]:.3f} nm), with a mixed paired direction.

## Interpretation boundary

The six short trajectories support stable maintenance of the common B261-A728 aromatic contact and do not show a replica-consistent global-stability or active-loop-mobility shift between the reconstructed states. The B293-A722 distance tendency is exploratory. These simulations do not estimate binding free energy, catalytic activity, imprinting, organismal fitness, or adaptation. Frames are repeated observations within trajectories; the biological comparison uses three independent replica summaries per endpoint.
"""
    (args.output / "SUMMARY.md").write_text(summary_text)

    manifest = {
        "analysis": "DNMT3L_FF_THERIAN_MD_20261104_V1_phase1_summary",
        "canonical_replicas": 6,
        "replicas_per_state": 3,
        "frames_per_replica": 500,
        "total_frames": 3000,
        "production_ns_per_replica": 5.0,
        "total_production_ns": 30.0,
        "validation": "PASS",
        "inference_unit": "replica summary",
        "excluded": ["amniote_YEE_seed1103_duplicate_corrupted_20261111"],
    }
    (args.output / "SUMMARY.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
