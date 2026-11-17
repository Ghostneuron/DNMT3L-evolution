#!/usr/bin/env python3
"""Build the molecular-dynamics figure with node-state terminology."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "analysis/therian_ff_simulation_20261104/analysis_summary/replica_metrics.tsv"
OUT = PROJECT / "figures/current_v08"

PANELS = [
    ("min_B261_A728_nm", "DNMT3L-261 to DNMT3A1-F732\n(2QRV F728) distance (nm)"),
    ("min_B291_A725_R729ref_nm", "DNMT3L-291 to DNMT3A1-R729\n(2QRV R725) distance (nm)"),
    ("min_B293_A722_G726ref_nm", "DNMT3L-293 to DNMT3A1-G726\n(2QRV G722) distance (nm)"),
    ("active_loop_mean_ca_rmsf_nm", "Adjacent DNMT3A-loop\nmean CA RMSF (nm)"),
]


def main() -> None:
    data = pd.read_csv(INPUT, sep="\t")
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2), constrained_layout=True)
    variants = [
        ("amniote_YEE", "Amniota node\nY/E/E", "#2563EB"),
        ("therian_FTD", "Therian node\nF/T/D", "#D97706"),
    ]
    jitter = np.array([-0.09, 0.0, 0.09])

    for index, (axis, (metric, ylabel)) in enumerate(zip(axes.flat, PANELS)):
        upper = 0.0
        for x, (variant, label, color) in enumerate(variants):
            values = data.loc[data["variant"].eq(variant), metric].to_numpy()
            axis.scatter(
                x + jitter[: len(values)],
                values,
                s=42,
                color=color,
                edgecolor="white",
                linewidth=0.5,
                zorder=3,
            )
            axis.errorbar(
                x,
                values.mean(),
                yerr=values.std(ddof=1),
                fmt="_",
                markersize=18,
                color="#111827",
                capsize=4,
                linewidth=1.2,
                zorder=4,
            )
            upper = max(upper, float(values.max()))
        axis.set_xticks([0, 1], [item[1] for item in variants])
        axis.set_ylabel(ylabel)
        axis.set_ylim(0, max(upper * 1.16, 0.7 if metric == "active_loop_mean_ca_rmsf_nm" else 0.55))
        if metric.startswith("min_"):
            axis.axhline(0.5, color="#64748B", linestyle="--", linewidth=0.8)
            axis.text(
                1.03,
                0.5,
                "0.5-nm contact threshold",
                fontsize=6.5,
                color="#475569",
                va="bottom",
                ha="right",
            )
        axis.grid(axis="y", color="#E5E7EB", linewidth=0.7)
        axis.set_title(chr(97 + index), loc="left", fontweight="bold", fontsize=11)

    fig.suptitle(
        "Independent 5-ns replicas of reconstructed DNMT3L node-state combinations",
        fontsize=11,
        fontweight="bold",
    )
    stem = OUT / "Figure_6_node_state_interface_md"
    fig.savefig(stem.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(stem)


if __name__ == "__main__":
    main()
