#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results"
OUTDIR = PROJECT / "figures/current_v12"
LINEAGES = ["Placentalia", "Marsupialia", "Sauropsida", "Amphibia"]
POSITIONS = [88, 90, 105, 106, 107, 204, 349, 351, 352, 355, 358, 377]
FEATURES = {
    88: "H3K4ME0_D88",
    90: "H3K4ME0_D90_D124",
    105: "G9A_INTERFACE_L139",
    106: "G9A_INTERFACE_F140",
    107: "H3K4ME0_I107",
    204: "PTM_STABILITY_K238",
    349: "SWITCH_Q349",
    351: "SWITCH_K351",
    352: "SWITCH_Q352",
    355: "SWITCH_K355",
    358: "H358_CANDIDATE",
    377: "PTM_STABILITY_K412",
}


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.16,
        1.12,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        zorder=10,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.4},
    )


def selected_rows(atlas: pd.DataFrame) -> pd.DataFrame:
    parts = []
    for position, feature in FEATURES.items():
        part = atlas.loc[(atlas["human_np_position"] == position) & (atlas["feature_id"] == feature)].copy()
        if part.empty:
            raise ValueError(f"Missing atlas feature {feature} at NP{position}")
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def summarize_pass(atlas: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (position, lineage), group in atlas.groupby(["human_np_position", "broad_lineage"], sort=True):
        calls = group.loc[group["mapping_status"].eq("PASS"), "query_residue"].dropna().astype(str)
        counts = Counter(calls)
        human = str(group["human_np_residue"].iloc[0])
        rows.append(
            {
                "human_np_position": int(position),
                "broad_lineage": lineage,
                "records": len(group),
                "pass_records": len(calls),
                "pass_fraction": len(calls) / len(group),
                "human_state": human,
                "human_state_count_pass": counts.get(human, 0),
                "human_state_frequency_pass": counts.get(human, 0) / len(calls) if len(calls) else np.nan,
                "residue_counts_pass": ";".join(f"{state}:{count}" for state, count in sorted(counts.items())),
            }
        )
    return pd.DataFrame(rows)


def count_maps(atlas: pd.DataFrame, position: int) -> tuple[list[Counter], np.ndarray, np.ndarray]:
    maps: list[Counter] = []
    passed = []
    records = []
    for lineage in LINEAGES:
        group = atlas.loc[(atlas["human_np_position"] == position) & (atlas["broad_lineage"] == lineage)]
        calls = group.loc[group["mapping_status"].eq("PASS"), "query_residue"].dropna().astype(str)
        maps.append(Counter(calls))
        passed.append(len(calls))
        records.append(len(group))
    return maps, np.asarray(passed, dtype=int), np.asarray(records, dtype=int)


def main() -> None:
    atlas = pd.read_csv(RESULTS / "functional_site_states_by_record.tsv", sep="\t")
    atlas = selected_rows(atlas)
    atlas = atlas.loc[atlas["broad_lineage"].isin(LINEAGES)].copy()
    summary = summarize_pass(atlas)
    summary.to_csv(RESULTS / "functional_site_summary_by_lineage_pass_only.tsv", sep="\t", index=False)

    pivot = summary.pivot(index="broad_lineage", columns="human_np_position", values="human_state_frequency_pass")
    pivot = pivot.reindex(index=LINEAGES, columns=POSITIONS)
    n_pass = summary.pivot(index="broad_lineage", columns="human_np_position", values="pass_records").reindex(index=LINEAGES, columns=POSITIONS)
    n_total = summary.pivot(index="broad_lineage", columns="human_np_position", values="records").reindex(index=LINEAGES, columns=POSITIONS)

    plt.rcParams.update({"font.family": "Arial", "axes.edgecolor": "#334155", "text.color": "#111827"})
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 7.2))

    ax = axes[0, 0]
    cmap = LinearSegmentedColormap.from_list("white_blue", ["#F8FAFC", "#93C5FD", "#1D4ED8"])
    image = ax.imshow(pivot.to_numpy(), aspect="auto", vmin=0, vmax=1, cmap=cmap)
    ax.set_xticks(range(len(POSITIONS)), [str(position) for position in POSITIONS], rotation=45, ha="right")
    ax.set_yticks(range(len(LINEAGES)), ["Eutheria", "Marsupialia", "Sauropsida", "Amphibia"])
    ax.set_xlabel("Human-referenced position")
    ax.set_title("Human state among accepted mappings")
    for row in range(pivot.shape[0]):
        for column in range(pivot.shape[1]):
            value = pivot.iloc[row, column]
            passed = int(n_pass.iloc[row, column])
            total = int(n_total.iloc[row, column])
            if pd.notna(value):
                color = "white" if value > 0.65 else "#0F172A"
                ax.text(column, row, f"{value:.2f}\n{passed}/{total}", ha="center", va="center", fontsize=4.7, color=color)
            else:
                ax.text(column, row, f"--\n0/{total}", ha="center", va="center", fontsize=4.7, color="#64748B")
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Frequency")
    panel_label(ax, "a")

    residue_colors = {
        "Q": "#2563EB",
        "K": "#DC2626",
        "E": "#F59E0B",
        "D": "#7C3AED",
        "L": "#0F766E",
        "C": "#64748B",
        "other": "#CBD5E1",
    }
    for axis, position, title, residues in [
        (axes[0, 1], 349, "NP349 states among accepted mappings", ["Q", "K", "E", "D"]),
        (axes[1, 0], 352, "NP352 states among accepted mappings", ["Q", "L", "C"]),
    ]:
        maps, passed, records = count_maps(atlas, position)
        bottom = np.zeros(len(LINEAGES))
        for residue in residues:
            values = np.asarray([counts.get(residue, 0) for counts in maps], dtype=float)
            fractions = np.divide(values, passed, out=np.zeros_like(values), where=passed > 0)
            axis.bar(range(len(LINEAGES)), fractions, bottom=bottom, color=residue_colors[residue], label=residue)
            bottom += fractions
        other = np.where(passed > 0, 1 - bottom, 0)
        axis.bar(range(len(LINEAGES)), other, bottom=bottom, color=residue_colors["other"], label="other")
        for index, (pcount, total) in enumerate(zip(passed, records)):
            axis.text(index, 1.02, f"accepted {pcount}/{total}", ha="center", va="bottom", fontsize=5.8, color="#475569")
        axis.set_xticks(range(len(LINEAGES)), ["Eutherian", "Marsupial", "Sauropsid", "Amphibian"], rotation=25, ha="right")
        axis.set_ylim(0, 1.13)
        axis.set_ylabel("Fraction")
        axis.set_title(title, pad=25)
        axis.legend(frameon=False, fontsize=6.1, ncol=len(residues) + 1, loc="lower center", bbox_to_anchor=(0.5, 1.0), borderaxespad=0)
    panel_label(axes[0, 1], "b")
    panel_label(axes[1, 0], "c")

    ax = axes[1, 1]
    focal = summary.loc[summary["human_np_position"].isin([349, 351, 352, 355])]
    mean_pass = focal.groupby("broad_lineage")["pass_fraction"].mean().reindex(LINEAGES)
    regions = pd.read_csv(RESULTS / "functional_region_summary_by_lineage.tsv", sep="\t").set_index("broad_lineage").reindex(LINEAGES)
    x = np.arange(len(LINEAGES))
    width = 0.34
    ax.bar(x - width / 2, mean_pass, width, color="#7C3AED", label="Mean accepted-mapping rate at NP349/351/352/355")
    ax.bar(x + width / 2, regions["coverage_pass_fraction"], width, color="#0F766E", label="15-residue region coverage")
    ax.set_xticks(x, ["Eutherian", "Marsupial", "Sauropsid", "Amphibian"], rotation=25, ha="right")
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Fraction of records")
    ax.set_title("Local-context criterion versus region coverage", pad=9)
    ax.legend(
        frameon=False,
        fontsize=6.2,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.30),
        ncol=1,
        borderaxespad=0,
    )
    panel_label(ax, "d")

    fig.suptitle("Experimentally defined DNMT3L site-state atlas", y=1.01, fontsize=11, fontweight="bold")
    fig.tight_layout(h_pad=2.4, w_pad=1.8, rect=[0, 0.05, 1, 1])
    OUTDIR.mkdir(parents=True, exist_ok=True)
    png = OUTDIR / "Figure_5_functional_site_atlas_accepted_mappings.png"
    pdf = OUTDIR / "Figure_5_functional_site_atlas_accepted_mappings.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
