#!/usr/bin/env python3
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results"
OUT = PROJECT / "figures"
METHOD_ORDER = ["CONSENSUS", "MACSE", "MAFFT", "PRANK", "MUSCLE5", "COBALT"]
METHOD_COLORS = {
    "CONSENSUS": "#1F2937",
    "MACSE": "#2563EB",
    "MAFFT": "#0F766E",
    "PRANK": "#7C3AED",
    "MUSCLE5": "#D97706",
    "COBALT": "#DC2626",
}
LINEAGE_ORDER = ["Placentalia", "Marsupialia", "Sauropsida", "Amphibia"]


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.2,
            "axes.titlesize": 9.2,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def panel_label(ax, label: str, x: float = -0.17, y: float = 1.06) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
        zorder=10,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.4},
    )


def save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.png", dpi=350, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def parse_counts(value: str) -> dict[str, int]:
    result: dict[str, int] = {}
    if not isinstance(value, str) or not value:
        return result
    for part in value.split(";"):
        residue, count = part.split(":")
        result[residue] = int(count)
    return result


def figure_1() -> None:
    audit = pd.read_csv(RESULTS / "ncbi_337_package_audit_summary.tsv", sep="\t")
    audit_values = dict(zip(audit["metric"], audit["value"]))
    representative = pd.read_csv(
        RESULTS / "ncbi_337_representative_selection_summary.tsv", sep="\t"
    )
    representative_values = dict(zip(representative["metric"], representative["value"]))
    family = pd.read_csv(
        RESULTS / "DNMT3_family_337_group_separation_summary.tsv", sep="\t"
    )
    aligners = pd.read_csv(RESULTS / "selection_337_input_summary.tsv", sep="\t")
    mask = pd.read_csv(
        RESULTS / "placental_five_aligner_consensus_mask_summary.tsv", sep="\t"
    )
    mask_values = dict(zip(mask["metric"], mask["value"]))

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2))

    ax = axes[0, 0]
    stages = ["NCBI\nGeneIDs", "Exact\npairs", "Sequence\nQC", "Domain\nQC", "Analysis\nready"]
    values = [
        int(audit_values["unique_gene_ids"]),
        int(audit_values["genes_with_exact_protein_cds_pair"]),
        int(representative_values["Sequence_QC_pass"]),
        int(representative_values["Domain_QC_pass"]),
        int(representative_values["Preliminary_alignment_ready"]),
    ]
    bars = ax.bar(range(len(values)), values, color=["#334155", "#475569", "#2563EB", "#0F766E", "#7C3AED"])
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 7, str(value), ha="center", fontweight="bold")
    ax.set_xticks(range(len(values)), stages)
    ax.set_ylim(0, 370)
    ax.set_ylabel("Records")
    ax.set_title("Independent criteria and their intersection")
    panel_label(ax, "a")

    ax = axes[0, 1]
    x = np.arange(len(family))
    bars = ax.bar(x, family["SH_aLRT_support"], color=["#0F766E", "#2563EB", "#7C3AED", "#334155"])
    ax.axhline(95, color="#64748B", lw=0.8, ls="--")
    family_labels = ["DNMT3L", "DNMT3A", "DNMT3B", "DNMT3A+B"]
    ax.set_xticks(x, family_labels)
    ax.set_ylim(0, 105)
    ax.set_ylabel("SH-aLRT support")
    ax.set_title("Exact root-independent family splits")
    for bar, tips in zip(bars, family["n_tips"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.2, f"n={tips}", ha="center", fontsize=7)
    panel_label(ax, "b")

    ax = axes[1, 0]
    aligners = aligners.set_index("analysis_id").loc[METHOD_ORDER[1:]].reset_index()
    x = np.arange(len(aligners))
    ax.bar(
        x,
        aligners["retained_human_positions"],
        color=[METHOD_COLORS[m] for m in aligners["analysis_id"]],
    )
    ax.set_xticks(x, aligners["analysis_id"], rotation=25, ha="right")
    ax.set_ylim(0, 400)
    ax.set_ylabel("Positions retained at ≥80% occupancy")
    ax.set_title("Retained positions differ among workflows")
    for index, row in aligners.iterrows():
        ax.text(index, row["retained_human_positions"] + 6, str(int(row["retained_human_positions"])), ha="center", fontsize=7)
    panel_label(ax, "c", y=1.20)

    ax = axes[1, 1]
    retained = int(float(mask_values["retained_human_positions"]))
    removed = int(float(mask_values["removed_low_occupancy_positions"]))
    masked = int(float(mask_values["taxon_position_calls_masked_in_all_methods"]))
    total = int(float(mask_values["retained_taxon_position_calls"]))
    wedges, _ = ax.pie(
        [retained, removed],
        colors=["#2563EB", "#E2E8F0"],
        startangle=90,
        counterclock=False,
        wedgeprops={"width": 0.34, "edgecolor": "white"},
    )
    ax.text(0, 0.08, f"{retained}/387", ha="center", va="center", fontsize=12, fontweight="bold")
    ax.text(0, -0.14, "shared positions", ha="center", va="center", fontsize=7.5)
    ax.text(
        0,
        -1.25,
        f"{masked:,}/{total:,} retained taxon-position calls masked\nfor missing or discordant exact codons ({100*masked/total:.2f}%)",
        ha="center",
        fontsize=7.3,
    )
    ax.set_title("Shared occupancy and exact-codon criterion", pad=12)
    panel_label(ax, "d", x=-0.20, y=1.25)

    fig.suptitle("DNMT3L data foundation and alignment audit", y=1.01, fontsize=11, fontweight="bold")
    fig.tight_layout(h_pad=2.2, w_pad=1.6, rect=[0.03, 0, 1, 1])
    save(fig, "Figure_1_updated_data_foundation")


def figure_2() -> None:
    sites = pd.read_csv(
        RESULTS / "placental_five_aligner_constraint_by_site.tsv", sep="\t"
    )
    tests = pd.read_csv(
        RESULTS / "placental_five_aligner_constraint_tests.tsv", sep="\t"
    )
    by_method = pd.read_csv(
        RESULTS / "placental_five_aligner_constraint_by_method_site.tsv", sep="\t"
    )

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2))

    ax = axes[0, 0]
    methods = METHOD_ORDER[1:]
    for index, method in enumerate(methods):
        row = tests.loc[
            tests["analysis"].eq(method)
            & tests["comparison"].eq("ADD_vs_C_terminal_MTase_like")
        ].iloc[0]
        ax.plot(
            [0, 1],
            [row["first_median_entropy"], row["second_median_entropy"]],
            marker="o",
            color=METHOD_COLORS[method],
            lw=1.2,
            label=method,
        )
    ax.set_xticks([0, 1], ["ADD", "C-terminal"])
    ax.set_ylabel("Median amino-acid entropy")
    ax.set_title("All five aligners show lower ADD entropy")
    ax.legend(frameon=False, fontsize=6.8, ncol=2)
    panel_label(ax, "a")

    ax = axes[0, 1]
    ax.scatter(sites["human_position"], sites["median_entropy"], s=8, color="#475569", alpha=0.75)
    ax.axvspan(41, 173, color="#BFDBFE", alpha=0.65, label="ADD")
    ax.axvspan(178, 380, color="#E2E8F0", alpha=0.45, label="C-terminal")
    ax.axvspan(341, 355, color="#F59E0B", alpha=0.28, label="switching helix")
    for position, label in [(289, "NP289"), (347, "NP347"), (358, "NP358")]:
        point = sites.loc[sites["human_position"].eq(position)]
        if not point.empty:
            value = float(point.iloc[0]["median_entropy"])
            ax.annotate(label, (position, value), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=7)
    ax.set_xlim(1, 387)
    ax.set_xlabel("Human NP_037501.2 position")
    ax.set_ylabel("Median entropy across aligners")
    ax.set_title("Position-level amino-acid variability")
    ax.legend(frameon=False, fontsize=6.4, ncol=2, loc="upper left")
    panel_label(ax, "b")

    ax = axes[1, 0]
    add = sites.loc[sites["domain"].eq("ADD"), "median_entropy"]
    cterminal = sites.loc[sites["domain"].eq("C_terminal_MTase_like"), "median_entropy"]
    box = ax.boxplot(
        [add, cterminal],
        tick_labels=[f"ADD\n(n={len(add)})", f"C-terminal\n(n={len(cterminal)})"],
        patch_artist=True,
        showfliers=False,
        widths=0.55,
    )
    for patch, color in zip(box["boxes"], ["#93C5FD", "#CBD5E1"]):
        patch.set_facecolor(color)
    row = tests.loc[
        tests["analysis"].eq("FIVE_ALIGNER_MEDIAN")
        & tests["comparison"].eq("ADD_vs_C_terminal_MTase_like")
    ].iloc[0]
    ax.text(
        0.5,
        0.91,
        f"Cliff's delta = {row['cliffs_delta']:.3f}; p = {row['p']:.2e}",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=7.4,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.5},
    )
    ax.set_ylabel("Median amino-acid entropy")
    ax.set_title("Eutherian domain contrast")
    panel_label(ax, "c")

    ax = axes[1, 1]
    background = sites.loc[
        sites["domain"].eq("C_terminal_MTase_like")
        & sites["structural_region"].eq("other"),
        "median_entropy",
    ]
    helix = sites.loc[sites["structural_region"].eq("switching_helix"), "median_entropy"]
    loop = sites.loc[sites["structural_region"].eq("adjacent_loop"), "median_entropy"]
    box = ax.boxplot(
        [background, helix, loop],
        tick_labels=[f"Other C-terminal\n(n={len(background)})", f"Helix\n(n={len(helix)})", f"Loop\n(n={len(loop)})"],
        patch_artist=True,
        showfliers=False,
        widths=0.55,
    )
    for patch, color in zip(box["boxes"], ["#CBD5E1", "#FBBF24", "#FDBA74"]):
        patch.set_facecolor(color)
    ax.set_ylabel("Median amino-acid entropy")
    ax.set_title("Local variability within the C-terminal region")
    panel_label(ax, "d")

    fig.suptitle("Five-aligner eutherian amino-acid variability", y=1.01, fontsize=11, fontweight="bold")
    fig.tight_layout(h_pad=2.0, w_pad=1.6)
    save(fig, "Figure_2_five_aligner_variability")


def figure_3() -> None:
    sites = pd.read_csv(RESULTS / "functional_site_summary_by_lineage.tsv", sep="\t")
    regions = pd.read_csv(RESULTS / "functional_region_summary_by_lineage.tsv", sep="\t")
    sites = sites.loc[sites["broad_lineage"].isin(LINEAGE_ORDER)].copy()
    positions = [88, 90, 105, 106, 107, 204, 349, 351, 352, 355, 358, 377]
    unique = sites.sort_values(["human_np_position", "atlas_role"]).drop_duplicates(
        ["human_np_position", "broad_lineage"]
    )
    pivot = unique.pivot(
        index="broad_lineage",
        columns="human_np_position",
        values="human_state_frequency_all_mapped",
    ).reindex(index=LINEAGE_ORDER, columns=positions)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.3))

    ax = axes[0, 0]
    cmap = LinearSegmentedColormap.from_list("white_blue", ["#F8FAFC", "#93C5FD", "#1D4ED8"])
    image = ax.imshow(pivot.to_numpy(), aspect="auto", vmin=0, vmax=1, cmap=cmap)
    ax.set_xticks(range(len(positions)), [str(position) for position in positions], rotation=45, ha="right")
    ax.set_yticks(range(len(LINEAGE_ORDER)), LINEAGE_ORDER)
    ax.set_xlabel("Human-referenced position")
    ax.set_title("Frequency of the human residue state")
    for row in range(pivot.shape[0]):
        for column in range(pivot.shape[1]):
            value = pivot.iloc[row, column]
            if pd.notna(value):
                ax.text(column, row, f"{value:.2f}", ha="center", va="center", fontsize=5.7, color="white" if value > 0.65 else "#0F172A")
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
    for axis, feature, title, residues in [
        (axes[0, 1], "SWITCH_Q349", "NP349 residue states", ["Q", "K", "E", "D"]),
        (axes[1, 0], "SWITCH_Q352", "NP352 residue states", ["Q", "L", "C"]),
    ]:
        subset = sites.loc[sites["feature_id"].eq(feature)].set_index("broad_lineage").reindex(LINEAGE_ORDER)
        bottom = np.zeros(len(LINEAGE_ORDER))
        count_maps = [parse_counts(value) for value in subset["residue_counts_all_mapped"]]
        for residue in residues:
            values = np.array([counts.get(residue, 0) for counts in count_maps], dtype=float)
            totals = np.array([sum(counts.values()) for counts in count_maps], dtype=float)
            fractions = np.divide(values, totals, out=np.zeros_like(values), where=totals > 0)
            axis.bar(range(len(LINEAGE_ORDER)), fractions, bottom=bottom, color=residue_colors[residue], label=residue)
            bottom += fractions
        other = 1 - bottom
        axis.bar(range(len(LINEAGE_ORDER)), other, bottom=bottom, color=residue_colors["other"], label="other")
        axis.set_xticks(range(len(LINEAGE_ORDER)), ["Placental", "Marsupial", "Sauropsid", "Amphibian"], rotation=25, ha="right")
        axis.set_ylim(0, 1)
        axis.set_ylabel("Resolved-state fraction")
        axis.set_title(title, pad=24)
        axis.legend(
            frameon=False,
            fontsize=6.2,
            ncol=len(residues) + 1,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.0),
            borderaxespad=0,
        )
    panel_label(axes[0, 1], "b")
    panel_label(axes[1, 0], "c")

    ax = axes[1, 1]
    regions = regions.set_index("broad_lineage").reindex(LINEAGE_ORDER)
    x = np.arange(len(LINEAGE_ORDER))
    width = 0.34
    ax.bar(
        x - width / 2,
        regions["median_identity_fraction_resolved"],
        width,
        color="#0F766E",
        label="Identity to human helix",
    )
    ax.bar(
        x + width / 2,
        regions["median_polar_or_charged_fraction"],
        width,
        color="#F59E0B",
        label="Polar or charged fraction",
    )
    ax.set_xticks(x, ["Placental", "Marsupial", "Sauropsid", "Amphibian"], rotation=25, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Median fraction")
    ax.set_title("Switching-helix sequence properties")
    ax.legend(frameon=False, fontsize=6.5)
    panel_label(ax, "d")

    fig.suptitle("Experimentally defined DNMT3L site-state atlas", y=1.01, fontsize=11, fontweight="bold")
    fig.tight_layout(h_pad=2.0, w_pad=1.6)
    save(fig, "Figure_3_functional_site_atlas")


def figure_4() -> None:
    analyses = pd.read_csv(RESULTS / "selection_337_analysis_summary.tsv", sep="\t").set_index("analysis_id").loc[METHOD_ORDER].reset_index()
    all_sites = pd.read_csv(RESULTS / "selection_337_all_site_summary.tsv", sep="\t")
    fubar = pd.read_csv(RESULTS / "selection_337_fubar_candidate_decisions.tsv", sep="\t")

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2))

    ax = axes[0, 0]
    x = np.arange(len(analyses))
    ax.bar(x, analyses["MEME_BH_q_lt_0.05"], color=[METHOD_COLORS[m] for m in analyses["analysis_id"]])
    ax.set_xticks(x, analyses["analysis_id"], rotation=25, ha="right")
    ax.set_ylabel("Sites with MEME BH q < 0.05")
    ax.set_title("Episodic calls vary by alignment")
    for index, value in enumerate(analyses["MEME_BH_q_lt_0.05"]):
        ax.text(index, value + 0.4, str(int(value)), ha="center", fontsize=7)
    panel_label(ax, "a")

    ax = axes[0, 1]
    h358 = all_sites.loc[all_sites["human_position"].eq(358)].set_index("analysis_id").loc[METHOD_ORDER]
    x = np.arange(len(h358))
    bars = ax.bar(x, h358["MEME_BH_q"], color=[METHOD_COLORS[m] for m in h358.index])
    ax.axhline(0.05, color="#B91C1C", lw=1, ls="--", label="BH q = 0.05")
    ax.set_xticks(x, h358.index, rotation=25, ha="right")
    ax.set_ylim(0, max(0.72, float(h358["MEME_BH_q"].max()) + 0.05))
    ax.set_ylabel("NP358 MEME BH q")
    ax.set_title("NP358 is alignment-sensitive")
    ax.legend(frameon=False, fontsize=6.8)
    panel_label(ax, "b")

    ax = axes[1, 0]
    focus = fubar.loc[fubar["human_position"].isin([251, 289, 347, 358])].copy()
    x = np.arange(len(focus))
    width = 0.35
    ax.bar(x - width / 2, focus["gene_tree_consensus_posterior"], width, color="#2563EB", label="Codon gene tree")
    ax.bar(x + width / 2, focus["ncbi_taxonomy_posterior"], width, color="#F59E0B", label="NCBI taxonomy")
    ax.axhline(0.90, color="#64748B", lw=0.9, ls="--")
    ax.set_xticks(x, [f"NP{int(p)}" for p in focus["human_position"]])
    ax.set_ylim(0.65, 1.01)
    ax.set_ylabel("FUBAR diversifying posterior")
    ax.set_title("Topology sensitivity of pervasive candidates")
    ax.legend(frameon=False, fontsize=6.7)
    panel_label(ax, "c")

    ax = axes[1, 1]
    ax.hlines(0, 1, 387, color="#475569", lw=4)
    ax.axvspan(41, 173, ymin=0.40, ymax=0.60, color="#93C5FD")
    ax.axvspan(178, 380, ymin=0.40, ymax=0.60, color="#CBD5E1")
    ax.axvspan(341, 355, ymin=0.30, ymax=0.70, color="#F59E0B", alpha=0.8)
    points = [(289, "NP289V", "#2563EB", 0.35), (347, "NP347L", "#0F766E", -0.35), (358, "NP358A", "#DC2626", 0.35)]
    for position, label, color, y in points:
        ax.vlines(position, 0, y, color=color, lw=1.4)
        ax.scatter([position], [y], s=30, color=color, zorder=3)
        ax.text(position, y + (0.09 if y > 0 else -0.09), label, ha="center", va="bottom" if y > 0 else "top", fontsize=7)
    ax.text(107, -0.08, "ADD", ha="center", va="top", fontsize=7)
    ax.text(260, -0.08, "C-terminal", ha="center", va="top", fontsize=7)
    ax.text(348, 0.08, "switching helix", ha="center", va="bottom", fontsize=6.5)
    ax.set_xlim(1, 387)
    ax.set_ylim(-0.65, 0.65)
    ax.set_yticks([])
    ax.set_xlabel("Human NP_037501.2 position")
    ax.set_title("Candidate locations and evidentiary status")
    panel_label(ax, "d")

    fig.suptitle("Selection-model sensitivity in the eutherian panel", y=1.01, fontsize=11, fontweight="bold")
    fig.tight_layout(h_pad=2.0, w_pad=1.6)
    save(fig, "Figure_4_selection_sensitivity")


def main() -> None:
    style()
    figure_1()
    figure_2()
    figure_3()
    figure_4()
    print(OUT)


if __name__ == "__main__":
    main()
