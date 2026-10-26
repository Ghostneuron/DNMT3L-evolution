from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle
from PIL import Image


INK = "#20262E"
MUTED = "#687683"
GRID = "#D8DEE4"
TEAL = "#168C88"
CORAL = "#D2644F"
GOLD = "#D5A027"
BLUE = "#3976A8"
LINEAGE_DISPLAY = {"Mammalia": "Sampled Theria", "Sauropsida": "Sauropsida", "Amphibia": "Amphibia"}
RESIDUE_PALETTE = {
    "A": "#8EC07C", "C": "#D6B656", "D": "#D65F5F", "E": "#C94B4B", "F": "#8F79B8",
    "G": "#A8B0B8", "H": "#6F87B8", "I": "#4F9D8D", "K": "#3C78A8", "L": "#2E8B78",
    "M": "#448F83", "N": "#D89C48", "P": "#B8739B", "Q": "#E2A83B", "R": "#245E91",
    "S": "#73A9C2", "T": "#7AB7B2", "V": "#60A37B", "W": "#6F5B9B", "Y": "#7B67A7", "-": "#E4E8EB",
}


def configure() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 10.0,
            "axes.labelsize": 8.5,
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "legend.fontsize": 7.0,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def panel_label(ax, label: str) -> None:
    ax.text(-0.08, 1.05, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")


def parse_counts(value: str) -> dict[str, int]:
    if not isinstance(value, str) or not value:
        return {}
    return {field.split(":")[0]: int(field.split(":")[1]) for field in value.split(";")}


def parse_posterior(value: str) -> dict[str, float]:
    return {field.split(":")[0]: float(field.split(":")[1]) for field in str(value).split(";")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--revision-root", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    revision = args.revision_root.resolve()
    outdir = revision / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    summaries = project / "analysis_v2" / "summaries"
    lineage = pd.read_csv(summaries / "structure_informed_residue_summary_by_lineage.csv")
    counterpart = pd.read_csv(summaries / "acidic_patch_counterpart_residue_summary_by_lineage.csv")
    ancestral = pd.read_csv(revision / "analysis" / "ancestral_state" / "ancestral_state_summary.csv")
    ancestral = ancestral.loc[ancestral["Branch_Length_Scheme"] == "taxonomy_step_lengths"]

    image_path = project / "analysis_v2" / "structure" / "9MPP_DNMT3L_switching_helix.png"
    image = Image.open(image_path).convert("RGB")
    pixels = np.asarray(image)
    nonwhite = np.any(pixels < 245, axis=2)
    yy, xx = np.where(nonwhite)
    if len(xx):
        pad = 45
        image = image.crop(
            (
                max(0, int(xx.min()) - pad),
                max(0, int(yy.min()) - pad),
                min(image.width, int(xx.max()) + pad),
                min(image.height, int(yy.max()) + pad),
            )
        )

    configure()
    fig = plt.figure(figsize=(12, 8.4), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1.02, 0.98])

    ax = fig.add_subplot(grid[0, 0])
    panel_label(ax, "a")
    ax.set_title("Human DNMT3L switching helix in 9MPP", loc="left", fontweight="bold")
    ax.imshow(image)
    ax.axis("off")
    ax.legend(
        handles=[
            Patch(facecolor="#F28E00", label="Switching helix, structural 340-354"),
            Patch(facecolor="#D100C9", label="Q348 and Q351"),
            Patch(facecolor="#00B9BE", label="A357 (NP h358)"),
        ],
        loc="lower left",
        frameon=False,
    )

    ax = fig.add_subplot(grid[0, 1])
    panel_label(ax, "b")
    ax.set_title("Major residue across NP341-360", loc="left", fontweight="bold")
    groups = [group for group in ["Mammalia", "Sauropsida", "Amphibia"] if group in set(lineage["Lineage_Group"])]
    positions = list(range(341, 361))
    for yi, group in enumerate(groups):
        subset = lineage.loc[lineage["Lineage_Group"] == group].set_index("Human_Position")
        for xi, position in enumerate(positions):
            if position in subset.index:
                row = subset.loc[position]
                residue = str(row["Major_Residue"])
                frequency = float(row["Major_Residue_Frequency"])
            else:
                residue, frequency = "-", 0.0
            ax.add_patch(
                Rectangle(
                    (xi, yi), 1, 1,
                    color=RESIDUE_PALETTE.get(residue, "#CCCCCC"),
                    alpha=0.35 + 0.65 * frequency,
                    ec="white",
                    lw=0.8,
                )
            )
            ax.text(xi + 0.5, yi + 0.53, residue, ha="center", va="center", fontsize=7.2, fontweight="bold")
    ax.set_xlim(0, len(positions))
    ax.set_ylim(len(groups), 0)
    ax.set_xticks(np.arange(len(positions)) + 0.5, positions, rotation=90, fontsize=6.8)
    ax.set_yticks(np.arange(len(groups)) + 0.5, [LINEAGE_DISPLAY[group] for group in groups])
    ax.axvline(15, color=INK, lw=1.1)
    ax.text(7.5, -0.23, "switching helix", transform=ax.transData, ha="center", fontsize=7, color=MUTED)
    ax.text(17.5, -0.23, "adjacent loop", transform=ax.transData, ha="center", fontsize=7, color=MUTED)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax = fig.add_subplot(grid[1, 0])
    panel_label(ax, "c")
    ax.set_title("Observed states at acidic-patch counterpart positions", loc="left", fontweight="bold")
    groups = [group for group in ["Mammalia", "Sauropsida", "Amphibia"] if group in set(counterpart["Lineage_Group"])]
    rows = []
    for group in groups:
        for position in (349, 352):
            hit = counterpart.loc[
                (counterpart["Lineage_Group"] == group) & (counterpart["Human_Position"] == position)
            ].iloc[0]
            counts = parse_counts(hit["Amino_Acid_Counts"])
            total = sum(counts.values())
            for residue, count in counts.items():
                rows.append({"Group": group, "Position": position, "Residue": residue, "Fraction": count / total})
    state = pd.DataFrame(rows)
    combinations = [(group, 349) for group in groups] + [(group, 352) for group in groups]
    labels = [f"{LINEAGE_DISPLAY[group]}\nNP{position}" for group, position in combinations]
    x = np.arange(len(combinations))
    bottom = np.zeros(len(combinations))
    residues = sorted(state["Residue"].unique())
    for residue in residues:
        values = np.array(
            [
                state.loc[
                    (state["Group"] == group)
                    & (state["Position"] == position)
                    & (state["Residue"] == residue),
                    "Fraction",
                ].sum()
                for group, position in combinations
            ]
        )
        ax.bar(x, values, bottom=bottom, color=RESIDUE_PALETTE.get(residue, "#BBBBBB"), edgecolor="white", linewidth=0.4, label=residue)
        bottom += values
    ax.axvline(len(groups) - 0.5, color=INK, lw=0.8)
    ax.set_xticks(x, labels, rotation=22, ha="right", fontsize=7)
    ax.set_ylim(0, 1.03)
    ax.set_ylabel("Fraction of resolved sequences")
    ax.legend(title="Residue", frameon=False, ncol=min(6, len(residues)), loc="upper center", bbox_to_anchor=(0.5, -0.27))
    ax.grid(axis="y", color=GRID, linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = fig.add_subplot(grid[1, 1])
    panel_label(ax, "d")
    ax.set_title("Ancestral-state reconstruction on the taxonomy topology", loc="left", fontweight="bold")
    combinations = [(349, "Therian"), (349, "Sauropsid"), (352, "Therian"), (352, "Sauropsid")]
    posterior_rows = []
    for position, group in combinations:
        row = ancestral.loc[ancestral["Human_Position"] == position].iloc[0]
        posterior = parse_posterior(row[f"{group}_MRCA_Posterior"])
        for residue, probability in posterior.items():
            posterior_rows.append({"Position": position, "Group": group, "Residue": residue, "Probability": probability})
    posterior_frame = pd.DataFrame(posterior_rows)
    x = np.arange(len(combinations))
    bottom = np.zeros(len(combinations))
    residues = sorted(posterior_frame["Residue"].unique())
    for residue in residues:
        values = np.array(
            [
                posterior_frame.loc[
                    (posterior_frame["Position"] == position)
                    & (posterior_frame["Group"] == group)
                    & (posterior_frame["Residue"] == residue),
                    "Probability",
                ].sum()
                for position, group in combinations
            ]
        )
        previous = bottom.copy()
        ax.bar(x, values, bottom=bottom, color=RESIDUE_PALETTE.get(residue, "#BBBBBB"), edgecolor="white", linewidth=0.5, label=residue)
        for xi, value, base in zip(x, values, previous):
            if value >= 0.05:
                ax.text(xi, base + value / 2, residue, ha="center", va="center", fontsize=7.2, fontweight="bold", color=INK)
        bottom += values
    labels = [f"{group} MRCA\nNP{position}" for position, group in combinations]
    ax.set_xticks(x, labels, rotation=18, ha="right", fontsize=7.2)
    ax.set_ylim(0, 1.03)
    ax.set_ylabel("Marginal ancestral-state probability")
    ax.text(
        0.0, -0.24,
        "Equal-rates model; branch lengths are compressed NCBI taxonomy steps.\nNP349 is stable under a unit-branch sensitivity; NP352 is not resolved in Sauropsida.",
        transform=ax.transAxes,
        color=MUTED,
        fontsize=7,
    )
    ax.grid(axis="y", color=GRID, linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for suffix in ("png", "pdf"):
        path = outdir / f"Figure_3_structural_state_evolution_revised.{suffix}"
        fig.savefig(path, dpi=350 if suffix == "png" else None, bbox_inches="tight", pad_inches=0.08)
        print(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
