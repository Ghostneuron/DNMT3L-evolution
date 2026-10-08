from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
from PIL import Image


PROJECT = Path(__file__).resolve().parents[1]
V2 = PROJECT / "analysis_v2"
SUMMARIES = V2 / "summaries"
OUT = PROJECT / "figures_v2"
STRUCTURE_IMAGE = V2 / "structure" / "9MPP_DNMT3L_switching_helix.png"
SOURCE_ORTHOLOGY = (
    PROJECT
    / "source_package"
    / "analysis"
    / "orthology_rebuild_20260826"
    / "orthology_audit"
    / "DNMT3L_orthology_audit_summary.csv"
)
ORTHOLOGY_SENSITIVITY = V2 / "orthology_sensitivity" / "group_separation_summary.csv"

INK = "#20262E"
MUTED = "#687683"
GRID = "#D8DEE4"
BLUE = "#3976A8"
TEAL = "#168C88"
CORAL = "#D2644F"
GOLD = "#D5A027"
GREEN = "#4D8D5A"
PURPLE = "#765AA6"
LIGHT = "#EEF2F4"
SCOPE_ORDER = ["Placentalia", "Euarchontoglires", "Laurasiatheria", "Sauropsida"]
ALIGNER_ORDER = ["MACSE", "MAFFT", "PRANK"]
SCOPE_COLORS = {
    "Placentalia": BLUE,
    "Euarchontoglires": TEAL,
    "Laurasiatheria": CORAL,
    "Sauropsida": GOLD,
}
LINEAGE_DISPLAY = {"Mammalia": "Sampled Theria", "Sauropsida": "Sauropsida", "Amphibia": "Amphibia"}


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def configure() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9,
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "legend.fontsize": 7.5,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def style_axes(ax, *, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.6, alpha=0.8)
        ax.set_axisbelow(True)


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.09,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
    )


def save(fig, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        path = OUT / f"{stem}.{suffix}"
        fig.savefig(path, dpi=350 if suffix == "png" else None, bbox_inches="tight", pad_inches=0.08)
        print(path.relative_to(PROJECT))
    plt.close(fig)


def figure_1_data_foundation() -> None:
    filtered = pd.read_csv(SUMMARIES / "filtered_alignment_summary.csv")
    concordance = pd.read_csv(SUMMARIES / "cross_aligner_concordance_by_site.csv")
    orth = pd.read_csv(SOURCE_ORTHOLOGY).iloc[0]
    family = pd.read_csv(ORTHOLOGY_SENSITIVITY)

    fig = plt.figure(figsize=(12, 8.1), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[0.88, 1.12])

    ax = fig.add_subplot(grid[0, 0])
    panel_label(ax, "a")
    ax.set_title("Record-level quality filters", loc="left", fontweight="bold")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    ax.axis("off")
    both_qc = int(orth["Selection_Ready"] + orth["Hybrid_Taxa"])
    stages = [
        ("Archived\ncandidates", int(orth["Candidate_Records"]), BLUE, 0.15, 1.24),
        ("Sequence QC", int(orth["Sequence_QC_Pass"]), TEAL, 2.55, 2.08),
        ("Coverage +\nh358 mapping", int(orth["Domain_QC_Pass"]), GOLD, 2.55, 0.40),
        ("Both QC\ngates", both_qc, MUTED, 5.30, 1.24),
        ("Selection-\nready", int(orth["Selection_Ready"]), CORAL, 8.05, 1.24),
    ]
    for label, value, color, x, y in stages:
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                1.72,
                1.30,
                boxstyle="round,pad=0.04,rounding_size=0.07",
                facecolor=color,
                edgecolor="none",
            )
        )
        ax.text(x + 0.86, y + 0.86, str(value), color="white", fontsize=16, fontweight="bold", ha="center")
        ax.text(x + 0.86, y + 0.34, label, color="white", fontsize=7.2, ha="center", va="center")
    arrow = {"arrowstyle": "-|>", "color": MUTED, "lw": 1.25}
    ax.annotate("", xy=(2.47, 2.70), xytext=(1.96, 2.12), arrowprops=arrow)
    ax.annotate("", xy=(2.47, 1.05), xytext=(1.96, 1.65), arrowprops=arrow)
    ax.annotate("", xy=(5.22, 2.10), xytext=(4.35, 2.55), arrowprops=arrow)
    ax.annotate("", xy=(5.22, 1.67), xytext=(4.35, 1.05), arrowprops=arrow)
    ax.annotate("", xy=(7.97, 1.89), xytext=(7.10, 1.89), arrowprops=arrow)
    ax.text(7.54, 2.72, "remove hybrid (n=1)", color=MUTED, fontsize=6.2, ha="center")
    ax.text(
        0.15,
        0.12,
        "Sequence and coverage gates were parallel filters; selection-ready records passed both and were nonhybrid.",
        color=MUTED,
        fontsize=7.1,
    )

    ax = fig.add_subplot(grid[0, 1])
    panel_label(ax, "b")
    ax.set_title("DNMT3-family placement including DNMT3C", loc="left", fontweight="bold")
    labels = {
        "DNMT3L": "DNMT3L",
        "DNMT3A": "DNMT3A",
        "DNMT3B_plus_DNMT3C": "DNMT3B + DNMT3C",
        "DNMT3A_plus_DNMT3B_plus_DNMT3C": "All catalytic DNMT3",
    }
    frame = family.copy()
    frame["Label"] = frame["Group"].map(labels)
    frame = frame.dropna(subset=["Label"])
    y = np.arange(len(frame))[::-1]
    colors = [TEAL, BLUE, CORAL, GOLD][: len(frame)]
    ax.barh(y, frame["SH_aLRT_Support"], color=colors, height=0.55)
    ax.set_yticks(y, [f"{label} (n={int(n)})" for label, n in zip(frame["Label"], frame["N_Tips"])])
    ax.set_xlim(0, 104)
    ax.set_xlabel("SH-aLRT support")
    for yi, value in zip(y, frame["SH_aLRT_Support"]):
        ax.text(value - 1.1, yi, f"{value:.1f}", va="center", ha="right", color="white", fontweight="bold")
    ax.text(0.0, -0.24, "No DNMT3L candidate was nearest to DNMT3A, DNMT3B, or DNMT3C.", transform=ax.transAxes, color=MUTED)
    style_axes(ax, grid_axis="x")

    ax = fig.add_subplot(grid[1, 0])
    panel_label(ax, "c")
    ax.set_title("Each aligner was filtered independently", loc="left", fontweight="bold")
    x = np.arange(len(SCOPE_ORDER))
    width = 0.23
    for offset, aligner, color in zip((-width, 0, width), ALIGNER_ORDER, (BLUE, TEAL, CORAL)):
        values = [
            int(filtered.loc[(filtered["Scope"] == scope) & (filtered["Aligner"] == aligner), "Retained_Codons"].iloc[0])
            for scope in SCOPE_ORDER
        ]
        ax.bar(x + offset, values, width, color=color, label=aligner)
        for xi, value in zip(x + offset, values):
            ax.text(xi, value + 4, str(value), ha="center", va="bottom", fontsize=6.6, rotation=90)
    ax.set_xticks(x, SCOPE_ORDER, rotation=15, ha="right")
    ax.set_ylabel("Retained human-mapped codons")
    ax.set_ylim(0, 420)
    ax.legend(frameon=False, ncol=3, loc="lower left")
    style_axes(ax, grid_axis="y")

    ax = fig.add_subplot(grid[1, 1])
    panel_label(ax, "d")
    ax.set_title("Cross-aligner amino-acid concordance", loc="left", fontweight="bold")
    rows = []
    for scope in SCOPE_ORDER:
        subset = concordance.loc[concordance["Scope"] == scope]
        eligible = subset["All_Three_Resolved_Fraction"] >= 0.80
        rows.append(
            {
                "Scope": scope,
                "Resolved": int(eligible.sum()),
                "Stable": int((eligible & (subset["All_Three_Exact_Among_Resolved"] >= 0.90)).sum()),
            }
        )
    stable = pd.DataFrame(rows)
    y = np.arange(len(stable))[::-1]
    ax.barh(y, stable["Resolved"], color=LIGHT, edgecolor="#AEB9C2", label=">=80% jointly resolved")
    ax.barh(y, stable["Stable"], color=[SCOPE_COLORS[s] for s in stable["Scope"]], label="plus >=90% exact")
    ax.set_yticks(y, stable["Scope"])
    ax.set_xlabel("Human reference positions")
    for yi, kept, total in zip(y, stable["Stable"], stable["Resolved"]):
        ax.text(total + 4, yi, f"{kept}/{total}", va="center", fontsize=7.5)
    ax.set_xlim(0, max(stable["Resolved"]) * 1.15)
    ax.legend(frameon=False, loc="lower right")
    style_axes(ax, grid_axis="x")

    save(fig, "Figure_1_data_foundation")


def figure_2_constraint_architecture() -> None:
    sites = pd.read_csv(SUMMARIES / "amino_acid_constraint_by_site.csv")
    domains = pd.read_csv(SUMMARIES / "amino_acid_constraint_by_domain.csv")
    helix = pd.read_csv(SUMMARIES / "switching_helix_constraint_tests.csv")
    window = pd.read_csv(SUMMARIES / "switching_helix_window_sensitivity.csv")
    tree_window = pd.read_csv(SUMMARIES / "tree_aware_parsimony_window_sensitivity.csv")

    fig = plt.figure(figsize=(12, 8.3), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1.08, 0.92])

    ax = fig.add_subplot(grid[0, :])
    panel_label(ax, "a")
    ax.set_title("Lineage-specific amino-acid variability along DNMT3L", loc="left", fontweight="bold")
    for scope in SCOPE_ORDER:
        subset = sites.loc[sites["Scope"] == scope].sort_values("Human_Position")
        rolling = subset["Median_AA_Entropy"].rolling(11, center=True, min_periods=4).median()
        ax.plot(subset["Human_Position"], rolling, color=SCOPE_COLORS[scope], linewidth=1.7, label=scope)
    ax.axvspan(41, 173, color=TEAL, alpha=0.08)
    ax.axvspan(341, 355, color=GOLD, alpha=0.15)
    ax.axvspan(356, 360, color=CORAL, alpha=0.11)
    ax.text(107, 2.75, "ADD", ha="center", color=TEAL, fontweight="bold")
    ax.text(348, 2.75, "switching helix", ha="center", color="#8B6713", fontweight="bold")
    ax.set_xlim(1, 387)
    ax.set_ylim(0, 3.0)
    ax.set_xlabel("Human NP_037501.2 position")
    ax.set_ylabel("11-site median Shannon entropy")
    ax.legend(frameon=False, ncol=4, loc="upper center")
    style_axes(ax, grid_axis="y")

    ax = fig.add_subplot(grid[1, 0])
    panel_label(ax, "b")
    ax.set_title("Domain-level constraint", loc="left", fontweight="bold")
    wanted = domains.loc[domains["Domain"].isin(["ADD", "C_terminal_MTase_like", "N_terminal"])].copy()
    domain_order = ["ADD", "C_terminal_MTase_like", "N_terminal"]
    display = ["ADD", "C-terminal", "N-terminal"]
    x = np.arange(len(SCOPE_ORDER))
    width = 0.23
    colors = [TEAL, BLUE, CORAL]
    for index, (domain, label, color) in enumerate(zip(domain_order, display, colors)):
        values = []
        for scope in SCOPE_ORDER:
            hit = wanted.loc[(wanted["Scope"] == scope) & (wanted["Domain"] == domain), "Median_AA_Entropy"]
            values.append(float(hit.iloc[0]) if len(hit) else np.nan)
        ax.bar(x + (index - 1) * width, values, width, color=color, label=label)
    ax.set_xticks(x, SCOPE_ORDER, rotation=15, ha="right")
    ax.set_ylabel("Median site entropy")
    ax.legend(frameon=False, ncol=3)
    style_axes(ax, grid_axis="y")

    ax = fig.add_subplot(grid[1, 1])
    panel_label(ax, "c")
    ax.set_title("Structural region versus C-terminal background", loc="left", fontweight="bold")
    region_order = ["switching_helix", "switching_helix_adjacent_loop"]
    labels = ["Switching helix", "Adjacent loop"]
    y = np.arange(len(SCOPE_ORDER))
    for offset, region, label, color in zip((-0.12, 0.12), region_order, labels, (GOLD, CORAL)):
        values = []
        qvals = []
        for scope in SCOPE_ORDER:
            row = helix.loc[(helix["Scope"] == scope) & (helix["Structural_Region"] == region)].iloc[0]
            values.append(row["Region_Median_Entropy"] / row["Background_Median_Entropy"])
            qvals.append(row["q_within_scope"])
        ax.scatter(values, y + offset, s=55, color=color, edgecolor="white", linewidth=0.8, label=label, zorder=3)
        for value, yi, q in zip(values, y + offset, qvals):
            if q <= 0.05:
                ax.text(value + 0.08, yi, f"q={q:.3g}", va="center", fontsize=6.6, color=MUTED)
    ax.axvline(1.0, color=INK, linewidth=1, linestyle="--")
    ax.set_yticks(y, SCOPE_ORDER)
    ax.set_xlabel("Median entropy ratio (region / background)")
    ax.set_xlim(0, max(4.35, ax.get_xlim()[1]))
    ax.legend(frameon=False, loc="upper left")
    evaluable_window_q = window["BH_q_across_evaluable_tests"].dropna()
    evaluable_tree_window_q = tree_window["BH_q_across_evaluable_tests"].dropna()
    if not evaluable_window_q.empty and not evaluable_tree_window_q.empty:
        ax.text(
            0.0,
            -0.22,
            "Same-length-window sensitivity: no corrected support "
            f"(entropy q>={evaluable_window_q.min():.3f}; tree-aware q>={evaluable_tree_window_q.min():.3f}).",
            transform=ax.transAxes,
            color=MUTED,
            fontsize=7,
        )
    style_axes(ax, grid_axis="x")

    save(fig, "Figure_2_constraint_architecture")


def _parse_counts(value: str) -> dict[str, int]:
    result: dict[str, int] = {}
    if not isinstance(value, str) or not value:
        return result
    for field in value.split(";"):
        residue, count = field.split(":")
        result[residue] = int(count)
    return result


def figure_3_structural_states() -> None:
    lineage = pd.read_csv(SUMMARIES / "structure_informed_residue_summary_by_lineage.csv")
    counterpart = pd.read_csv(SUMMARIES / "acidic_patch_counterpart_residue_summary_by_lineage.csv")
    image = Image.open(STRUCTURE_IMAGE).convert("RGB")
    pixels = np.asarray(image)
    nonwhite = np.any(pixels < 245, axis=2)
    yy, xx = np.where(nonwhite)
    if len(xx):
        pad = 45
        left = max(0, int(xx.min()) - pad)
        right = min(image.width, int(xx.max()) + pad)
        top = max(0, int(yy.min()) - pad)
        bottom = min(image.height, int(yy.max()) + pad)
        image = image.crop((left, top, right, bottom))

    fig = plt.figure(figsize=(12, 8.2), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=[1.05, 0.95], height_ratios=[1.1, 0.9])

    ax = fig.add_subplot(grid[:, 0])
    panel_label(ax, "a")
    ax.set_title("DNMT3L C-terminal switching helix in 9MPP", loc="left", fontweight="bold")
    ax.imshow(image)
    ax.axis("off")
    ax.legend(
        handles=[
            Patch(facecolor="#F28E00", label="Switching helix, structural 340-354"),
            Patch(facecolor="#D100C9", label="Q348 and Q351"),
            Patch(facecolor="#00B9BE", label="A357 (NP h358)"),
            Patch(facecolor="#BEBEBE", label="Remaining DNMT3L C-terminal domain"),
        ],
        loc="lower left",
        frameon=False,
        fontsize=7.5,
    )

    ax = fig.add_subplot(grid[0, 1])
    panel_label(ax, "b")
    ax.set_title("Major residue across structural positions 341-360", loc="left", fontweight="bold")
    groups = [group for group in ["Mammalia", "Sauropsida", "Amphibia"] if group in set(lineage["Lineage_Group"])]
    positions = list(range(341, 361))
    residue_palette = {
        "A": "#8EC07C", "C": "#D6B656", "D": "#D65F5F", "E": "#C94B4B", "F": "#8F79B8",
        "G": "#A8B0B8", "H": "#6F87B8", "I": "#4F9D8D", "K": "#3C78A8", "L": "#2E8B78",
        "M": "#448F83", "N": "#D89C48", "P": "#B8739B", "Q": "#E2A83B", "R": "#245E91",
        "S": "#73A9C2", "T": "#7AB7B2", "V": "#60A37B", "W": "#6F5B9B", "Y": "#7B67A7", "-": "#E4E8EB",
    }
    for yi, group in enumerate(groups):
        subset = lineage.loc[lineage["Lineage_Group"] == group].set_index("Human_Position")
        for xi, position in enumerate(positions):
            if position not in subset.index:
                residue, frequency = "-", 0.0
            else:
                row = subset.loc[position]
                residue, frequency = str(row["Major_Residue"]), float(row["Major_Residue_Frequency"])
            ax.add_patch(Rectangle((xi, yi), 1, 1, color=residue_palette.get(residue, "#CCCCCC"), alpha=0.35 + 0.65 * frequency, ec="white", lw=0.8))
            ax.text(xi + 0.5, yi + 0.53, residue, ha="center", va="center", fontsize=7.5, fontweight="bold")
    ax.set_xlim(0, len(positions))
    ax.set_ylim(len(groups), 0)
    ax.set_xticks(np.arange(len(positions)) + 0.5, positions, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(len(groups)) + 0.5, [LINEAGE_DISPLAY[group] for group in groups])
    ax.axvline(15, color=INK, lw=1.2)
    ax.text(7.5, -0.23, "switching helix", transform=ax.transData, ha="center", fontsize=7, color=MUTED)
    ax.text(17.5, -0.23, "adjacent loop", transform=ax.transData, ha="center", fontsize=7, color=MUTED)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax = fig.add_subplot(grid[1, 1])
    panel_label(ax, "c")
    ax.set_title("Lineage states at acidic-patch counterpart positions", loc="left", fontweight="bold")
    groups = [group for group in ["Mammalia", "Sauropsida", "Amphibia"] if group in set(counterpart["Lineage_Group"])]
    rows = []
    for group in groups:
        for position in (349, 352):
            hit = counterpart.loc[(counterpart["Lineage_Group"] == group) & (counterpart["Human_Position"] == position)].iloc[0]
            counts = _parse_counts(hit["Amino_Acid_Counts"])
            total = sum(counts.values())
            for residue, count in counts.items():
                rows.append({"Group": group, "Position": position, "Residue": residue, "Fraction": count / total})
    state = pd.DataFrame(rows)
    labels = [f"{LINEAGE_DISPLAY[group]}\nNP349 / Q348" for group in groups] + [f"{LINEAGE_DISPLAY[group]}\nNP352 / Q351" for group in groups]
    combinations = [(group, 349) for group in groups] + [(group, 352) for group in groups]
    x = np.arange(len(combinations))
    bottom = np.zeros(len(combinations))
    residues = sorted(state["Residue"].unique())
    for residue in residues:
        values = np.array([
            state.loc[(state["Group"] == group) & (state["Position"] == position) & (state["Residue"] == residue), "Fraction"].sum()
            for group, position in combinations
        ])
        ax.bar(x, values, bottom=bottom, color=residue_palette.get(residue, "#BBBBBB"), edgecolor="white", linewidth=0.4, label=residue)
        bottom += values
    ax.axvline(len(groups) - 0.5, color=INK, lw=0.8)
    ax.set_xticks(x, labels, rotation=25, ha="right", fontsize=7)
    ax.set_ylim(0, 1.03)
    ax.set_ylabel("Fraction of resolved sequences")
    ax.legend(title="Residue", frameon=False, ncol=min(6, len(residues)), loc="upper center", bbox_to_anchor=(0.5, -0.28))
    style_axes(ax, grid_axis="y")

    save(fig, "Figure_3_structural_state_evolution")


def figure_4_selection_robustness() -> None:
    robustness = pd.read_csv(SUMMARIES / "site_robustness_by_scope.csv")
    counts = pd.read_csv(SUMMARIES / "selection_call_counts_by_alignment.csv")
    busted = pd.read_csv(SUMMARIES / "BUSTED_gene_wide_summary.csv")
    h358 = pd.read_csv(SUMMARIES / "h358_expanded_audit.csv")
    topology_path = SUMMARIES / "h358_tree_topology_sensitivity.csv"
    topology = pd.read_csv(topology_path) if topology_path.exists() else pd.DataFrame()

    fig = plt.figure(figsize=(12, 8.4), layout="constrained")
    grid = fig.add_gridspec(2, 2)

    ax = fig.add_subplot(grid[0, 0])
    panel_label(ax, "a")
    ax.set_title("Robust site classes across three independent aligners", loc="left", fontweight="bold")
    categories = ["Robust_Pervasive_q10", "Robust_Episodic_q10", "Robust_Purifying_q05"]
    labels = ["Pervasive\n(FUBAR + FEL)", "Episodic\n(MEME)", "Purifying\n(FUBAR + FEL)"]
    x = np.arange(len(SCOPE_ORDER))
    width = 0.23
    colors = [CORAL, PURPLE, TEAL]
    for offset, category, label, color in zip((-width, 0, width), categories, labels, colors):
        values = [int(as_bool(robustness.loc[robustness["Scope"] == scope, category]).sum()) for scope in SCOPE_ORDER]
        ax.bar(x + offset, values, width, color=color, label=label)
        for xi, value in zip(x + offset, values):
            if value:
                ax.text(xi, value + 0.4, str(value), ha="center", fontsize=7)
    ax.set_xticks(x, SCOPE_ORDER, rotation=15, ha="right")
    ax.set_ylabel("Robust human-mapped sites")
    ax.set_ylim(0, 270)
    ax.legend(frameon=False, ncol=3, loc="upper center")
    style_axes(ax, grid_axis="y")

    ax = fig.add_subplot(grid[0, 1])
    panel_label(ax, "b")
    ax.set_title("Gene-wide BUSTED evidence", loc="left", fontweight="bold")
    busted = busted.copy()
    busted_cap = 16.0
    busted["minus_log10_q"] = np.minimum(
        -np.log10(busted["BUSTED_q_across_12"].clip(lower=10 ** -busted_cap)),
        busted_cap,
    )
    x = np.arange(len(SCOPE_ORDER))
    width = 0.23
    for offset, aligner, color in zip((-width, 0, width), ALIGNER_ORDER, (BLUE, TEAL, CORAL)):
        values = []
        for scope in SCOPE_ORDER:
            hit = busted.loc[(busted["Scope"] == scope) & (busted["Aligner"] == aligner), "minus_log10_q"]
            values.append(float(hit.iloc[0]) if len(hit) else np.nan)
        ax.bar(x + offset, values, width, color=color, label=aligner)
    ax.axhline(-math.log10(0.05), color=INK, linestyle="--", linewidth=1, label="q = 0.05")
    ax.set_xticks(x, SCOPE_ORDER, rotation=15, ha="right")
    ax.set_ylabel("-log10(BH q across 12 tests), capped at 16")
    ax.set_ylim(0, busted_cap + 2.0)
    ax.legend(frameon=False, ncol=2)
    style_axes(ax, grid_axis="y")

    ax = fig.add_subplot(grid[1, 0])
    panel_label(ax, "c")
    ax.set_title("h358 evidence across scopes and aligners", loc="left", fontweight="bold")
    display = h358.copy()
    display["Scope"] = pd.Categorical(display["Scope"], categories=SCOPE_ORDER, ordered=True)
    display["Aligner"] = pd.Categorical(display["Aligner"], categories=ALIGNER_ORDER, ordered=True)
    display = display.sort_values(["Scope", "Aligner"]).reset_index(drop=True)
    display["Label"] = (
        display["Scope"].astype(str).str.replace("Euarchontoglires", "Euarch.", regex=False)
        + "\n"
        + display["Aligner"].astype(str)
    )
    metric_values = display[["FUBAR_PP_positive", "FEL_q", "MEME_q"]].to_numpy(dtype=float)
    threshold_pass = np.column_stack(
        [
            display["FUBAR_PP_positive"].to_numpy(dtype=float) >= 0.90,
            display["FEL_q"].to_numpy(dtype=float) <= 0.10,
            display["MEME_q"].to_numpy(dtype=float) <= 0.10,
        ]
    ).astype(float)
    evidence_cmap = LinearSegmentedColormap.from_list("threshold_pass", ["#F1F3F4", "#2A9D8F"])
    ax.imshow(threshold_pass, aspect="auto", cmap=evidence_cmap, vmin=0, vmax=1)
    ax.set_xticks(range(3), ["FUBAR PP", "FEL q", "MEME q"])
    ax.set_yticks(range(len(display)), display["Label"], fontsize=6.5)
    for yi in range(metric_values.shape[0]):
        for xi in range(metric_values.shape[1]):
            value = metric_values[yi, xi]
            label = f"{value:.3f}" if value >= 0.001 else f"{value:.1e}"
            ax.text(
                xi,
                yi,
                label,
                ha="center",
                va="center",
                fontsize=6.2,
                color="white" if threshold_pass[yi, xi] else INK,
            )
    ax.text(
        0.0,
        -0.13,
        "Shading marks PP >= 0.90 or q <= 0.10; printed values retain each model's native scale.",
        transform=ax.transAxes,
        color=MUTED,
        fontsize=7,
    )
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax = fig.add_subplot(grid[1, 1])
    panel_label(ax, "d")
    ax.set_title("Alignment and tree sensitivity", loc="left", fontweight="bold")
    criteria = ["Joint_Pervasive_q10", "MEME_q10", "Joint_Purifying_q05"]
    names = ["Pervasive", "Episodic", "Purifying"]
    matrix = []
    row_labels = []
    scope_labels = {
        "Placentalia": "Placentalia",
        "Euarchontoglires": "Euarch.",
        "Laurasiatheria": "Laurasiatheria",
        "Sauropsida": "Sauropsida",
    }
    for scope in SCOPE_ORDER:
        for aligner in ALIGNER_ORDER:
            hit = counts.loc[(counts["Scope"] == scope) & (counts["Aligner"] == aligner)]
            if hit.empty:
                continue
            matrix.append([int(hit.iloc[0][criterion]) for criterion in criteria])
            row_labels.append(f"{scope_labels[scope]}\n{aligner}")
    array = np.asarray(matrix, dtype=float)
    cmap = LinearSegmentedColormap.from_list("calls", ["#F4F6F7", GOLD, CORAL])
    image_plot = ax.imshow(array, aspect="auto", cmap=cmap)
    ax.set_xticks(range(len(names)), names)
    ax.set_yticks(range(len(row_labels)), row_labels, fontsize=6.2)
    for yi in range(array.shape[0]):
        for xi in range(array.shape[1]):
            ax.text(xi, yi, str(int(array[yi, xi])), ha="center", va="center", fontsize=7)
    if not topology.empty:
        changed = ((topology["FUBAR_PP_positive"] >= 0.90) != (topology["SpeciesTree_FUBAR_PP_positive"] >= 0.90)).sum()
        ax.text(0.0, -0.16, f"At h358, {int(changed)}/4 MACSE scopes cross the FUBAR threshold when tree topology changes.", transform=ax.transAxes, color=MUTED, fontsize=7)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.colorbar(image_plot, ax=ax, shrink=0.7, label="Sites called in one alignment")

    save(fig, "Figure_4_selection_robustness")


def figure_s1_saturation_diagnostics() -> None:
    summary = pd.read_csv(V2 / "saturation_diagnostics" / "codon_saturation_summary.csv")
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.25), layout="constrained")
    x = np.arange(len(SCOPE_ORDER))
    offsets = {"MACSE": -0.22, "MAFFT": 0.0, "PRANK": 0.22}
    aligner_colors = {"MACSE": BLUE, "MAFFT": TEAL, "PRANK": CORAL}

    ax = axes[0]
    panel_label(ax, "a")
    ax.set_title("Pairwise synonymous divergence", loc="left", fontweight="bold")
    for aligner in ALIGNER_ORDER:
        rows = summary.set_index(["Scope", "Aligner"])
        medians = np.array([rows.loc[(scope, aligner), "NG86_dS_Median"] for scope in SCOPE_ORDER], dtype=float)
        upper = np.array([rows.loc[(scope, aligner), "NG86_dS_P90"] for scope in SCOPE_ORDER], dtype=float)
        ax.errorbar(
            x + offsets[aligner],
            medians,
            yerr=np.vstack([np.zeros_like(medians), np.maximum(upper - medians, 0)]),
            fmt="o",
            color=aligner_colors[aligner],
            capsize=3,
            linewidth=1.2,
            markersize=5,
            label=aligner,
        )
    ax.axhline(1, color=INK, linestyle="--", linewidth=0.9)
    ax.axhline(2, color=INK, linestyle=":", linewidth=0.9)
    ax.set_xticks(x, SCOPE_ORDER, rotation=18, ha="right")
    ax.set_ylabel("NG86 dS (point: median; whisker: P90)")
    ax.legend(frameon=False, ncol=3, loc="upper left")
    style_axes(ax, grid_axis="y")

    ax = axes[1]
    panel_label(ax, "b")
    ax.set_title("Upper-tail multiple-hit indicators", loc="left", fontweight="bold")
    width = 0.20
    for aligner in ALIGNER_ORDER:
        rows = summary.set_index(["Scope", "Aligner"])
        ge1 = np.array([rows.loc[(scope, aligner), "NG86_dS_ge_1_Pct"] for scope in SCOPE_ORDER], dtype=float)
        ge2 = np.array([rows.loc[(scope, aligner), "NG86_dS_ge_2_Pct"] for scope in SCOPE_ORDER], dtype=float)
        xpos = x + offsets[aligner]
        ax.bar(xpos, ge1, width, color=aligner_colors[aligner], alpha=0.42)
        ax.bar(xpos, ge2, width, color=aligner_colors[aligner], alpha=0.95)
    ax.set_xticks(x, SCOPE_ORDER, rotation=18, ha="right")
    ax.set_ylabel("Finite pairwise estimates (%)")
    ax.legend(
        handles=[
            Patch(facecolor=INK, alpha=0.42, label="dS >= 1"),
            Patch(facecolor=INK, alpha=0.95, label="dS >= 2"),
        ],
        frameon=False,
        loc="upper left",
    )
    style_axes(ax, grid_axis="y")

    ax = axes[2]
    panel_label(ax, "c")
    ax.set_title("Third-codon-position divergence", loc="left", fontweight="bold")
    for aligner in ALIGNER_ORDER:
        rows = summary.set_index(["Scope", "Aligner"])
        medians = np.array([rows.loc[(scope, aligner), "Pos3_PDistance_Median"] for scope in SCOPE_ORDER], dtype=float)
        upper = np.array([rows.loc[(scope, aligner), "Pos3_PDistance_P90"] for scope in SCOPE_ORDER], dtype=float)
        ax.errorbar(
            x + offsets[aligner],
            medians,
            yerr=np.vstack([np.zeros_like(medians), np.maximum(upper - medians, 0)]),
            fmt="o",
            color=aligner_colors[aligner],
            capsize=3,
            linewidth=1.2,
            markersize=5,
        )
        for xi, scope in zip(x + offsets[aligner], SCOPE_ORDER):
            undefined = float(rows.loc[(scope, aligner), "Pos3_K80_Undefined_Pct"])
            if undefined > 0:
                ax.text(xi, upper[list(SCOPE_ORDER).index(scope)] + 0.015, f"{undefined:.1f}% undef.", ha="center", fontsize=6, color=MUTED)
    ax.set_xticks(x, SCOPE_ORDER, rotation=18, ha="right")
    ax.set_ylabel("P distance (point: median; whisker: P90)")
    style_axes(ax, grid_axis="y")

    save(fig, "Figure_S1_saturation_diagnostics")


def figure_s2_tree_aware_constraint() -> None:
    entropy_domains = pd.read_csv(SUMMARIES / "amino_acid_constraint_by_domain.csv")
    tree_domains = pd.read_csv(SUMMARIES / "tree_aware_parsimony_domain_tests.csv")
    tree_structural = pd.read_csv(SUMMARIES / "tree_aware_parsimony_structural_tests.csv")
    entropy_windows = pd.read_csv(SUMMARIES / "switching_helix_window_sensitivity.csv")
    tree_windows = pd.read_csv(SUMMARIES / "tree_aware_parsimony_window_sensitivity.csv")

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.6), layout="constrained")
    x = np.arange(len(SCOPE_ORDER))

    ax = axes[0]
    panel_label(ax, "a")
    ax.set_title("ADD versus C-terminal constraint", loc="left", fontweight="bold")
    entropy_ratios = []
    tree_ratios = []
    for scope in SCOPE_ORDER:
        entropy = entropy_domains.loc[entropy_domains["Scope"] == scope].set_index("Domain")
        entropy_ratios.append(
            entropy.loc["ADD", "Median_AA_Entropy"]
            / entropy.loc["C_terminal_MTase_like", "Median_AA_Entropy"]
        )
        tree = tree_domains.loc[tree_domains["Scope"] == scope].iloc[0]
        tree_ratios.append(
            tree.ADD_Median_Changes_Per_100 / tree.C_terminal_Median_Changes_Per_100
        )
    width = 0.32
    ax.bar(x - width / 2, entropy_ratios, width, color=BLUE, label="Entropy ratio")
    ax.bar(x + width / 2, tree_ratios, width, color=TEAL, label="Tree-aware change ratio")
    ax.axhline(1.0, color=INK, linestyle="--", linewidth=0.9)
    ax.set_xticks(x, SCOPE_ORDER, rotation=18, ha="right")
    ax.set_ylabel("ADD / C-terminal median")
    ax.legend(frameon=False, fontsize=7)
    style_axes(ax, grid_axis="y")

    ax = axes[1]
    panel_label(ax, "b")
    ax.set_title("Tree-aware structural-region contrast", loc="left", fontweight="bold")
    y = np.arange(len(SCOPE_ORDER))
    for offset, region, label, color in zip(
        (-0.12, 0.12),
        ("switching_helix", "switching_helix_adjacent_loop"),
        ("Switching helix", "Adjacent loop"),
        (GOLD, CORAL),
    ):
        values = []
        for scope in SCOPE_ORDER:
            row = tree_structural.loc[
                (tree_structural["Scope"] == scope)
                & (tree_structural["Structural_Region"] == region)
            ].iloc[0]
            values.append(
                row.Region_Median_Changes_Per_100 / row.Background_Median_Changes_Per_100
            )
        ax.scatter(values, y + offset, s=55, color=color, edgecolor="white", linewidth=0.8, label=label)
    ax.axvline(1.0, color=INK, linestyle="--", linewidth=0.9)
    ax.set_yticks(y, SCOPE_ORDER)
    ax.set_xlabel("Median minimum-change ratio")
    ax.legend(frameon=False, fontsize=7)
    style_axes(ax, grid_axis="x")

    ax = axes[2]
    panel_label(ax, "c")
    ax.set_title("Equal-length-window sensitivity", loc="left", fontweight="bold")
    entropy = entropy_windows.merge(
        tree_windows,
        on=["Scope", "Structural_Region"],
        suffixes=("_entropy", "_tree"),
    )
    labels = []
    entropy_q = []
    tree_q = []
    scope_short = {
        "Placentalia": "Plac.",
        "Euarchontoglires": "Euarch.",
        "Laurasiatheria": "Laurasia.",
        "Sauropsida": "Saurop.",
    }
    for scope in SCOPE_ORDER:
        for region, region_label in (
            ("switching_helix", "helix"),
            ("switching_helix_adjacent_loop", "loop"),
        ):
            row = entropy.loc[
                (entropy["Scope"] == scope) & (entropy["Structural_Region"] == region)
            ].iloc[0]
            labels.append(f"{scope_short[scope]} {region_label}")
            entropy_q.append(row.BH_q_across_evaluable_tests_entropy)
            tree_q.append(row.BH_q_across_evaluable_tests_tree)
    y = np.arange(len(labels))[::-1]
    entropy_values = -np.log10(pd.Series(entropy_q, dtype=float))
    tree_values = -np.log10(pd.Series(tree_q, dtype=float))
    ax.scatter(entropy_values, y - 0.12, color=BLUE, s=36, label="Entropy")
    ax.scatter(tree_values, y + 0.12, color=TEAL, marker="s", s=34, label="Tree-aware")
    for yi, entropy_value, tree_value in zip(y, entropy_values, tree_values):
        if np.isnan(entropy_value) and np.isnan(tree_value):
            ax.text(0.55, yi, "NA", ha="center", va="center", color=MUTED, fontsize=7)
    ax.axvline(-math.log10(0.05), color=INK, linestyle="--", linewidth=0.9, label="q = 0.05")
    ax.set_yticks(y, labels, fontsize=7)
    ax.set_xlabel("-log10(corrected q)")
    ax.legend(frameon=False, fontsize=7)
    style_axes(ax, grid_axis="x")

    save(fig, "Figure_S2_tree_aware_constraint")


def main() -> None:
    configure()
    figure_1_data_foundation()
    figure_2_constraint_architecture()
    figure_3_structural_states()
    required = [
        SUMMARIES / "site_robustness_by_scope.csv",
        SUMMARIES / "selection_call_counts_by_alignment.csv",
        SUMMARIES / "BUSTED_gene_wide_summary.csv",
        SUMMARIES / "h358_expanded_audit.csv",
    ]
    if all(path.exists() for path in required):
        figure_4_selection_robustness()
    else:
        print("Selection summaries are not complete; Figure 4 was deferred.")
    saturation_path = V2 / "saturation_diagnostics" / "codon_saturation_summary.csv"
    if saturation_path.exists():
        figure_s1_saturation_diagnostics()
    else:
        print("Saturation summaries are not complete; Figure S1 was deferred.")
    tree_required = [
        SUMMARIES / "tree_aware_parsimony_domain_tests.csv",
        SUMMARIES / "tree_aware_parsimony_structural_tests.csv",
        SUMMARIES / "tree_aware_parsimony_window_sensitivity.csv",
        SUMMARIES / "switching_helix_window_sensitivity.csv",
    ]
    if all(path.exists() for path in tree_required):
        figure_s2_tree_aware_constraint()
    else:
        print("Tree-aware constraint summaries are not complete; Figure S2 was deferred.")


if __name__ == "__main__":
    main()
