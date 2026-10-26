#!/usr/bin/env python3
"""Build a selection-summary figure that reflects the revised evidence hierarchy."""

from pathlib import Path
import math

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
SUMMARIES = ROOT / "analysis_v2" / "summaries"
OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(parents=True, exist_ok=True)

SCOPE_ORDER = ["Placentalia", "Euarchontoglires", "Laurasiatheria", "Sauropsida"]
ALIGNER_ORDER = ["MACSE", "MAFFT", "PRANK"]

INK = "#1f2933"
MUTED = "#718096"
BLUE = "#3f7cac"
TEAL = "#2a9d8f"
CORAL = "#cf4f4f"
PURPLE = "#7566a8"
GOLD = "#e6aa35"
PALE = "#f3f5f7"


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def panel_label(ax, label: str) -> None:
    ax.text(-0.12, 1.04, label, transform=ax.transAxes, fontsize=16,
            fontweight="bold", va="bottom", ha="left", color=INK)


def style_axes(ax, grid_axis: str = "y") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(colors=INK, labelcolor=INK)
    ax.grid(axis=grid_axis, color="#d8dee5", linewidth=0.65, alpha=0.75)
    ax.set_axisbelow(True)


def main() -> None:
    robustness = pd.read_csv(SUMMARIES / "site_robustness_by_scope.csv")
    busted = pd.read_csv(SUMMARIES / "BUSTED_gene_wide_summary.csv")
    error_sink = pd.read_csv(SUMMARIES / "BUSTED_error_sink_summary.csv")
    multihit = pd.read_csv(SUMMARIES / "BUSTED_multihit_moderate_risk_summary.csv")
    fel_mnv = pd.read_csv(SUMMARIES / "FEL_multihit_primary_sites_summary.csv")
    meme_mnv = pd.read_csv(SUMMARIES / "MEME_multihit_primary_sites_summary.csv")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": INK,
    })
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.6), layout="constrained")

    # Panel a: make small diversifying classes legible without hiding purifying dominance.
    ax = axes[0, 0]
    panel_label(ax, "a")
    ax.set_title("Alignment-robust site classes", loc="left", fontweight="bold")
    categories = [
        ("Robust_Pervasive_q10", "Pervasive\nq <= 0.10"),
        ("Robust_Episodic_q10", "Episodic\nq <= 0.10"),
        ("Robust_Purifying_q05", "Purifying\nq <= 0.05"),
    ]
    matrix = np.zeros((len(SCOPE_ORDER), len(categories)), dtype=int)
    for i, scope in enumerate(SCOPE_ORDER):
        rows = robustness.loc[robustness["Scope"] == scope]
        for j, (column, _) in enumerate(categories):
            matrix[i, j] = int(as_bool(rows[column]).sum())
    heat = ax.imshow(np.log10(matrix + 1), aspect="auto",
                     cmap=LinearSegmentedColormap.from_list("counts", [PALE, GOLD, CORAL]))
    ax.set_xticks(range(len(categories)), [label for _, label in categories])
    ax.set_yticks(range(len(SCOPE_ORDER)), SCOPE_ORDER)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, str(matrix[i, j]), ha="center", va="center",
                    fontsize=10, fontweight="bold",
                    color="white" if np.log10(matrix[i, j] + 1) > 1.6 else INK)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(heat, ax=ax, shrink=0.76, pad=0.025)
    cb.set_label("log10(count + 1)")

    # Panel b: distinguish standard alignment agreement from model sensitivity.
    ax = axes[0, 1]
    panel_label(ax, "b")
    ax.set_title("Gene-wide episodic-selection model tiers", loc="left", fontweight="bold")
    cap = 16.0
    x = np.arange(len(SCOPE_ORDER))
    offsets = {"MACSE": -0.18, "MAFFT": 0.0, "PRANK": 0.18}
    colors = {"MACSE": BLUE, "MAFFT": TEAL, "PRANK": CORAL}
    for aligner in ALIGNER_ORDER:
        vals = []
        for scope in SCOPE_ORDER:
            hit = busted.loc[(busted["Scope"] == scope) & (busted["Aligner"] == aligner)]
            q = float(hit.iloc[0]["BUSTED_q_across_12"])
            vals.append(min(-math.log10(max(q, 10 ** -cap)), cap))
        ax.scatter(x + offsets[aligner], vals, s=45, color=colors[aligner],
                   edgecolor="white", linewidth=0.7, label=f"Standard {aligner}", zorder=3)
    for i, scope in enumerate(SCOPE_ORDER):
        hit = error_sink.loc[error_sink["Scope"] == scope]
        if not hit.empty:
            q = float(hit.iloc[0]["BUSTED_ErrorSink_q"])
            ax.scatter(i - 0.07, min(-math.log10(max(q, 10 ** -cap)), cap),
                       marker="D", s=42, facecolor="none", edgecolor=PURPLE,
                       linewidth=1.4, zorder=4)
        hit = multihit.loc[multihit["Scope"] == scope]
        if not hit.empty:
            q_col = next(c for c in hit.columns if "q" in c.lower() and "BUSTED" in c)
            q = float(hit.iloc[0][q_col])
            ax.scatter(i + 0.07, min(-math.log10(max(q, 10 ** -cap)), cap),
                       marker="s", s=42, facecolor="none", edgecolor=INK,
                       linewidth=1.4, zorder=4)
    ax.axhline(-math.log10(0.05), color=INK, linestyle="--", linewidth=1)
    ax.set_xticks(x, SCOPE_ORDER, rotation=15, ha="right")
    ax.set_ylabel("-log10(corrected q), capped at 16")
    ax.set_ylim(-0.25, 16.7)
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=colors[a],
               markeredgecolor="white", markersize=7, label=f"Standard {a}")
        for a in ALIGNER_ORDER
    ] + [
        Line2D([0], [0], marker="D", color="none", markerfacecolor="none",
               markeredgecolor=PURPLE, markersize=7, label="MACSE error sink"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="none",
               markeredgecolor=INK, markersize=7, label="MACSE multinucleotide"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=7, ncol=2, loc="upper right")
    style_axes(ax, grid_axis="y")

    # Panel c: the most important sitewise robustness result.
    ax = axes[1, 0]
    panel_label(ax, "c")
    ax.set_title("Strict sitewise candidates under multinucleotide models", loc="left", fontweight="bold")
    # Restrict the conditional tables to strict primary scope-site candidates. h358 was
    # also evaluated in every scope as a prespecified sensitivity and must not inflate
    # the denominator unless it crossed the primary strict threshold in that scope.
    def retained_count(df: pd.DataFrame, primary_column: str) -> tuple[int, int]:
        candidates = robustness.loc[
            as_bool(robustness[primary_column]), ["Scope", "Human_Position"]
        ].drop_duplicates()
        evaluated = candidates.merge(df, on=["Scope", "Human_Position"], how="left")
        groups = evaluated.groupby(["Scope", "Human_Position"], dropna=False)
        kept = sum(bool(as_bool(grp["Conditional_Robust_2of3"]).any()) for _, grp in groups)
        return len(candidates), int(kept)

    pervasive_total, pervasive_kept = retained_count(fel_mnv, "Robust_Pervasive_q05")
    episodic_total, episodic_kept = retained_count(meme_mnv, "Robust_Episodic_q05")
    totals = [pervasive_total, episodic_total]
    kept = [pervasive_kept, episodic_kept]
    labels = ["Pervasive\n(FEL)", "Episodic\n(MEME)"]
    xx = np.arange(2)
    ax.bar(xx, totals, width=0.56, color="#dce3e9", label="Strict primary candidates")
    ax.bar(xx, kept, width=0.56, color=[TEAL, PURPLE], label="Retained in >=2 aligners")
    for xi, total, retain in zip(xx, totals, kept):
        ax.text(xi, total + 0.2, f"{retain}/{total} retained", ha="center", fontweight="bold")
    ax.set_xticks(xx, labels)
    ax.set_ylabel("Scope-site candidates")
    ax.set_ylim(0, max(totals) + 1.5)
    ax.legend(handles=[Patch(facecolor="#dce3e9", label="Strict primary candidates"),
                       Patch(facecolor=TEAL, label="Retained after matched model")],
              frameon=False, fontsize=8, loc="upper left")
    style_axes(ax, grid_axis="y")

    # Panel d: show candidate positions and scopes; h358 is one labeled point.
    ax = axes[1, 1]
    panel_label(ax, "d")
    ax.set_title("Exploratory diversifying candidates along DNMT3L", loc="left", fontweight="bold")
    ax.axvspan(41, 173, color="#dceee8", alpha=0.75, label="ADD")
    ax.axvspan(178, 380, color="#e8edf5", alpha=0.75, label="C-terminal region")
    ax.axvspan(341, 355, color=GOLD, alpha=0.35, label="Switching helix")
    ax.axvspan(356, 360, color=CORAL, alpha=0.25, label="Adjacent loop")
    for yi, scope in enumerate(SCOPE_ORDER):
        rows = robustness.loc[robustness["Scope"] == scope]
        per = rows.loc[as_bool(rows["Robust_Pervasive_q10"]), "Human_Position"].astype(int)
        epi = rows.loc[as_bool(rows["Robust_Episodic_q10"]), "Human_Position"].astype(int)
        ax.scatter(per, np.full(len(per), yi - 0.11), marker="o", s=42,
                   color=TEAL, edgecolor="white", linewidth=0.6, zorder=3)
        ax.scatter(epi, np.full(len(epi), yi + 0.11), marker="^", s=48,
                   color=PURPLE, edgecolor="white", linewidth=0.6, zorder=3)
        if 358 in set(per) | set(epi):
            yv = yi - 0.11 if 358 in set(per) else yi + 0.11
            ax.scatter([358], [yv], s=105, facecolor="none", edgecolor=INK,
                       linewidth=1.1, zorder=4)
    ax.annotate("h358", xy=(358, -0.11), xytext=(322, -0.72),
                arrowprops=dict(arrowstyle="-", color=INK, lw=0.8), fontsize=7)
    ax.set_xlim(0, 390)
    ax.set_ylim(-0.85, len(SCOPE_ORDER) - 0.55)
    ax.set_yticks(range(len(SCOPE_ORDER)), SCOPE_ORDER)
    ax.set_xlabel("Human reference position (NP_037501.2)")
    ax.invert_yaxis()
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=TEAL,
               markeredgecolor="white", markersize=7, label="Pervasive q <= 0.10"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=PURPLE,
               markeredgecolor="white", markersize=7, label="Episodic q <= 0.10"),
    ], frameon=False, fontsize=7, loc="upper left")
    style_axes(ax, grid_axis="x")

    stem = OUT / "Figure_4_selection_evidence_revised"
    fig.savefig(stem.with_suffix(".png"), dpi=350, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
