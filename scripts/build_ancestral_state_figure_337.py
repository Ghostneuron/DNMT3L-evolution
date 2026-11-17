#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D


PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results"
EXTENDED = PROJECT / "analysis/regulatory_coevolution_ptm_20261104/results"
OUT_PNG = PROJECT / "figures/Figure_4_ancestral_state_context.png"
OUT_PDF = PROJECT / "figures/Figure_4_ancestral_state_context.pdf"


def node_state(table: pd.DataFrame, lineage: str) -> tuple[str, float, float]:
    tax = table.loc[(table["topology"] == "taxonomy") & (table["lineage"] == lineage)].iloc[0]
    gene = table.loc[(table["topology"] == "gene_tree") & (table["lineage"] == lineage)].iloc[0]
    return str(tax["maximum_posterior_state"]), float(tax["maximum_posterior"]), float(gene["maximum_posterior"])


def draw_tree(ax: plt.Axes, nodes: pd.DataFrame) -> None:
    y = {"Amphibia": 4.6, "Sauropsida": 3.45, "Marsupialia": 2.15, "Placentalia": 1.05}
    x_root, x_amniota, x_theria, x_tip = 0.25, 1.55, 2.85, 4.3
    y_theria = (y["Marsupialia"] + y["Placentalia"]) / 2
    y_amniota = (y["Sauropsida"] + y_theria) / 2
    y_root = (y["Amphibia"] + y_amniota) / 2

    ax.plot([x_root, x_root], [y_amniota, y["Amphibia"]], color="#334155", lw=1.8)
    ax.plot([x_root, x_tip], [y["Amphibia"], y["Amphibia"]], color="#334155", lw=1.8)
    ax.plot([x_root, x_amniota], [y_amniota, y_amniota], color="#334155", lw=1.8)
    ax.plot([x_amniota, x_amniota], [y_theria, y["Sauropsida"]], color="#334155", lw=1.8)
    ax.plot([x_amniota, x_tip], [y["Sauropsida"], y["Sauropsida"]], color="#334155", lw=1.8)
    ax.plot([x_amniota, x_theria], [y_theria, y_theria], color="#334155", lw=1.8)
    ax.plot([x_theria, x_theria], [y["Placentalia"], y["Marsupialia"]], color="#334155", lw=1.8)
    ax.plot([x_theria, x_tip], [y["Marsupialia"], y["Marsupialia"]], color="#334155", lw=1.8)
    ax.plot([x_theria, x_tip], [y["Placentalia"], y["Placentalia"]], color="#334155", lw=1.8)

    positions = {
        "Amphibia": (x_root, y["Amphibia"]),
        "Amniota": (x_amniota, y_amniota),
        "Sauropsida": (x_amniota, y["Sauropsida"]),
        "Theria": (x_theria, y_theria),
        "Marsupialia": (x_theria, y["Marsupialia"]),
        "Placentalia": (x_theria, y["Placentalia"]),
    }
    colors = {"K": "#0f766e", "Q": "#2563eb"}
    for lineage, (xv, yv) in positions.items():
        state, pp_tax, pp_gene = node_state(nodes, lineage)
        ax.scatter([xv], [yv], s=115, color=colors.get(state, "#64748b"), edgecolor="white", linewidth=1.2, zorder=4)
        dx = -0.08 if lineage in {"Amphibia", "Sauropsida", "Marsupialia", "Placentalia"} else 0.08
        ha = "right" if dx < 0 else "left"
        display_lineage = "Therian node" if lineage == "Theria" else lineage
        ax.text(xv + dx, yv + 0.18, f"{display_lineage}  {state}", ha=ha, va="bottom", fontsize=9.0, fontweight="bold")
        ax.text(xv + dx, yv - 0.17, f"{pp_tax:.3f} / {pp_gene:.3f}", ha=ha, va="top", fontsize=7.4, color="#475569")

    ax.annotate(
        "K-to-Q node-state difference\n(monotreme unsampled)",
        xy=((x_amniota + x_theria) / 2, y_theria),
        xytext=((x_amniota + x_theria) / 2, y_theria + 0.52),
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
        color="#b45309",
        arrowprops=dict(arrowstyle="-|>", color="#b45309", lw=1.3),
    )
    ax.text(0.05, 0.1, "Posterior values: taxonomy / gene-tree topology", transform=ax.transAxes, fontsize=7.7, color="#475569")
    ax.set_xlim(-0.25, 4.8)
    ax.set_ylim(0.45, 5.1)
    ax.axis("off")
    ax.set_title("a  NP349 alignment-sensitive ancestral reconstruction", loc="left", fontsize=11.5, fontweight="bold", pad=7)


def draw_substitutions(ax: plt.Axes, inventory: pd.DataFrame) -> None:
    rows = inventory.loc[inventory["robust_therian_stem_change"].astype(bool)].copy()
    rows = rows.sort_values("human_np_position")
    add = rows.loc[rows["add_domain"].astype(bool)]
    cterm = rows.loc[rows["c_terminal"].astype(bool)]

    ax.axis("off")
    ax.set_title("b  Robust Amniota-to-therian-node substitutions", loc="left", fontsize=11.5, fontweight="bold", pad=7)
    ax.text(0.03, 0.91, "ADD domain  6", fontsize=9.5, fontweight="bold", color="#0f766e", transform=ax.transAxes)
    ax.text(0.52, 0.91, "C-terminal region  12", fontsize=9.5, fontweight="bold", color="#7c3aed", transform=ax.transAxes)

    def render(group: pd.DataFrame, x: float, y0: float) -> None:
        for i, row in enumerate(group.itertuples(index=False)):
            y = y0 - 0.057 * i
            contact = bool(row.interface_union_5A)
            color = "#b45309" if contact else "#1f2937"
            weight = "bold" if contact else "normal"
            ax.text(x, y, f"NP{int(row.human_np_position)}  {row.state_pair}", transform=ax.transAxes, fontsize=8.5, color=color, fontweight=weight)
            if contact:
                ax.scatter([x - 0.018], [y + 0.008], transform=ax.transAxes, s=25, color="#f59e0b", edgecolor="white", linewidth=0.6)

    render(add, 0.05, 0.84)
    render(cterm, 0.54, 0.84)
    ax.text(0.03, 0.015, "Orange marks DNMT3A- or nucleosome-contacting positions in PDB 9MPP.", transform=ax.transAxes, fontsize=7.7, color="#475569")


def draw_interface(ax: plt.Axes, summary: pd.DataFrame) -> None:
    ax.axis("off")
    ax.set_title("c  Matched DNMT3A-DNMT3L interface states", loc="left", fontsize=11.5, fontweight="bold", pad=7)
    ax.text(0.02, 0.88, "265 species with matched orthologues", transform=ax.transAxes, fontsize=9.5, fontweight="bold")

    labels = {
        "DNMT3A1_725_state": "DNMT3A E725",
        "DNMT3A1_726_state": "DNMT3A G726",
        "DNMT3A1_729_state": "DNMT3A R729",
        "DNMT3A1_732_state": "DNMT3A F732",
        "DNMT3A1_772_state": "DNMT3A F772",
        "DNMT3L_NP_261_state": "DNMT3L NP261",
        "DNMT3L_NP_291_state": "DNMT3L NP291",
        "DNMT3L_NP_293_state": "DNMT3L NP293",
    }
    y = 0.78
    for key, label in labels.items():
        row = summary.loc[summary["gene_coordinate"] == key].iloc[0]
        n_states = int(row["distinct_resolved_states"])
        counts = str(row["state_counts"])
        color = "#0f766e" if key.startswith("DNMT3A") else "#7c3aed"
        ax.text(0.03, y, label, transform=ax.transAxes, fontsize=8.5, color=color, fontweight="bold")
        ax.text(0.34, y, f"{n_states} state" + ("" if n_states == 1 else "s"), transform=ax.transAxes, fontsize=8.3)
        ax.text(0.53, y, counts, transform=ax.transAxes, fontsize=8.0, color="#334155")
        y -= 0.085
    ax.text(0.03, 0.04, "DNMT3A contact states are invariant in this sampled panel; DNMT3L states vary.", transform=ax.transAxes, fontsize=7.7, color="#475569")


def main() -> None:
    nodes = pd.read_csv(RESULTS / "switching_helix_337_asr_major_nodes.tsv", sep="\t")
    inventory = pd.read_csv(RESULTS / "therian_stem_robust_substitutions.tsv", sep="\t")
    interface = pd.read_csv(EXTENDED / "dnmt3a_dnmt3l_matched_interface_summary_vertebrate.tsv", sep="\t")

    plt.rcParams.update({"font.family": "Arial", "axes.edgecolor": "#334155", "text.color": "#111827"})
    fig = plt.figure(figsize=(12, 9.7), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[0.85, 1.15], width_ratios=[1.06, 0.94])
    draw_tree(fig.add_subplot(grid[0, :]), nodes)
    draw_substitutions(fig.add_subplot(grid[1, 0]), inventory)
    draw_interface(fig.add_subplot(grid[1, 1]), interface)
    fig.suptitle("DNMT3L ancestral sequence change and matched interface context", fontsize=15, fontweight="bold")

    legend = [
        Line2D([0], [0], marker="o", color="w", label="K", markerfacecolor="#0f766e", markersize=8),
        Line2D([0], [0], marker="o", color="w", label="Q", markerfacecolor="#2563eb", markersize=8),
        Line2D([0], [0], marker="o", color="w", label="Interface-proximal substitution", markerfacecolor="#f59e0b", markersize=8),
    ]
    fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, -0.012), ncol=3, frameon=False, fontsize=8.5)

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(OUT_PNG)
    print(OUT_PDF)


if __name__ == "__main__":
    main()
