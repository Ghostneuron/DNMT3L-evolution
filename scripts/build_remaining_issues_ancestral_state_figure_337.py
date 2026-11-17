#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "scripts/build_ancestral_state_figure_337.py"
SPEC = importlib.util.spec_from_file_location("ancestral_state_v09", SOURCE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not load {SOURCE}")
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)

OUTDIR = PROJECT / "figures/current_v09"
BASE.OUT_PNG = OUTDIR / "Figure_4_ancestral_state_context_corrected.png"
BASE.OUT_PDF = OUTDIR / "Figure_4_ancestral_state_context_corrected.pdf"


def draw_tree_without_overlap(ax, nodes) -> None:
    y = {"Amphibia": 4.6, "Sauropsida": 3.45, "Marsupialia": 2.15, "Placentalia": 1.05}
    x_root, x_amniota, x_theria, x_tip = 0.25, 1.55, 2.85, 4.3
    y_theria = (y["Marsupialia"] + y["Placentalia"]) / 2
    y_amniota = (y["Sauropsida"] + y_theria) / 2

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
    text_positions = {
        "Amphibia": (x_root - 0.08, y["Amphibia"] + 0.18, y["Amphibia"] - 0.17, "right"),
        "Sauropsida": (x_amniota - 0.08, y["Sauropsida"] + 0.18, y["Sauropsida"] - 0.17, "right"),
        "Amniota": (x_amniota - 0.08, y_amniota + 0.18, y_amniota - 0.17, "right"),
        "Marsupialia": (x_theria + 0.08, y["Marsupialia"] + 0.37, y["Marsupialia"] + 0.17, "left"),
        "Theria": (x_theria - 0.08, y_theria + 0.16, y_theria - 0.14, "right"),
        "Placentalia": (x_theria + 0.08, y["Placentalia"] + 0.16, y["Placentalia"] - 0.17, "left"),
    }
    for lineage, (xv, yv) in positions.items():
        state, pp_tax, pp_gene = BASE.node_state(nodes, lineage)
        ax.scatter([xv], [yv], s=115, color=colors.get(state, "#64748b"), edgecolor="white", linewidth=1.2, zorder=4)
        text_x, state_y, posterior_y, ha = text_positions[lineage]
        display_lineage = "Therian node" if lineage == "Theria" else lineage
        ax.text(text_x, state_y, f"{display_lineage}  {state}", ha=ha, va="bottom", fontsize=9.0, fontweight="bold")
        ax.text(text_x, posterior_y, f"{pp_tax:.3f} / {pp_gene:.3f}", ha=ha, va="top", fontsize=7.4, color="#475569")

    ax.annotate(
        "K-to-Q node-state difference\n(monotreme unsampled)",
        xy=((x_amniota + x_theria) / 2, y_theria),
        xytext=((x_amniota + x_theria) / 2, y_theria + 0.72),
        ha="center",
        va="center",
        fontsize=9.3,
        fontweight="bold",
        color="#b45309",
        arrowprops=dict(arrowstyle="-|>", color="#b45309", lw=1.3),
    )
    ax.text(
        0.05,
        0.018,
        "Mapping caveat: all four amphibian NP349 calls fail the local-context identity gate;\n"
        "Geotrypetes K lies beside a seven-residue deletion and two masked residues, and its omission lowers Amniota K to 0.518 / 0.499.",
        transform=ax.transAxes,
        fontsize=6.9,
        color="#991B1B",
        fontweight="bold",
        va="bottom",
    )
    ax.text(
        0.70,
        0.92,
        "Posterior values: taxonomy / gene-tree topology",
        transform=ax.transAxes,
        fontsize=7.7,
        color="#475569",
        ha="center",
    )
    ax.set_xlim(-0.25, 4.8)
    ax.set_ylim(-0.05, 5.1)
    ax.axis("off")
    ax.set_title("a  Exploratory NP349 reconstruction depends on Geotrypetes", loc="left", fontsize=11.5, fontweight="bold", pad=7)


BASE.draw_tree = draw_tree_without_overlap


if __name__ == "__main__":
    OUTDIR.mkdir(parents=True, exist_ok=True)
    BASE.main()
