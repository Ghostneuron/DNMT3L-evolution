#!/usr/bin/env python3
"""Plot the accepted DNMT3L G9a-contact triplet states on the rodent taxonomy tree."""
from __future__ import annotations

import copy
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo


PROJECT = Path(__file__).resolve().parents[1]
ATLAS = PROJECT / "results/functional_site_states_by_record.tsv"
TREE = (
    PROJECT
    / "analysis/selection_reanalysis_337/trees/Placentalia_CONSENSUS_selection_ncbi_taxonomy_topology.nwk"
)
OUTDIR = PROJECT / "figures/current_v12"
STATE_OUTPUT = PROJECT / "results/g9a_rodent_tree_states.tsv"

RODENTIA_NODE = "NCBI_9989"
MUROIDEA_NODE = "NCBI_337687"
REQUIRED_POSITIONS = [105, 106, 107]

MURIDAE = {
    "Acomys minous", "Acomys russatus", "Apodemus sylvaticus", "Arvicanthis niloticus",
    "Mastomys coucha", "Meriones unguiculatus", "Micromys minutus", "Mus musculus",
    "Mus pahari", "Psammomys obesus", "Rattus norvegicus", "Rattus rattus",
}
CRICETIDAE = {
    "Arvicola amphibius", "Chionomys nivalis", "Cricetulus griseus", "Mesocricetus auratus",
    "Microtus ochrogaster", "Microtus oregoni", "Microtus pennsylvanicus", "Myodes glareolus",
    "Onychomys torridus", "Peromyscus californicus insignis", "Peromyscus eremicus",
    "Peromyscus leucopus", "Peromyscus maniculatus bairdii", "Phodopus roborovskii",
}
SPALACIDAE = {"Nannospalax galili"}

COLORS = {"LFI": "#D95F02", "LLI": "#2C7FB8"}


def accepted_triplets() -> pd.DataFrame:
    data = pd.read_csv(ATLAS, sep="\t")
    data = data.loc[data["feature_id"].str.startswith("G9A_INTERFACE")].copy()
    pivot = data.pivot(
        index=["safe_id", "species", "broad_lineage"],
        columns="human_np_position",
        values=["query_residue", "mapping_status"],
    )
    accepted = pivot["mapping_status"][REQUIRED_POSITIONS].eq("PASS").all(axis=1)
    states = pivot.loc[accepted, "query_residue"][REQUIRED_POSITIONS].astype(str).agg("".join, axis=1)
    frame = pivot.loc[accepted].index.to_frame(index=False)
    frame["triplet_state"] = states.to_numpy()
    return frame


def family(species: str) -> str:
    if species in MURIDAE:
        return "Muridae"
    if species in CRICETIDAE:
        return "Cricetidae"
    if species in SPALACIDAE:
        return "Spalacidae"
    return "non-muroid Rodentia"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    tree = Phylo.read(TREE, "newick")
    rodent_clade = next(clade for clade in tree.find_clades() if clade.name == RODENTIA_NODE)
    rodent_tree = Phylo.BaseTree.Tree(root=copy.deepcopy(rodent_clade), rooted=True)
    terminals = rodent_tree.get_terminals()
    rodent_ids = {terminal.name for terminal in terminals}

    states = accepted_triplets()
    states = states.loc[states["safe_id"].isin(rodent_ids)].copy()
    if len(terminals) != 43 or len(states) != 43:
        raise ValueError(f"Expected 43 rodent tips with accepted triplets, got {len(terminals)} and {len(states)}")
    state_by_id = states.set_index("safe_id")["triplet_state"].to_dict()
    species_by_id = states.set_index("safe_id")["species"].to_dict()
    if set(state_by_id) != rodent_ids:
        raise ValueError("Rodent tree tips and accepted triplet records do not match")

    states["family"] = states["species"].map(family)
    states["muroid"] = states["family"].ne("non-muroid Rodentia")
    if not states.loc[states["muroid"], "triplet_state"].eq("LFI").all():
        raise ValueError("A sampled muroid lacks LFI")
    if not states.loc[~states["muroid"], "triplet_state"].eq("LLI").all():
        raise ValueError("A sampled non-muroid rodent lacks LLI")
    states.sort_values(["muroid", "family", "species"]).to_csv(STATE_OUTPUT, sep="\t", index=False)

    y_by_terminal = {terminal: len(terminals) - 1 - index for index, terminal in enumerate(terminals)}
    y_by_clade: dict[object, float] = {}

    def assign_y(clade) -> float:
        if clade.is_terminal():
            y_by_clade[clade] = y_by_terminal[clade]
        else:
            child_y = [assign_y(child) for child in clade.clades]
            y_by_clade[clade] = (min(child_y) + max(child_y)) / 2
        return y_by_clade[clade]

    assign_y(rodent_tree.root)
    depths: dict[object, int] = {}

    def assign_depth(clade, depth: int = 0) -> None:
        depths[clade] = depth
        for child in clade.clades:
            assign_depth(child, depth + 1)

    assign_depth(rodent_tree.root)
    max_depth = max(depths.values())

    fig, ax = plt.subplots(figsize=(11.7, 12.2))
    fig.subplots_adjust(left=0.055, right=0.985, top=0.94, bottom=0.085)
    line_color = "#475569"
    for clade in rodent_tree.find_clades(order="preorder"):
        if clade.is_terminal():
            continue
        child_y = [y_by_clade[child] for child in clade.clades]
        ax.plot([depths[clade], depths[clade]], [min(child_y), max(child_y)], color=line_color, lw=0.8)
        for child in clade.clades:
            ax.plot(
                [depths[clade], depths[child]],
                [y_by_clade[child], y_by_clade[child]],
                color=line_color,
                lw=0.8,
            )

    muroid = next(clade for clade in rodent_tree.find_clades() if clade.name == MUROIDEA_NODE)
    muroid_y = [y_by_clade[tip] for tip in muroid.get_terminals()]
    ax.axhspan(min(muroid_y) - 0.45, max(muroid_y) + 0.45, color="#FFF3E6", zorder=-5)
    ax.text(
        0.05,
        max(muroid_y) + 0.62,
        "Muroidea: 27/27 LFI",
        color=COLORS["LFI"],
        fontsize=10.5,
        fontweight="bold",
        va="bottom",
        zorder=20,
        bbox={"facecolor": "#FFF9F2", "edgecolor": "none", "pad": 1.5},
    )

    label_x = max_depth + 0.25
    state_x = max_depth + 4.55
    family_x = max_depth + 5.75
    for terminal in terminals:
        y = y_by_clade[terminal]
        safe_id = terminal.name
        species = species_by_id[safe_id]
        triplet = state_by_id[safe_id]
        fam = family(species)
        ax.text(label_x, y, species, fontsize=8.2, fontstyle="italic", va="center", color="#111827")
        ax.scatter(state_x, y, s=34, marker="s", color=COLORS[triplet], edgecolor="white", lw=0.45, zorder=4)
        ax.text(state_x + 0.18, y, triplet, fontsize=8.1, fontweight="bold", va="center", color=COLORS[triplet])
        if fam != "non-muroid Rodentia":
            ax.text(family_x, y, fam, fontsize=7.2, va="center", color="#6B7280")

    ax.text(label_x, len(terminals) + 0.2, "Species", fontsize=9, fontweight="bold", color="#374151")
    ax.text(state_x - 0.05, len(terminals) + 0.2, "Triplet", fontsize=9, fontweight="bold", color="#374151")
    ax.text(family_x, len(terminals) + 0.2, "Muroid family", fontsize=9, fontweight="bold", color="#374151")

    ax.text(
        0.02,
        -1.25,
        "All 16 sampled non-muroid rodents carried LLI. Across the full eutherian panel, LFI occurred in 27/27 "
        "sampled muroids and 0/192 other species\nwith accepted mappings.",
        fontsize=9.2,
        color="#374151",
        va="top",
    )
    ax.set_xlim(-0.25, family_x + 1.35)
    ax.set_ylim(-2.35, len(terminals) + 0.8)
    ax.axis("off")
    fig.suptitle(
        "A single G9a-contact substitution distinguishes sampled muroids from other rodents",
        fontsize=15.5,
        fontweight="bold",
        y=0.978,
    )
    ax.set_title(
        "DNMT3L triplet state at NP105/NP106/NP107 on the NCBI Taxonomy topology",
        fontsize=11,
        color="#475569",
        pad=12,
    )

    png = OUTDIR / "Figure_1_muroid_g9a_triplet_tree.png"
    pdf = OUTDIR / "Figure_1_muroid_g9a_triplet_tree.pdf"
    fig.savefig(png, dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(png)
    print(pdf)
    print(STATE_OUTPUT)


if __name__ == "__main__":
    main()
