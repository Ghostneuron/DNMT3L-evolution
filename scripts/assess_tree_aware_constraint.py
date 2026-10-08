from __future__ import annotations

from math import inf
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO
from Bio.Data import CodonTable
from scipy.stats import kruskal, mannwhitneyu


PROJECT = Path(__file__).resolve().parents[1]
V2 = PROJECT / "analysis_v2"
SUMMARIES = V2 / "summaries"
SCOPES = ("Placentalia", "Euarchontoglires", "Laurasiatheria", "Sauropsida")
CTERMINAL_START = 178
CTERMINAL_END = 380
REGIONS = {
    "switching_helix": (341, 355),
    "switching_helix_adjacent_loop": (356, 360),
}
STOP_CODONS = set(CodonTable.unambiguous_dna_by_name["Standard"].stop_codons)
FORWARD = CodonTable.unambiguous_dna_by_name["Standard"].forward_table


def bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float)
    if valid.empty:
        return result
    order = valid.sort_values().index
    ranked = valid.loc[order].to_numpy()
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result.loc[order] = np.minimum(adjusted, 1.0)
    return result


def translate_codon(codon: str) -> str | None:
    codon = codon.upper().replace("U", "T")
    if len(codon) != 3 or set(codon) - set("ACGT") or codon in STOP_CODONS:
        return None
    return FORWARD[codon]


def sankoff_minimum_changes(tree, tip_states: dict[str, str | None]) -> int:
    observed = sorted({state for state in tip_states.values() if state is not None})
    if len(observed) <= 1:
        return 0

    def costs(clade) -> dict[str, float]:
        if clade.is_terminal():
            state = tip_states.get(str(clade.name))
            if state is None:
                return {candidate: 0.0 for candidate in observed}
            return {candidate: 0.0 if candidate == state else inf for candidate in observed}
        child_costs = [costs(child) for child in clade.clades]
        return {
            candidate: sum(
                min(value + (child_state != candidate) for child_state, value in child.items())
                for child in child_costs
            )
            for candidate in observed
        }

    return int(min(costs(tree.root).values()))


def main() -> None:
    site_map = pd.read_csv(SUMMARIES / "retained_site_map.csv")
    constraint_sites = pd.read_csv(SUMMARIES / "amino_acid_constraint_by_site.csv")
    rows: list[dict[str, object]] = []

    for scope in SCOPES:
        dataset = f"{scope}_MACSE"
        alignment_path = V2 / "selection_inputs" / f"{dataset}_selection_input.fasta"
        tree_path = V2 / "species_trees" / f"{dataset}_NCBI_taxonomy_tree.nwk"
        sequences = {record.id: str(record.seq).upper() for record in SeqIO.parse(alignment_path, "fasta")}
        tree = Phylo.read(tree_path, "newick")
        tree_tips = {str(tip.name) for tip in tree.get_terminals()}
        if tree_tips != set(sequences):
            missing_tree = sorted(set(sequences) - tree_tips)
            missing_alignment = sorted(tree_tips - set(sequences))
            raise ValueError(
                f"Tree/alignment mismatch for {dataset}: missing from tree={missing_tree[:5]}, "
                f"missing from alignment={missing_alignment[:5]}"
            )
        lengths = {len(sequence) for sequence in sequences.values()}
        if len(lengths) != 1 or next(iter(lengths)) % 3:
            raise ValueError(f"Invalid codon alignment dimensions for {dataset}: {lengths}")

        mapping = site_map.loc[
            (site_map["Scope"] == scope) & (site_map["Aligner"] == "MACSE")
        ].sort_values("Filtered_Site")
        for row in mapping.itertuples(index=False):
            offset = (int(row.Filtered_Site) - 1) * 3
            tip_states = {
                name: translate_codon(sequence[offset : offset + 3])
                for name, sequence in sequences.items()
            }
            resolved = sum(state is not None for state in tip_states.values())
            changes = sankoff_minimum_changes(tree, tip_states)
            rows.append(
                {
                    "Scope": scope,
                    "Human_Position": int(row.Human_Position),
                    "Domain": row.Domain,
                    "Structural_Region": row.Structural_Region,
                    "Resolved_Tips": resolved,
                    "Minimum_AA_Changes": changes,
                    "Minimum_Changes_Per_100_Resolved": (
                        100.0 * changes / resolved if resolved else np.nan
                    ),
                }
            )

    sites = pd.DataFrame(rows).merge(
        constraint_sites[
            [
                "Scope",
                "Human_Position",
                "Aligners_Retaining_Site",
                "All_Three_Resolved_N",
                "All_Three_Resolved_Fraction",
                "All_Three_Exact_Among_Resolved",
                "Alignment_Stable",
            ]
        ],
        on=["Scope", "Human_Position"],
        how="left",
    )
    sites.to_csv(SUMMARIES / "tree_aware_parsimony_constraint_by_site.csv", index=False)

    stable = sites.loc[
        (sites["Aligners_Retaining_Site"] == 3)
        & sites["Alignment_Stable"].astype(str).str.lower().isin({"true", "1"})
    ].copy()
    domain_rows: list[dict[str, object]] = []
    structural_rows: list[dict[str, object]] = []
    for scope, frame in stable.groupby("Scope"):
        groups = {
            domain: values["Minimum_Changes_Per_100_Resolved"].dropna().to_numpy()
            for domain, values in frame.groupby("Domain")
            if len(values) >= 2
        }
        statistic, p_value = kruskal(*groups.values()) if len(groups) >= 2 else (np.nan, np.nan)
        add = groups.get("ADD", np.array([]))
        cterminal = groups.get("C_terminal_MTase_like", np.array([]))
        pair_stat, pair_p = (
            mannwhitneyu(add, cterminal, alternative="two-sided")
            if len(add) and len(cterminal)
            else (np.nan, np.nan)
        )
        domain_rows.append(
            {
                "Scope": scope,
                "ADD_Sites": len(add),
                "C_terminal_Sites": len(cterminal),
                "ADD_Median_Changes_Per_100": np.median(add) if len(add) else np.nan,
                "C_terminal_Median_Changes_Per_100": np.median(cterminal) if len(cterminal) else np.nan,
                "Kruskal_All_Domains_p": p_value,
                "ADD_vs_C_terminal_U": pair_stat,
                "ADD_vs_C_terminal_p": pair_p,
            }
        )

        cterminal_frame = frame.loc[frame["Domain"] == "C_terminal_MTase_like"]
        for region in ("switching_helix", "switching_helix_adjacent_loop"):
            focal = cterminal_frame.loc[
                cterminal_frame["Structural_Region"] == region,
                "Minimum_Changes_Per_100_Resolved",
            ].dropna().to_numpy()
            background = cterminal_frame.loc[
                cterminal_frame["Structural_Region"] == "other",
                "Minimum_Changes_Per_100_Resolved",
            ].dropna().to_numpy()
            region_stat, region_p = (
                mannwhitneyu(focal, background, alternative="two-sided")
                if len(focal) and len(background)
                else (np.nan, np.nan)
            )
            structural_rows.append(
                {
                    "Scope": scope,
                    "Structural_Region": region,
                    "Region_Sites": len(focal),
                    "Background_Sites": len(background),
                    "Region_Median_Changes_Per_100": np.median(focal) if len(focal) else np.nan,
                    "Background_Median_Changes_Per_100": np.median(background) if len(background) else np.nan,
                    "Mann_Whitney_U": region_stat,
                    "p": region_p,
                }
            )

    domains = pd.DataFrame(domain_rows).sort_values("Scope")
    domains["ADD_vs_C_terminal_q_across_scopes"] = bh_adjust(domains["ADD_vs_C_terminal_p"])
    domains.to_csv(SUMMARIES / "tree_aware_parsimony_domain_tests.csv", index=False)

    structural = pd.DataFrame(structural_rows).sort_values(["Scope", "Structural_Region"])
    structural["BH_q_across_scope_region_tests"] = bh_adjust(structural["p"])
    structural.to_csv(SUMMARIES / "tree_aware_parsimony_structural_tests.csv", index=False)

    window_rows: list[dict[str, object]] = []
    for scope, frame in stable.groupby("Scope"):
        cterminal = frame.loc[
            frame["Human_Position"].between(CTERMINAL_START, CTERMINAL_END)
        ].set_index("Human_Position")["Minimum_Changes_Per_100_Resolved"]
        for region, (start, end) in REGIONS.items():
            length = end - start + 1
            focal_positions = list(range(start, end + 1))
            focal_complete = all(position in cterminal.index for position in focal_positions)
            focal_median = float(cterminal.loc[focal_positions].median()) if focal_complete else np.nan
            windows: list[float] = []
            for window_start in range(CTERMINAL_START, CTERMINAL_END - length + 2):
                positions = list(range(window_start, window_start + length))
                if all(position in cterminal.index for position in positions):
                    windows.append(float(cterminal.loc[positions].median()))
            at_least_as_variable = (
                int(sum(value >= focal_median for value in windows))
                if focal_complete and windows
                else np.nan
            )
            empirical_p = (
                (1 + at_least_as_variable) / (1 + len(windows))
                if focal_complete and windows
                else np.nan
            )
            window_rows.append(
                {
                    "Scope": scope,
                    "Structural_Region": region,
                    "Region_Length": length,
                    "Focal_Window_Complete": focal_complete,
                    "Focal_Median_Changes_Per_100": focal_median,
                    "Eligible_Comparable_Windows": len(windows),
                    "Windows_At_Least_As_Variable": at_least_as_variable,
                    "Empirical_Upper_Tail_p": empirical_p,
                }
            )
    windows = pd.DataFrame(window_rows).sort_values(["Scope", "Structural_Region"])
    windows["BH_q_across_evaluable_tests"] = bh_adjust(windows["Empirical_Upper_Tail_p"])
    windows.to_csv(SUMMARIES / "tree_aware_parsimony_window_sensitivity.csv", index=False)

    print(domains.to_string(index=False))
    print("\nStructural-region sensitivity")
    print(structural.to_string(index=False))
    print("\nTree-aware same-length-window sensitivity")
    print(windows.to_string(index=False))


if __name__ == "__main__":
    main()
