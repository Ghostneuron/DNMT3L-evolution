from __future__ import annotations

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu


PROJECT = Path(__file__).resolve().parents[1]
SUMMARIES = PROJECT / "analysis_v2" / "summaries"


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


def main() -> None:
    sites = pd.read_csv(SUMMARIES / "retained_site_map.csv")
    concordance = pd.read_csv(SUMMARIES / "cross_aligner_concordance_by_site.csv")

    alignment_domains = (
        sites.loc[sites["Domain"].isin(["ADD", "C_terminal_MTase_like"])]
        .groupby(["Scope", "Aligner", "Domain"])
        .agg(
            Retained_Sites=("Human_Position", "size"),
            Median_AA_Entropy=("Amino_Acid_Shannon_Entropy", "median"),
        )
        .reset_index()
    )
    medians = alignment_domains.pivot(
        index=["Scope", "Aligner"], columns="Domain", values="Median_AA_Entropy"
    )
    counts = alignment_domains.pivot(
        index=["Scope", "Aligner"], columns="Domain", values="Retained_Sites"
    )
    alignment_contrasts = medians.join(counts, lsuffix="_Median", rsuffix="_Sites").reset_index()
    alignment_contrasts = alignment_contrasts.rename(
        columns={
            "ADD_Median": "ADD_Median_AA_Entropy",
            "C_terminal_MTase_like_Median": "C_terminal_Median_AA_Entropy",
            "ADD_Sites": "ADD_Retained_Sites",
            "C_terminal_MTase_like_Sites": "C_terminal_Retained_Sites",
        }
    )
    alignment_contrasts["ADD_to_C_terminal_Median_Ratio"] = (
        alignment_contrasts["ADD_Median_AA_Entropy"]
        / alignment_contrasts["C_terminal_Median_AA_Entropy"]
    )
    alignment_contrasts["ADD_Lower_Entropy"] = (
        alignment_contrasts["ADD_Median_AA_Entropy"]
        < alignment_contrasts["C_terminal_Median_AA_Entropy"]
    )
    alignment_contrasts.to_csv(
        SUMMARIES / "ADD_vs_C_terminal_entropy_by_alignment.csv", index=False
    )

    consensus = (
        sites.groupby(["Scope", "Human_Position", "Domain", "Structural_Region"])
        .agg(
            Aligners_Retaining_Site=("Aligner", "nunique"),
            Median_AA_Entropy=("Amino_Acid_Shannon_Entropy", "median"),
            Minimum_Major_AA_Frequency=("Major_Amino_Acid_Frequency", "min"),
            Maximum_AA_State_Count=("Amino_Acid_State_Count", "max"),
            Minimum_Occupancy=("Occupancy", "min"),
        )
        .reset_index()
    )
    consensus = consensus.merge(
        concordance[
            [
                "Scope",
                "Human_Position",
                "All_Three_Resolved_N",
                "All_Three_Resolved_Fraction",
                "All_Three_Exact_Among_Resolved",
            ]
        ],
        on=["Scope", "Human_Position"],
        how="left",
    )
    consensus["Alignment_Stable"] = (
        (consensus["All_Three_Resolved_Fraction"] >= 0.80)
        & (consensus["All_Three_Exact_Among_Resolved"] >= 0.90)
    )
    consensus["Invariant_Across_Each_Aligner"] = consensus["Maximum_AA_State_Count"] == 1
    consensus.to_csv(SUMMARIES / "amino_acid_constraint_by_site.csv", index=False)

    eligible = consensus.loc[(consensus["Aligners_Retaining_Site"] == 3) & consensus["Alignment_Stable"]].copy()
    domain_summary = (
        eligible.groupby(["Scope", "Domain"])
        .agg(
            Eligible_Sites=("Human_Position", "size"),
            Median_AA_Entropy=("Median_AA_Entropy", "median"),
            Mean_AA_Entropy=("Median_AA_Entropy", "mean"),
            Median_Major_AA_Frequency=("Minimum_Major_AA_Frequency", "median"),
            Invariant_Sites=("Invariant_Across_Each_Aligner", "sum"),
        )
        .reset_index()
    )
    domain_summary["Invariant_Fraction"] = domain_summary["Invariant_Sites"] / domain_summary["Eligible_Sites"]
    domain_summary.to_csv(SUMMARIES / "amino_acid_constraint_by_domain.csv", index=False)

    omnibus_rows: list[dict[str, object]] = []
    pairwise_rows: list[dict[str, object]] = []
    for scope, data in eligible.groupby("Scope"):
        domain_groups = {
            domain: frame["Median_AA_Entropy"].to_numpy()
            for domain, frame in data.groupby("Domain")
            if len(frame) >= 2
        }
        statistic, p_value = kruskal(*domain_groups.values()) if len(domain_groups) >= 2 else (np.nan, np.nan)
        omnibus_rows.append(
            {
                "Scope": scope,
                "Domains_Tested": ";".join(domain_groups),
                "Kruskal_Statistic": statistic,
                "Kruskal_p": p_value,
            }
        )
        for first, second in combinations(domain_groups, 2):
            stat, pair_p = mannwhitneyu(domain_groups[first], domain_groups[second], alternative="two-sided")
            pairwise_rows.append(
                {
                    "Scope": scope,
                    "Domain_1": first,
                    "Domain_2": second,
                    "Domain_1_N": len(domain_groups[first]),
                    "Domain_2_N": len(domain_groups[second]),
                    "Mann_Whitney_U": stat,
                    "p": pair_p,
                }
            )
    pd.DataFrame(omnibus_rows).to_csv(SUMMARIES / "constraint_domain_omnibus_tests.csv", index=False)
    pairwise = pd.DataFrame(pairwise_rows)
    pairwise["q_within_scope"] = pairwise.groupby("Scope")["p"].transform(bh_adjust)
    pairwise.to_csv(SUMMARIES / "constraint_domain_pairwise_tests.csv", index=False)

    structural_rows: list[dict[str, object]] = []
    for scope, data in eligible.groupby("Scope"):
        cterminal = data.loc[data["Domain"] == "C_terminal_MTase_like"]
        for region in ("switching_helix", "switching_helix_adjacent_loop"):
            focal = cterminal.loc[cterminal["Structural_Region"] == region, "Median_AA_Entropy"].to_numpy()
            background = cterminal.loc[cterminal["Structural_Region"] == "other", "Median_AA_Entropy"].to_numpy()
            if len(focal) and len(background):
                statistic, p_value = mannwhitneyu(focal, background, alternative="two-sided")
            else:
                statistic, p_value = np.nan, np.nan
            structural_rows.append(
                {
                    "Scope": scope,
                    "Structural_Region": region,
                    "Region_Sites": len(focal),
                    "Background_C_terminal_Sites": len(background),
                    "Region_Median_Entropy": np.median(focal) if len(focal) else np.nan,
                    "Background_Median_Entropy": np.median(background) if len(background) else np.nan,
                    "Mann_Whitney_U": statistic,
                    "p": p_value,
                }
            )
    structural = pd.DataFrame(structural_rows)
    structural["q_within_scope"] = structural.groupby("Scope")["p"].transform(bh_adjust)
    structural.to_csv(SUMMARIES / "switching_helix_constraint_tests.csv", index=False)

    print(domain_summary.to_string(index=False))
    print("\nKruskal-Wallis tests")
    print(pd.DataFrame(omnibus_rows).to_string(index=False))


if __name__ == "__main__":
    main()
