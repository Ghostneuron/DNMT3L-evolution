from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


PROJECT = Path(__file__).resolve().parents[1]
V2 = PROJECT / "analysis_v2"
HYPHY = V2 / "hyphy"
SUMMARIES = V2 / "summaries"
SOURCE = PROJECT / "source_package" / "analysis" / "orthology_rebuild_20260826"


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


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def parse_site_json(dataset: str, method: str, site_map: pd.DataFrame) -> pd.DataFrame:
    path = HYPHY / dataset / f"{method}_all_sites_standard.json"
    data = load_json(path)
    rows = data["MLE"]["content"]["0"]
    mapping = site_map.loc[site_map["Dataset"] == dataset].sort_values("Filtered_Site").copy()
    if len(rows) != len(mapping):
        raise ValueError(f"{dataset} {method}: {len(rows)} JSON rows != {len(mapping)} mapped sites")

    base = mapping[["Scope", "Aligner", "Dataset", "Filtered_Site", "Human_Position", "Domain", "Structural_Region", "Resolved_N", "Occupancy", "Amino_Acid_State_Count", "Major_Amino_Acid_Frequency", "Amino_Acid_Shannon_Entropy", "Amino_Acid_Counts"]].reset_index(drop=True)
    if method == "FUBAR":
        values = pd.DataFrame(
            [row[:6] for row in rows],
            columns=["FUBAR_alpha", "FUBAR_beta", "FUBAR_beta_minus_alpha", "FUBAR_PP_negative", "FUBAR_PP_positive", "FUBAR_BF_positive"],
        )
    elif method == "FEL":
        values = pd.DataFrame(
            [row[:6] for row in rows],
            columns=["FEL_alpha", "FEL_beta", "FEL_neutral_rate", "FEL_LRT", "FEL_p", "FEL_total_branch_length"],
        )
        values["FEL_q"] = bh_adjust(values["FEL_p"])
        values["FEL_direction"] = np.where(
            values["FEL_beta"] > values["FEL_alpha"],
            "diversifying",
            np.where(values["FEL_beta"] < values["FEL_alpha"], "purifying", "neutral"),
        )
    elif method == "MEME":
        values = pd.DataFrame(
            [row[:13] for row in rows],
            columns=[
                "MEME_alpha",
                "MEME_beta_background",
                "MEME_weight_background",
                "MEME_beta_positive",
                "MEME_weight_positive",
                "MEME_LRT",
                "MEME_p",
                "MEME_branches_EBF100",
                "MEME_total_branch_length",
                "MEME_logL",
                "MEME_FEL_logL",
                "MEME_FEL_alpha",
                "MEME_FEL_beta",
            ],
        )
        values["MEME_q"] = bh_adjust(values["MEME_p"])
    else:
        raise ValueError(method)
    return pd.concat([base, values], axis=1)


def parse_busted(dataset: str) -> dict[str, object]:
    path = HYPHY / dataset / "BUSTED_all_sites_standard.json"
    data = load_json(path)
    tests = data.get("test results", {})
    return {
        "Dataset": dataset,
        "Scope": dataset.rsplit("_", 1)[0],
        "Aligner": dataset.rsplit("_", 1)[1],
        "BUSTED_LRT": tests.get("LRT"),
        "BUSTED_p": tests.get("p-value"),
        "BUSTED_background_proportion": tests.get("background proportion"),
    }


def build_masked_comparison(new_sites: pd.DataFrame) -> pd.DataFrame:
    old_site_file = SOURCE / "confidence_filtered_alignments" / "codon_confidence_filter_by_site.csv"
    old_map = pd.read_csv(old_site_file)
    old_map = old_map.loc[
        old_map["Alignment"].str.startswith("Placentalia_") & old_map["Retained"].astype(bool)
    ].copy()
    old_map["Alignment_Site"] = old_map.groupby("Alignment").cumcount() + 1
    old_map = old_map[["Alignment", "Alignment_Site", "Human_Position"]]

    fubar_rows: list[pd.DataFrame] = []
    for aligner in ("MACSE", "MAFFT", "PRANK"):
        label = f"Placentalia_{aligner}"
        path = SOURCE / "selection" / label / "FUBAR_deduplicated.json"
        rows = load_json(path)["MLE"]["content"]["0"]
        frame = old_map.loc[old_map["Alignment"] == label].sort_values("Alignment_Site").copy()
        if len(frame) != len(rows):
            raise ValueError(f"Old masked FUBAR map mismatch for {label}")
        frame["Aligner"] = aligner
        frame["Masked_FUBAR_PP_positive"] = [row[4] for row in rows]
        fubar_rows.append(frame[["Aligner", "Alignment_Site", "Human_Position", "Masked_FUBAR_PP_positive"]])
    old_fubar = pd.concat(fubar_rows, ignore_index=True)

    old_meme = pd.read_csv(SOURCE / "selection_audit" / "MEME_all_sites_with_BH.csv")
    old_meme = old_meme.loc[old_meme["Taxonomic_Scope"] == "Placentalia", ["Aligner", "Alignment_Site", "MEME_p", "MEME_BH_q"]].rename(
        columns={"MEME_p": "Masked_MEME_p", "MEME_BH_q": "Masked_MEME_q"}
    )
    old = old_fubar.merge(old_meme, on=["Aligner", "Alignment_Site"], how="left")

    new = new_sites.loc[
        new_sites["Scope"] == "Placentalia",
        ["Aligner", "Human_Position", "FUBAR_PP_positive", "FEL_p", "FEL_q", "MEME_p", "MEME_q"],
    ].rename(
        columns={
            "FUBAR_PP_positive": "Unmasked_FUBAR_PP_positive",
            "FEL_p": "Unmasked_FEL_p",
            "FEL_q": "Unmasked_FEL_q",
            "MEME_p": "Unmasked_MEME_p",
            "MEME_q": "Unmasked_MEME_q",
        }
    )
    comparison = old.merge(new, on=["Aligner", "Human_Position"], how="outer")
    comparison["FUBAR_PP_Shift_Unmasked_Minus_Masked"] = comparison["Unmasked_FUBAR_PP_positive"] - comparison["Masked_FUBAR_PP_positive"]
    comparison["MEME_q_Shift_Unmasked_Minus_Masked"] = comparison["Unmasked_MEME_q"] - comparison["Masked_MEME_q"]
    return comparison.sort_values(["Human_Position", "Aligner"])


def main() -> None:
    SUMMARIES.mkdir(parents=True, exist_ok=True)
    site_map = pd.read_csv(SUMMARIES / "retained_site_map.csv")
    site_map["Dataset"] = site_map["Scope"] + "_" + site_map["Aligner"]
    datasets = sorted(site_map["Dataset"].unique())
    if len(datasets) != 12:
        raise ValueError(f"Expected 12 datasets, found {len(datasets)}")

    method_tables = {
        method: pd.concat([parse_site_json(dataset, method, site_map) for dataset in datasets], ignore_index=True)
        for method in ("FUBAR", "FEL", "MEME")
    }
    keys = ["Scope", "Aligner", "Dataset", "Filtered_Site", "Human_Position", "Domain", "Structural_Region", "Resolved_N", "Occupancy", "Amino_Acid_State_Count", "Major_Amino_Acid_Frequency", "Amino_Acid_Shannon_Entropy", "Amino_Acid_Counts"]
    sites = method_tables["FUBAR"].merge(method_tables["FEL"], on=keys, how="inner").merge(method_tables["MEME"], on=keys, how="inner")

    concordance = pd.read_csv(SUMMARIES / "cross_aligner_concordance_by_site.csv")
    sites = sites.merge(concordance, on=["Scope", "Human_Position", "Domain", "Structural_Region"], how="left")
    sites["FUBAR_Positive"] = sites["FUBAR_PP_positive"] >= 0.90
    sites["FUBAR_Purifying"] = sites["FUBAR_PP_negative"] >= 0.90
    sites["FEL_Positive_q10"] = (sites["FEL_direction"] == "diversifying") & (sites["FEL_q"] <= 0.10)
    sites["FEL_Positive_q05"] = (sites["FEL_direction"] == "diversifying") & (sites["FEL_q"] <= 0.05)
    sites["FEL_Purifying_q05"] = (sites["FEL_direction"] == "purifying") & (sites["FEL_q"] <= 0.05)
    sites["MEME_q10"] = sites["MEME_q"] <= 0.10
    sites["MEME_q05"] = sites["MEME_q"] <= 0.05
    sites["Joint_Pervasive_q10"] = sites["FUBAR_Positive"] & sites["FEL_Positive_q10"]
    sites["Joint_Pervasive_q05"] = sites["FUBAR_Positive"] & sites["FEL_Positive_q05"]
    sites["Joint_Purifying_q05"] = sites["FUBAR_Purifying"] & sites["FEL_Purifying_q05"]
    sites.to_csv(SUMMARIES / "sitewise_selection_all_models.csv", index=False)

    call_counts = (
        sites.groupby(["Scope", "Aligner", "Dataset"])[
            [
                "FUBAR_Positive",
                "FUBAR_Purifying",
                "FEL_Positive_q10",
                "FEL_Positive_q05",
                "FEL_Purifying_q05",
                "MEME_q10",
                "MEME_q05",
                "Joint_Pervasive_q10",
                "Joint_Pervasive_q05",
                "Joint_Purifying_q05",
            ]
        ]
        .sum()
        .reset_index()
    )
    call_counts.to_csv(SUMMARIES / "selection_call_counts_by_alignment.csv", index=False)

    jaccard_rows: list[dict[str, object]] = []
    for scope, scope_data in sites.groupby("Scope"):
        for criterion in ("FUBAR_Positive", "FEL_Positive_q10", "MEME_q10", "Joint_Pervasive_q10", "Joint_Purifying_q05"):
            sets = {
                aligner: set(scope_data.loc[(scope_data["Aligner"] == aligner) & scope_data[criterion], "Human_Position"].astype(int))
                for aligner in ("MACSE", "MAFFT", "PRANK")
            }
            for first, second in (("MACSE", "MAFFT"), ("MACSE", "PRANK"), ("MAFFT", "PRANK")):
                intersection = sets[first] & sets[second]
                union = sets[first] | sets[second]
                jaccard_rows.append(
                    {
                        "Scope": scope,
                        "Criterion": criterion,
                        "Aligner_1": first,
                        "Aligner_2": second,
                        "Aligner_1_Sites": len(sets[first]),
                        "Aligner_2_Sites": len(sets[second]),
                        "Intersection_Sites": len(intersection),
                        "Union_Sites": len(union),
                        "Jaccard": len(intersection) / len(union) if union else 1.0,
                        "Intersection_Human_Positions": ";".join(str(value) for value in sorted(intersection)),
                    }
                )
    pd.DataFrame(jaccard_rows).to_csv(SUMMARIES / "cross_aligner_call_jaccard.csv", index=False)

    grouped_rows: list[dict[str, object]] = []
    for (scope, position), group in sites.groupby(["Scope", "Human_Position"]):
        concord = group["All_Three_Exact_Among_Resolved"].dropna()
        joint_coverage = group["All_Three_Resolved_Fraction"].dropna()
        stable_fraction = float(concord.iloc[0]) if not concord.empty else np.nan
        joint_coverage_fraction = float(joint_coverage.iloc[0]) if not joint_coverage.empty else np.nan
        stable = bool(
            not np.isnan(stable_fraction)
            and not np.isnan(joint_coverage_fraction)
            and stable_fraction >= 0.90
            and joint_coverage_fraction >= 0.80
        )
        grouped_rows.append(
            {
                "Scope": scope,
                "Human_Position": int(position),
                "Domain": group["Domain"].iloc[0],
                "Structural_Region": group["Structural_Region"].iloc[0],
                "Aligners_Tested": group["Aligner"].nunique(),
                "Alignment_Stable": stable,
                "Three_Aligner_AA_Concordance": stable_fraction,
                "Three_Aligner_Joint_Resolved_Fraction": joint_coverage_fraction,
                "Jointly_Resolved_Taxa": group["All_Three_Resolved_N"].dropna().iloc[0] if group["All_Three_Resolved_N"].notna().any() else np.nan,
                "FUBAR_Positive_Aligners": int(group["FUBAR_Positive"].sum()),
                "FEL_Positive_q10_Aligners": int(group["FEL_Positive_q10"].sum()),
                "FEL_Positive_q05_Aligners": int(group["FEL_Positive_q05"].sum()),
                "Joint_Pervasive_q10_Aligners": int(group["Joint_Pervasive_q10"].sum()),
                "Joint_Pervasive_q05_Aligners": int(group["Joint_Pervasive_q05"].sum()),
                "MEME_q10_Aligners": int(group["MEME_q10"].sum()),
                "MEME_q05_Aligners": int(group["MEME_q05"].sum()),
                "Joint_Purifying_q05_Aligners": int(group["Joint_Purifying_q05"].sum()),
            }
        )
    robustness = pd.DataFrame(grouped_rows)
    robustness["Robust_Pervasive_q10"] = robustness["Alignment_Stable"] & (robustness["Joint_Pervasive_q10_Aligners"] >= 2)
    robustness["Robust_Pervasive_q05"] = robustness["Alignment_Stable"] & (robustness["Joint_Pervasive_q05_Aligners"] >= 2)
    robustness["Robust_Episodic_q10"] = robustness["Alignment_Stable"] & (robustness["MEME_q10_Aligners"] >= 2)
    robustness["Robust_Episodic_q05"] = robustness["Alignment_Stable"] & (robustness["MEME_q05_Aligners"] >= 2)
    robustness["Robust_Purifying_q05"] = robustness["Alignment_Stable"] & (robustness["Joint_Purifying_q05_Aligners"] >= 2)
    robustness.to_csv(SUMMARIES / "site_robustness_by_scope.csv", index=False)

    candidate_mask = robustness[["Robust_Pervasive_q10", "Robust_Episodic_q10"]].any(axis=1) | (robustness["Human_Position"] == 358)
    robustness.loc[candidate_mask].to_csv(SUMMARIES / "diversifying_candidates_and_h358.csv", index=False)

    cross_rows: list[dict[str, object]] = []
    for position, group in robustness.groupby("Human_Position"):
        by_scope = group.set_index("Scope")
        def flag(scope: str, column: str) -> bool:
            return bool(by_scope.loc[scope, column]) if scope in by_scope.index else False
        cross_rows.append(
            {
                "Human_Position": int(position),
                "Domain": group["Domain"].iloc[0],
                "Structural_Region": group["Structural_Region"].iloc[0],
                "Euarchontoglires_Pervasive_q10": flag("Euarchontoglires", "Robust_Pervasive_q10"),
                "Laurasiatheria_Pervasive_q10": flag("Laurasiatheria", "Robust_Pervasive_q10"),
                "Sauropsida_Pervasive_q10": flag("Sauropsida", "Robust_Pervasive_q10"),
                "Euarchontoglires_Episodic_q10": flag("Euarchontoglires", "Robust_Episodic_q10"),
                "Laurasiatheria_Episodic_q10": flag("Laurasiatheria", "Robust_Episodic_q10"),
                "Sauropsida_Episodic_q10": flag("Sauropsida", "Robust_Episodic_q10"),
            }
        )
    cross = pd.DataFrame(cross_rows)
    cross["Cross_Placental_Pervasive_Recurrence"] = cross["Euarchontoglires_Pervasive_q10"] & cross["Laurasiatheria_Pervasive_q10"]
    cross["Cross_Placental_Episodic_Recurrence"] = cross["Euarchontoglires_Episodic_q10"] & cross["Laurasiatheria_Episodic_q10"]
    cross["Deep_Lineage_Pervasive_Recurrence"] = (cross["Euarchontoglires_Pervasive_q10"] | cross["Laurasiatheria_Pervasive_q10"]) & cross["Sauropsida_Pervasive_q10"]
    cross["Deep_Lineage_Episodic_Recurrence"] = (cross["Euarchontoglires_Episodic_q10"] | cross["Laurasiatheria_Episodic_q10"]) & cross["Sauropsida_Episodic_q10"]
    cross.to_csv(SUMMARIES / "cross_scope_recurrence.csv", index=False)

    domain_rows: list[dict[str, object]] = []
    for scope, scope_data in robustness.groupby("Scope"):
        eligible = scope_data.loc[(scope_data["Aligners_Tested"] >= 2) & scope_data["Alignment_Stable"]]
        for category in ("Robust_Pervasive_q10", "Robust_Episodic_q10", "Robust_Purifying_q05"):
            for domain in sorted(eligible["Domain"].unique()):
                in_domain = eligible["Domain"] == domain
                selected = eligible[category].astype(bool)
                table = [
                    [int((in_domain & selected).sum()), int((in_domain & ~selected).sum())],
                    [int((~in_domain & selected).sum()), int((~in_domain & ~selected).sum())],
                ]
                odds_ratio, p_value = fisher_exact(table, alternative="greater")
                domain_rows.append(
                    {
                        "Scope": scope,
                        "Category": category,
                        "Domain": domain,
                        "Eligible_Sites": int(in_domain.sum()),
                        "Selected_Sites": int((in_domain & selected).sum()),
                        "Other_Eligible_Sites": int((~in_domain).sum()),
                        "Other_Selected_Sites": int((~in_domain & selected).sum()),
                        "Fisher_Odds_Ratio": odds_ratio,
                        "Fisher_p": p_value,
                    }
                )
    domains = pd.DataFrame(domain_rows)
    domains["Fisher_q_within_scope_category"] = domains.groupby(["Scope", "Category"])["Fisher_p"].transform(bh_adjust)
    domains.to_csv(SUMMARIES / "domain_enrichment.csv", index=False)

    structural_rows: list[dict[str, object]] = []
    for scope, scope_data in robustness.groupby("Scope"):
        eligible = scope_data.loc[(scope_data["Aligners_Tested"] >= 2) & scope_data["Alignment_Stable"]]
        for category in ("Robust_Pervasive_q10", "Robust_Episodic_q10", "Robust_Purifying_q05"):
            for region in ("switching_helix", "switching_helix_adjacent_loop"):
                in_region = eligible["Structural_Region"] == region
                selected = eligible[category].astype(bool)
                table = [
                    [int((in_region & selected).sum()), int((in_region & ~selected).sum())],
                    [int((~in_region & selected).sum()), int((~in_region & ~selected).sum())],
                ]
                odds_ratio, p_value = fisher_exact(table, alternative="greater")
                structural_rows.append(
                    {
                        "Scope": scope,
                        "Category": category,
                        "Structural_Region": region,
                        "Eligible_Sites": int(in_region.sum()),
                        "Selected_Sites": int((in_region & selected).sum()),
                        "Other_Eligible_Sites": int((~in_region).sum()),
                        "Other_Selected_Sites": int((~in_region & selected).sum()),
                        "Fisher_Odds_Ratio": odds_ratio,
                        "Fisher_p": p_value,
                    }
                )
    structural = pd.DataFrame(structural_rows)
    structural["Fisher_q_within_scope_category"] = structural.groupby(["Scope", "Category"])["Fisher_p"].transform(bh_adjust)
    structural.to_csv(SUMMARIES / "switching_helix_enrichment.csv", index=False)

    busted = pd.DataFrame([parse_busted(dataset) for dataset in datasets])
    busted["BUSTED_q_across_12"] = bh_adjust(busted["BUSTED_p"])
    busted["BUSTED_q05"] = busted["BUSTED_q_across_12"] <= 0.05
    busted.to_csv(SUMMARIES / "BUSTED_gene_wide_summary.csv", index=False)
    busted_scope = busted.groupby("Scope").agg(
        Aligners_Tested=("Aligner", "nunique"),
        Significant_Aligners=("BUSTED_q05", "sum"),
        Minimum_p=("BUSTED_p", "min"),
        Maximum_p=("BUSTED_p", "max"),
        Minimum_q=("BUSTED_q_across_12", "min"),
        Maximum_q=("BUSTED_q_across_12", "max"),
    ).reset_index()
    busted_scope["Robust_Gene_Wide_Selection"] = busted_scope["Significant_Aligners"] >= 2
    busted_scope.to_csv(SUMMARIES / "BUSTED_robustness_by_scope.csv", index=False)

    comparison = build_masked_comparison(sites)
    comparison.to_csv(SUMMARIES / "masked_vs_unmasked_placental_comparison.csv", index=False)
    sites.loc[sites["Human_Position"] == 358].to_csv(SUMMARIES / "h358_expanded_audit.csv", index=False)

    print("Robust candidate counts by scope")
    print(
        robustness.groupby("Scope")[["Robust_Pervasive_q10", "Robust_Pervasive_q05", "Robust_Episodic_q10", "Robust_Episodic_q05", "Robust_Purifying_q05"]]
        .sum()
        .to_string()
    )
    print("\nCross-scope recurrence")
    print(cross[["Cross_Placental_Pervasive_Recurrence", "Cross_Placental_Episodic_Recurrence", "Deep_Lineage_Pervasive_Recurrence", "Deep_Lineage_Episodic_Recurrence"]].sum().to_string())
    print("\nBUSTED")
    print(busted[["Dataset", "BUSTED_p", "BUSTED_q_across_12"]].to_string(index=False))


if __name__ == "__main__":
    main()
