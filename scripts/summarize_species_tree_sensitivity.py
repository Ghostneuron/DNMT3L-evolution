from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
V2 = PROJECT / "analysis_v2"
OUT = V2 / "summaries"
JSON_DIR = V2 / "species_tree_hyphy"


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


def parse_dataset(dataset: str, mapping: pd.DataFrame) -> pd.DataFrame:
    site_map = mapping.loc[mapping["Dataset"] == dataset].sort_values("Filtered_Site").copy().reset_index(drop=True)
    result = site_map[["Scope", "Aligner", "Dataset", "Filtered_Site", "Human_Position", "Domain", "Structural_Region"]].copy()

    fubar = json.loads((JSON_DIR / dataset / "FUBAR_species_tree.json").read_text())["MLE"]["content"]["0"]
    fel = json.loads((JSON_DIR / dataset / "FEL_species_tree.json").read_text())["MLE"]["content"]["0"]
    meme = json.loads((JSON_DIR / dataset / "MEME_species_tree.json").read_text())["MLE"]["content"]["0"]
    if not (len(result) == len(fubar) == len(fel) == len(meme)):
        raise ValueError(f"Species-tree row mismatch for {dataset}")

    result["SpeciesTree_FUBAR_alpha"] = [row[0] for row in fubar]
    result["SpeciesTree_FUBAR_beta"] = [row[1] for row in fubar]
    result["SpeciesTree_FUBAR_PP_negative"] = [row[3] for row in fubar]
    result["SpeciesTree_FUBAR_PP_positive"] = [row[4] for row in fubar]
    result["SpeciesTree_FEL_alpha"] = [row[0] for row in fel]
    result["SpeciesTree_FEL_beta"] = [row[1] for row in fel]
    result["SpeciesTree_FEL_p"] = [row[4] for row in fel]
    result["SpeciesTree_FEL_q"] = bh_adjust(result["SpeciesTree_FEL_p"])
    result["SpeciesTree_MEME_p"] = [row[6] for row in meme]
    result["SpeciesTree_MEME_q"] = bh_adjust(result["SpeciesTree_MEME_p"])
    result["SpeciesTree_Joint_Pervasive_q10"] = (
        (result["SpeciesTree_FUBAR_PP_positive"] >= 0.90)
        & (result["SpeciesTree_FEL_beta"] > result["SpeciesTree_FEL_alpha"])
        & (result["SpeciesTree_FEL_q"] <= 0.10)
    )
    result["SpeciesTree_Joint_Pervasive_q05"] = (
        (result["SpeciesTree_FUBAR_PP_positive"] >= 0.90)
        & (result["SpeciesTree_FEL_beta"] > result["SpeciesTree_FEL_alpha"])
        & (result["SpeciesTree_FEL_q"] <= 0.05)
    )
    result["SpeciesTree_MEME_q10"] = result["SpeciesTree_MEME_q"] <= 0.10
    result["SpeciesTree_MEME_q05"] = result["SpeciesTree_MEME_q"] <= 0.05
    return result


def main() -> None:
    mapping = pd.read_csv(OUT / "retained_site_map.csv")
    mapping["Dataset"] = mapping["Scope"] + "_" + mapping["Aligner"]
    datasets = [f"{scope}_MACSE" for scope in ("Placentalia", "Euarchontoglires", "Laurasiatheria", "Sauropsida")]
    species_sites = pd.concat([parse_dataset(dataset, mapping) for dataset in datasets], ignore_index=True)

    gene_sites = pd.read_csv(OUT / "sitewise_selection_all_models.csv")
    gene_sites = gene_sites.loc[gene_sites["Aligner"] == "MACSE"].copy()
    gene_columns = [
        "Scope",
        "Aligner",
        "Dataset",
        "Filtered_Site",
        "Human_Position",
        "FUBAR_PP_positive",
        "FEL_p",
        "FEL_q",
        "MEME_p",
        "MEME_q",
        "Joint_Pervasive_q10",
        "Joint_Pervasive_q05",
        "MEME_q10",
        "MEME_q05",
    ]
    comparison = gene_sites[gene_columns].merge(
        species_sites,
        on=["Scope", "Aligner", "Dataset", "Filtered_Site", "Human_Position"],
        how="inner",
    )
    comparison.to_csv(OUT / "gene_tree_vs_species_tree_sitewise.csv", index=False)
    comparison.loc[comparison["Human_Position"] == 358].to_csv(OUT / "h358_tree_topology_sensitivity.csv", index=False)

    candidate = pd.read_csv(OUT / "diversifying_candidates_and_h358.csv")
    candidate = candidate.merge(
        species_sites,
        on=["Scope", "Human_Position", "Domain", "Structural_Region"],
        how="left",
    )
    candidate.to_csv(OUT / "candidate_species_tree_sensitivity.csv", index=False)

    busted_rows: list[dict[str, object]] = []
    for dataset in datasets:
        payload = json.loads((JSON_DIR / dataset / "BUSTED_species_tree.json").read_text())
        tests = payload.get("test results", {})
        busted_rows.append(
            {
                "Dataset": dataset,
                "Scope": dataset.rsplit("_", 1)[0],
                "BUSTED_SpeciesTree_LRT": tests.get("LRT"),
                "BUSTED_SpeciesTree_p": tests.get("p-value"),
            }
        )
    busted = pd.DataFrame(busted_rows)
    busted["BUSTED_SpeciesTree_q_across_4"] = bh_adjust(busted["BUSTED_SpeciesTree_p"])
    standard = pd.read_csv(OUT / "BUSTED_gene_wide_summary.csv")
    standard = standard.loc[standard["Aligner"] == "MACSE", ["Dataset", "BUSTED_LRT", "BUSTED_p", "BUSTED_q_across_12"]]
    busted = standard.merge(busted, on="Dataset", how="inner")
    busted.to_csv(OUT / "BUSTED_tree_topology_sensitivity.csv", index=False)

    print("BUSTED gene-tree versus species-tree sensitivity")
    print(busted.to_string(index=False))
    print("\nh358")
    print(comparison.loc[comparison["Human_Position"] == 358, ["Scope", "FUBAR_PP_positive", "SpeciesTree_FUBAR_PP_positive", "FEL_q", "SpeciesTree_FEL_q", "MEME_q", "SpeciesTree_MEME_q"]].to_string(index=False))


if __name__ == "__main__":
    main()
