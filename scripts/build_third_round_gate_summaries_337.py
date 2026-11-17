#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
ATLAS = PROJECT / "results/functional_site_states_by_record.tsv"
EXTENDED = PROJECT / "analysis/regulatory_coevolution_ptm_20261104/results"
G9A_OUTPUT = PROJECT / "results/g9a_triplet_pass_summary.tsv"
G9A_MUROID_OUTPUT = PROJECT / "results/g9a_muroid_clade_summary.tsv"
INTERFACE_OUTPUT = PROJECT / "results/matched_interface_panel_audit.tsv"

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


def build_g9a_summary() -> None:
    data = pd.read_csv(ATLAS, sep="\t")
    data = data.loc[data["feature_id"].str.startswith("G9A_INTERFACE")].copy()
    pivot = data.pivot(
        index=["safe_id", "species", "broad_lineage"],
        columns="human_np_position",
        values=["query_residue", "mapping_status"],
    )
    required = [105, 106, 107]
    pass_all = pivot["mapping_status"][required].eq("PASS").all(axis=1)
    states = pivot.loc[pass_all, "query_residue"][required].astype(str).agg("".join, axis=1)
    passed = pivot.loc[pass_all].index.to_frame(index=False)
    passed["triplet_state"] = states.to_numpy()

    totals = pivot.index.to_frame(index=False).groupby("broad_lineage").size().to_dict()
    pass_counts = passed.groupby("broad_lineage").size().to_dict()
    rows: list[dict[str, object]] = []
    for (lineage, state), count in passed.groupby(["broad_lineage", "triplet_state"]).size().items():
        rows.append(
            {
                "broad_lineage": lineage,
                "lineage_records": int(totals[lineage]),
                "pass_at_all_three_positions": int(pass_counts[lineage]),
                "triplet_state": state,
                "state_count_among_pass": int(count),
            }
        )
    rows.append(
        {
            "broad_lineage": "ALL",
            "lineage_records": len(pivot),
            "pass_at_all_three_positions": int(pass_all.sum()),
            "triplet_state": "ALL_PASS_RECORDS",
            "state_count_among_pass": int(pass_all.sum()),
        }
    )
    pd.DataFrame(rows).to_csv(G9A_OUTPUT, sep="\t", index=False)

    eutherians = passed.loc[passed["broad_lineage"].eq("Placentalia")].copy()
    muroids = MURIDAE | CRICETIDAE | SPALACIDAE
    observed_muroids = set(eutherians.loc[eutherians["species"].isin(muroids), "species"])
    if observed_muroids != muroids:
        raise ValueError(f"Unexpected sampled muroid set: {sorted(observed_muroids ^ muroids)}")
    if not eutherians.loc[eutherians["species"].isin(muroids), "triplet_state"].eq("LFI").all():
        raise ValueError("Not every sampled muroid carries LFI")
    nonmuroid = eutherians.loc[~eutherians["species"].isin(muroids)]
    if nonmuroid["triplet_state"].eq("LFI").any():
        raise ValueError("An LFI state was found outside sampled muroids among eutherians")

    family_sets = {
        "Muridae": MURIDAE,
        "Cricetidae": CRICETIDAE,
        "Spalacidae": SPALACIDAE,
    }
    muroid_rows: list[dict[str, object]] = []
    for family, members in family_sets.items():
        subset = eutherians.loc[eutherians["species"].isin(members)]
        muroid_rows.append(
            {
                "scope": "sampled_muroid_family",
                "taxon": family,
                "records_with_triplet": len(subset),
                "LFI": int(subset["triplet_state"].eq("LFI").sum()),
                "other_state": int(subset["triplet_state"].ne("LFI").sum()),
                "interpretation": "LFI fixed among sampled family members",
            }
        )
    muroid_rows.extend(
        [
            {
                "scope": "sampled_muroids",
                "taxon": "Muroidea",
                "records_with_triplet": len(muroids),
                "LFI": len(muroids),
                "other_state": 0,
                "interpretation": "all sampled muroids carry LFI",
            },
            {
                "scope": "sampled_nonmuroid_eutherians",
                "taxon": "Eutheria excluding sampled muroids",
                "records_with_triplet": len(nonmuroid),
                "LFI": int(nonmuroid["triplet_state"].eq("LFI").sum()),
                "other_state": int(nonmuroid["triplet_state"].ne("LFI").sum()),
                "interpretation": "no sampled nonmuroid eutherian carries LFI",
            },
            {
                "scope": "nearest_sampled_nonmuroid",
                "taxon": "Jaculus jaculus",
                "records_with_triplet": 1,
                "LFI": 0,
                "other_state": 1,
                "interpretation": "LLI",
            },
        ]
    )
    pd.DataFrame(muroid_rows).to_csv(G9A_MUROID_OUTPUT, sep="\t", index=False)


def build_interface_audit() -> None:
    original = pd.read_csv(EXTENDED / "dnmt3a_dnmt3l_matched_interface_states.tsv", sep="\t")
    expanded = pd.read_csv(EXTENDED / "dnmt3a_dnmt3l_matched_interface_states_vertebrate.tsv", sep="\t")
    therian_clades = {"Laurasiatheria", "Euarchontoglires", "Marsupialia", "Afrotheria", "Xenarthra"}
    expanded_therians = expanded.loc[expanded["clade"].isin(therian_clades)].copy()
    added = sorted(set(expanded_therians["safe_id"]) - set(original["safe_id"]))
    rows = [
        {
            "panel": "original_therian",
            "records": len(original),
            "placentals": int((original["clade"] != "Marsupialia").sum()),
            "marsupials": int((original["clade"] == "Marsupialia").sum()),
            "np261_F": int(original["DNMT3L_NP_261_state"].eq("F").sum()),
            "added_therian_safe_ids": "",
        },
        {
            "panel": "expanded_therian_subset",
            "records": len(expanded_therians),
            "placentals": int((expanded_therians["clade"] != "Marsupialia").sum()),
            "marsupials": int((expanded_therians["clade"] == "Marsupialia").sum()),
            "np261_F": int(expanded_therians["DNMT3L_NP_261_state"].eq("F").sum()),
            "added_therian_safe_ids": ";".join(added),
        },
        {
            "panel": "expanded_vertebrate",
            "records": len(expanded),
            "placentals": int(expanded["clade"].isin(therian_clades - {"Marsupialia"}).sum()),
            "marsupials": int(expanded["clade"].eq("Marsupialia").sum()),
            "np261_F": int(expanded["DNMT3L_NP_261_state"].eq("F").sum()),
            "added_therian_safe_ids": ";".join(added),
        },
    ]
    pd.DataFrame(rows).to_csv(INTERFACE_OUTPUT, sep="\t", index=False)


def main() -> None:
    build_g9a_summary()
    build_interface_audit()
    print(G9A_OUTPUT)
    print(G9A_MUROID_OUTPUT)
    print(INTERFACE_OUTPUT)


if __name__ == "__main__":
    main()
