from __future__ import annotations

import os

import json
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
V2 = PROJECT / "analysis_v2"
SUMMARIES = V2 / "summaries"
INPUTS = V2 / "selection_inputs"
TREES = V2 / "phylogenies"
OUT = V2 / "conditional_sensitivity"
LOGS = V2 / "logs" / "conditional_sensitivity"


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


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


def valid_json(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        json.loads(path.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError):
        return False
    return True


def run(command: list[str], log: Path, output: Path) -> float:
    if valid_json(output):
        return 0.0
    output.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    with log.open("w") as handle:
        handle.write("COMMAND: " + " ".join(command) + "\n\n")
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, text=True, check=False)
    if completed.returncode:
        raise subprocess.CalledProcessError(completed.returncode, command)
    if not valid_json(output):
        raise ValueError(f"Invalid HyPhy JSON: {output}")
    return time.monotonic() - start


def busted_error_sink(scope: str) -> dict[str, object]:
    dataset = f"{scope}_MACSE"
    output = OUT / dataset / "BUSTED_error_sink.json"
    command = [
        os.environ.get("HYPHY_BIN", "hyphy"),
        "busted",
        "--alignment",
        str(INPUTS / f"{dataset}_selection_input.fasta"),
        "--tree",
        str(TREES / dataset / f"{dataset}.treefile"),
        "--code",
        "Universal",
        "--branches",
        "All",
        "--multiple-hits",
        "None",
        "--error-sink",
        "Yes",
        "--output",
        str(output),
    ]
    elapsed = run(command, LOGS / f"{dataset}_busted_error_sink.log", output)
    tests = json.loads(output.read_text()).get("test results", {})
    return {
        "Scope": scope,
        "Dataset": dataset,
        "BUSTED_ErrorSink_LRT": tests.get("LRT"),
        "BUSTED_ErrorSink_p": tests.get("p-value"),
        "Elapsed_Seconds": round(elapsed, 3),
        "Output": str(output.relative_to(PROJECT)),
        "Command": " ".join(command),
    }


def run_busted_multihit(scope: str) -> dict[str, object]:
    dataset = f"{scope}_MACSE"
    output = OUT / dataset / "BUSTED_DoubleTriple.json"
    command = [
        os.environ.get("HYPHY_BIN", "hyphy"),
        "busted",
        "--alignment",
        str(INPUTS / f"{dataset}_selection_input.fasta"),
        "--tree",
        str(TREES / dataset / f"{dataset}.treefile"),
        "--code",
        "Universal",
        "--branches",
        "All",
        "--srv",
        "Yes",
        "--syn-rates",
        "3",
        "--multiple-hits",
        "Double+Triple",
        "--error-sink",
        "No",
        "--output",
        str(output),
    ]
    elapsed = run(command, LOGS / f"{dataset}_busted_multihit.log", output)
    tests = json.loads(output.read_text()).get("test results", {})
    return {
        "Scope": scope,
        "Dataset": dataset,
        "BUSTED_DoubleTriple_LRT": tests.get("LRT"),
        "BUSTED_DoubleTriple_p": tests.get("p-value"),
        "Elapsed_Seconds": round(elapsed, 3),
        "Output": str(output.relative_to(PROJECT)),
        "Command": " ".join(command),
    }


def site_filter_string(sites: list[int], total_sites: int) -> str:
    ordered = sorted({site for site in sites if 1 <= site <= total_sites})
    if not ordered:
        raise ValueError("No valid filtered sites were supplied")
    # HyPhy 2.5.94 requires comma-prefixed integers here. Its helper also
    # matches decimal prefixes, which are audited and excluded downstream.
    return "," + ",".join(str(site) for site in ordered)


def audit_prefix_expansion(
    targets: list[int], evaluated: list[int], dataset: str, model: str
) -> tuple[str, str, str]:
    target_set = set(targets)
    evaluated_set = set(evaluated)
    missing = sorted(target_set - evaluated_set)
    extras = sorted(evaluated_set - target_set)
    unrelated = [
        site
        for site in extras
        if not any(str(target).startswith(str(site)) for target in target_set)
    ]
    if unrelated:
        raise ValueError(
            f"{dataset}: {model} evaluated unrelated sites {unrelated}; "
            "only decimal-prefix expansion is expected from HyPhy 2.5.94"
        )
    return (
        ";".join(map(str, sorted(evaluated_set))),
        ";".join(map(str, extras)),
        ";".join(map(str, missing)),
    )


def meme_multihit(scope: str, aligner: str, human_positions: list[int], site_map: pd.DataFrame) -> list[dict[str, object]]:
    dataset = f"{scope}_{aligner}"
    mapping = site_map.loc[site_map["Dataset"] == dataset].set_index("Human_Position")
    targets = [(position, int(mapping.loc[position, "Filtered_Site"])) for position in human_positions if position in mapping.index]
    if not targets:
        return []
    total_sites = int(mapping["Filtered_Site"].max())
    site_filter = site_filter_string([site for _, site in targets], total_sites)
    output = OUT / dataset / "MEME_episodic_sites_DoubleTriple.json"
    command = [
        os.environ.get("HYPHY_BIN", "hyphy"),
        "meme",
        "--alignment",
        str(INPUTS / f"{dataset}_selection_input.fasta"),
        "--tree",
        str(TREES / dataset / f"{dataset}.treefile"),
        "--code",
        "Universal",
        "--branches",
        "All",
        "--pvalue",
        "0.1",
        "--multiple-hits",
        "Double+Triple",
        "--site-multihit",
        "Estimate",
        "--limit-to-sites",
        site_filter,
        "--output",
        str(output),
    ]
    elapsed = run(command, LOGS / f"{dataset}_meme_episodic_multihit.log", output)
    payload = json.loads(output.read_text())
    rows = payload["MLE"]["content"]["0"]
    headers = [str(header[0]) for header in payload["MLE"]["headers"]]
    recorded_filter = (
        payload.get("analysis", {})
        .get("settings", {})
        .get("site-filter", {})
        .get("site-filter")
    )
    if recorded_filter != site_filter:
        raise ValueError(f"{dataset}: HyPhy recorded site filter {recorded_filter!r}, expected {site_filter!r}")
    if len(rows) != total_sites:
        raise ValueError(f"{dataset}: {len(rows)} MEME rows != {total_sites} filtered sites")
    p_index = headers.index("p-value")
    logl_indices = (headers.index("MEME LogL"), headers.index("FEL LogL"))
    two_hit_index = headers.index("Relative rate estimate for 2-nucleotide substitutions")
    three_hit_index = headers.index("Relative rate estimate for 3-nucleotide substitutions")
    evaluated_sites = [
        index + 1
        for index, row in enumerate(rows)
        if any(np.isfinite(row[column]) and row[column] != 0 for column in logl_indices)
    ]
    evaluated_text, prefix_text, unevaluable_text = audit_prefix_expansion(
        [site for _, site in targets], evaluated_sites, dataset, "MEME"
    )
    evaluated_set = set(evaluated_sites)
    result: list[dict[str, object]] = []
    for human_position, site in targets:
        row = rows[site - 1]
        if len(row) != len(headers):
            raise ValueError(f"{dataset} site {site}: {len(row)} values != {len(headers)} headers")
        target_evaluable = site in evaluated_set
        if not target_evaluable and not (
            float(row[p_index]) == 1
            and all(np.isfinite(row[index]) and row[index] == 0 for index in logl_indices)
        ):
            raise ValueError(f"{dataset} site {site}: unexpected unevaluable MEME row")
        result.append(
            {
                "Scope": scope,
                "Aligner": aligner,
                "Dataset": dataset,
                "Human_Position": human_position,
                "Filtered_Site": site,
                "MEME_DoubleTriple_p": row[p_index],
                "Two_Nucleotide_Relative_Rate": row[two_hit_index],
                "Three_Nucleotide_Relative_Rate": row[three_hit_index],
                "Site_Filter": site_filter,
                "HyPhy_Evaluated_Sites": evaluated_text,
                "Prefix_Expanded_Sites": prefix_text,
                "Unevaluable_Target_Sites": unevaluable_text,
                "Target_Evaluable": target_evaluable,
                "Elapsed_Seconds": round(elapsed, 3),
                "Output": str(output.relative_to(PROJECT)),
                "Command": " ".join(command),
            }
        )
    return result


def fel_multihit(scope: str, aligner: str, human_positions: list[int], site_map: pd.DataFrame) -> list[dict[str, object]]:
    dataset = f"{scope}_{aligner}"
    mapping = site_map.loc[site_map["Dataset"] == dataset].set_index("Human_Position")
    targets = [(position, int(mapping.loc[position, "Filtered_Site"])) for position in human_positions if position in mapping.index]
    if not targets:
        return []
    total_sites = int(mapping["Filtered_Site"].max())
    site_filter = site_filter_string([site for _, site in targets], total_sites)
    output = OUT / dataset / "FEL_primary_sites_DoubleTriple.json"
    command = [
        os.environ.get("HYPHY_BIN", "hyphy"),
        "fel",
        "--alignment",
        str(INPUTS / f"{dataset}_selection_input.fasta"),
        "--tree",
        str(TREES / dataset / f"{dataset}.treefile"),
        "--code",
        "Universal",
        "--branches",
        "All",
        "--pvalue",
        "0.1",
        "--multiple-hits",
        "Double+Triple",
        "--site-multihit",
        "Estimate",
        "--limit-to-sites",
        site_filter,
        "--output",
        str(output),
    ]
    elapsed = run(command, LOGS / f"{dataset}_fel_multihit.log", output)
    payload = json.loads(output.read_text())
    rows = payload["MLE"]["content"]["0"]
    headers = [str(header[0]) for header in payload["MLE"]["headers"]]
    recorded_filter = (
        payload.get("analysis", {})
        .get("settings", {})
        .get("site-filter", {})
        .get("site-filter")
    )
    if recorded_filter != site_filter:
        raise ValueError(f"{dataset}: HyPhy recorded site filter {recorded_filter!r}, expected {site_filter!r}")
    if len(rows) != total_sites:
        raise ValueError(f"{dataset}: {len(rows)} FEL rows != {total_sites} filtered sites")
    alpha_index = headers.index("alpha")
    beta_index = headers.index("beta")
    p_index = headers.index("p-value")
    branch_length_index = headers.index("Total branch length")
    two_hit_index = headers.index("2H rate")
    three_hit_index = headers.index("3H rate")
    evaluated_sites = [
        index + 1
        for index, row in enumerate(rows)
        if np.isfinite(row[branch_length_index]) and row[branch_length_index] > 0
    ]
    evaluated_text, prefix_text, unevaluable_text = audit_prefix_expansion(
        [site for _, site in targets], evaluated_sites, dataset, "FEL"
    )
    evaluated_set = set(evaluated_sites)
    result: list[dict[str, object]] = []
    for human_position, site in targets:
        row = rows[site - 1]
        if len(row) != len(headers):
            raise ValueError(f"{dataset} site {site}: {len(row)} values != {len(headers)} headers")
        target_evaluable = site in evaluated_set
        if not target_evaluable and not (
            float(row[p_index]) == 1
            and np.isfinite(row[branch_length_index])
            and row[branch_length_index] == 0
        ):
            raise ValueError(f"{dataset} site {site}: unexpected unevaluable FEL row")
        result.append(
            {
                "Scope": scope,
                "Aligner": aligner,
                "Dataset": dataset,
                "Human_Position": human_position,
                "Filtered_Site": site,
                "FEL_DoubleTriple_alpha": row[alpha_index],
                "FEL_DoubleTriple_beta": row[beta_index],
                "FEL_DoubleTriple_p": row[p_index],
                "FEL_DoubleTriple_direction": "diversifying" if row[beta_index] > row[alpha_index] else "purifying_or_neutral",
                "Two_Nucleotide_Relative_Rate": row[two_hit_index],
                "Three_Nucleotide_Relative_Rate": row[three_hit_index],
                "Site_Filter": site_filter,
                "HyPhy_Evaluated_Sites": evaluated_text,
                "Prefix_Expanded_Sites": prefix_text,
                "Unevaluable_Target_Sites": unevaluable_text,
                "Target_Evaluable": target_evaluable,
                "Elapsed_Seconds": round(elapsed, 3),
                "Output": str(output.relative_to(PROJECT)),
                "Command": " ".join(command),
            }
        )
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    scope_summary = pd.read_csv(SUMMARIES / "BUSTED_robustness_by_scope.csv")
    robust_sites = pd.read_csv(SUMMARIES / "site_robustness_by_scope.csv")
    site_map = pd.read_csv(SUMMARIES / "retained_site_map.csv")
    site_map["Dataset"] = site_map["Scope"] + "_" + site_map["Aligner"]

    busted_summary_path = SUMMARIES / "BUSTED_error_sink_summary.csv"
    prior_busted_elapsed: dict[str, float] = {}
    if busted_summary_path.exists() and busted_summary_path.stat().st_size:
        try:
            prior_busted = pd.read_csv(busted_summary_path)
            prior_busted_elapsed = {
                str(row.Scope): float(row.Elapsed_Seconds)
                for row in prior_busted.itertuples(index=False)
                if pd.notna(row.Elapsed_Seconds) and float(row.Elapsed_Seconds) > 0
            }
        except (pd.errors.EmptyDataError, AttributeError, TypeError, ValueError):
            prior_busted_elapsed = {}

    busted_scopes = scope_summary.loc[as_bool(scope_summary["Robust_Gene_Wide_Selection"]), "Scope"].tolist()
    busted_rows: list[dict[str, object]] = []
    if busted_scopes:
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = {executor.submit(busted_error_sink, scope): scope for scope in busted_scopes}
            for future in as_completed(futures):
                row = future.result()
                busted_rows.append(row)
                print(f"BUSTED error-sink {row['Scope']}: {row['Elapsed_Seconds']} s", flush=True)
    busted_columns = [
        "Scope",
        "Dataset",
        "BUSTED_ErrorSink_LRT",
        "BUSTED_ErrorSink_p",
        "Elapsed_Seconds",
        "Output",
        "Command",
    ]
    busted = pd.DataFrame(busted_rows, columns=busted_columns)
    if not busted.empty:
        busted["BUSTED_ErrorSink_q"] = bh_adjust(busted["BUSTED_ErrorSink_p"])
        standard = pd.read_csv(SUMMARIES / "BUSTED_gene_wide_summary.csv")
        standard = standard.loc[standard["Aligner"] == "MACSE", ["Dataset", "BUSTED_p", "BUSTED_q_across_12"]]
        busted = standard.merge(busted, on="Dataset", how="right")
        reused = busted["Elapsed_Seconds"] == 0
        busted.loc[reused, "Elapsed_Seconds"] = busted.loc[reused, "Scope"].map(prior_busted_elapsed).fillna(0)
        busted = busted.sort_values("Scope")
    else:
        busted["BUSTED_ErrorSink_q"] = pd.Series(dtype=float)
    busted.to_csv(busted_summary_path, index=False)

    saturation = pd.read_csv(V2 / "saturation_diagnostics" / "codon_saturation_summary.csv")
    moderate_scopes = saturation.loc[
        (saturation["Aligner"] == "MACSE")
        & saturation["Diagnostic_Result"].astype(str).str.startswith("moderate"),
        "Scope",
    ].astype(str).tolist()
    multihit_busted_scopes = sorted(set(busted_scopes) & set(moderate_scopes))
    busted_multihit_summary_path = SUMMARIES / "BUSTED_multihit_moderate_risk_summary.csv"
    prior_multihit_elapsed: dict[str, float] = {}
    if busted_multihit_summary_path.exists() and busted_multihit_summary_path.stat().st_size:
        try:
            prior_multihit = pd.read_csv(busted_multihit_summary_path)
            prior_multihit_elapsed = {
                str(row.Scope): float(row.Elapsed_Seconds)
                for row in prior_multihit.itertuples(index=False)
                if pd.notna(row.Elapsed_Seconds) and float(row.Elapsed_Seconds) > 0
            }
        except (pd.errors.EmptyDataError, AttributeError, TypeError, ValueError):
            prior_multihit_elapsed = {}
    busted_multihit_rows: list[dict[str, object]] = []
    if multihit_busted_scopes:
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = {
                executor.submit(run_busted_multihit, scope): scope
                for scope in multihit_busted_scopes
            }
            for future in as_completed(futures):
                row = future.result()
                busted_multihit_rows.append(row)
                print(f"BUSTED Double+Triple {row['Scope']}: {row['Elapsed_Seconds']} s", flush=True)
    busted_multihit_columns = [
        "Scope",
        "Dataset",
        "BUSTED_DoubleTriple_LRT",
        "BUSTED_DoubleTriple_p",
        "Elapsed_Seconds",
        "Output",
        "Command",
    ]
    busted_multihit = pd.DataFrame(busted_multihit_rows, columns=busted_multihit_columns)
    if not busted_multihit.empty:
        busted_multihit["BUSTED_DoubleTriple_q_across_tested_scopes"] = bh_adjust(
            busted_multihit["BUSTED_DoubleTriple_p"]
        )
        reused = busted_multihit["Elapsed_Seconds"] == 0
        busted_multihit.loc[reused, "Elapsed_Seconds"] = (
            busted_multihit.loc[reused, "Scope"].map(prior_multihit_elapsed).fillna(0)
        )
        busted_multihit = busted_multihit.sort_values("Scope")
    else:
        busted_multihit["BUSTED_DoubleTriple_q_across_tested_scopes"] = pd.Series(dtype=float)
    busted_multihit.to_csv(busted_multihit_summary_path, index=False)

    strict = robust_sites.loc[as_bool(robust_sites["Robust_Episodic_q05"])]
    target_by_scope: dict[str, list[int]] = {}
    for scope in ("Placentalia", "Euarchontoglires", "Laurasiatheria", "Sauropsida"):
        positions = sorted(strict.loc[strict["Scope"] == scope, "Human_Position"].astype(int).unique().tolist())
        if 358 not in positions:
            positions.append(358)
        target_by_scope[scope] = sorted(positions)

    multihit_rows: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = {
            executor.submit(meme_multihit, scope, aligner, positions, site_map): (scope, aligner)
            for scope, positions in target_by_scope.items()
            for aligner in ("MACSE", "MAFFT", "PRANK")
        }
        for future in as_completed(futures):
            rows = future.result()
            multihit_rows.extend(rows)
            scope, aligner = futures[future]
            print(f"MEME Double+Triple {scope}_{aligner}: {len(rows)} focal sites", flush=True)
    multihit_columns = [
        "Scope",
        "Aligner",
        "Dataset",
        "Human_Position",
        "Filtered_Site",
        "MEME_DoubleTriple_p",
        "Two_Nucleotide_Relative_Rate",
        "Three_Nucleotide_Relative_Rate",
        "Site_Filter",
        "HyPhy_Evaluated_Sites",
        "Prefix_Expanded_Sites",
        "Unevaluable_Target_Sites",
        "Target_Evaluable",
        "Elapsed_Seconds",
        "Output",
        "Command",
    ]
    multihit = pd.DataFrame(multihit_rows, columns=multihit_columns)
    if not multihit.empty:
        multihit["MEME_DoubleTriple_q_across_focal_tests"] = bh_adjust(multihit["MEME_DoubleTriple_p"])
        multihit["Conditional_Significant_q05"] = (
            multihit["MEME_DoubleTriple_q_across_focal_tests"] <= 0.05
        )
        multihit["Conditional_Significant_Aligners"] = multihit.groupby(
            ["Scope", "Human_Position"]
        )["Conditional_Significant_q05"].transform("sum").astype(int)
        multihit["Conditional_Robust_2of3"] = multihit["Conditional_Significant_Aligners"] >= 2
        multihit = multihit.sort_values(["Scope", "Human_Position", "Aligner"])
    else:
        multihit["MEME_DoubleTriple_q_across_focal_tests"] = pd.Series(dtype=float)
        multihit["Conditional_Significant_q05"] = pd.Series(dtype=bool)
        multihit["Conditional_Significant_Aligners"] = pd.Series(dtype=int)
        multihit["Conditional_Robust_2of3"] = pd.Series(dtype=bool)
    multihit.to_csv(SUMMARIES / "MEME_multihit_primary_sites_summary.csv", index=False)

    strict_pervasive = robust_sites.loc[as_bool(robust_sites["Robust_Pervasive_q05"])]
    fel_target_by_scope: dict[str, list[int]] = {}
    for scope in ("Placentalia", "Euarchontoglires", "Laurasiatheria", "Sauropsida"):
        positions = sorted(
            strict_pervasive.loc[strict_pervasive["Scope"] == scope, "Human_Position"]
            .astype(int)
            .unique()
            .tolist()
        )
        if 358 not in positions:
            positions.append(358)
        fel_target_by_scope[scope] = sorted(positions)

    fel_rows: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = {
            executor.submit(fel_multihit, scope, aligner, positions, site_map): (scope, aligner)
            for scope, positions in fel_target_by_scope.items()
            for aligner in ("MACSE", "MAFFT", "PRANK")
        }
        for future in as_completed(futures):
            rows = future.result()
            fel_rows.extend(rows)
            scope, aligner = futures[future]
            print(f"FEL Double+Triple {scope}_{aligner}: {len(rows)} focal sites", flush=True)
    fel_columns = [
        "Scope",
        "Aligner",
        "Dataset",
        "Human_Position",
        "Filtered_Site",
        "FEL_DoubleTriple_alpha",
        "FEL_DoubleTriple_beta",
        "FEL_DoubleTriple_p",
        "FEL_DoubleTriple_direction",
        "Two_Nucleotide_Relative_Rate",
        "Three_Nucleotide_Relative_Rate",
        "Site_Filter",
        "HyPhy_Evaluated_Sites",
        "Prefix_Expanded_Sites",
        "Unevaluable_Target_Sites",
        "Target_Evaluable",
        "Elapsed_Seconds",
        "Output",
        "Command",
    ]
    fel = pd.DataFrame(fel_rows, columns=fel_columns)
    if not fel.empty:
        fel["FEL_DoubleTriple_q_across_focal_tests"] = bh_adjust(fel["FEL_DoubleTriple_p"])
        fel["Conditional_Diversifying_q05"] = (
            (fel["FEL_DoubleTriple_direction"] == "diversifying")
            & (fel["FEL_DoubleTriple_q_across_focal_tests"] <= 0.05)
        )
        fel["Conditional_Significant_Aligners"] = fel.groupby(
            ["Scope", "Human_Position"]
        )["Conditional_Diversifying_q05"].transform("sum").astype(int)
        fel["Conditional_Robust_2of3"] = fel["Conditional_Significant_Aligners"] >= 2
        fel = fel.sort_values(["Scope", "Human_Position", "Aligner"])
    else:
        fel["FEL_DoubleTriple_q_across_focal_tests"] = pd.Series(dtype=float)
        fel["Conditional_Diversifying_q05"] = pd.Series(dtype=bool)
        fel["Conditional_Significant_Aligners"] = pd.Series(dtype=int)
        fel["Conditional_Robust_2of3"] = pd.Series(dtype=bool)
    fel.to_csv(SUMMARIES / "FEL_multihit_primary_sites_summary.csv", index=False)
    print("Conditional sensitivity analyses complete")


if __name__ == "__main__":
    main()
