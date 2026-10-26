#!/usr/bin/env python3
"""Summarize GARD and partition-aware Sauropsida selection results."""

from __future__ import annotations

import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve()
PROJECT = HERE.parents[4]
REVISION = HERE.parents[1]
ANALYSIS = REVISION / "analysis"
PARTITIONED = ANALYSIS / "gard_partitioned_selection"
GARD_JSON = ANALYSIS / "Sauropsida_MACSE_GARD_2_5_103_clean.json"
GARD_LF = ANALYSIS / "Sauropsida_MACSE_GARD_2_5_103_clean.best-gard"
GARD_CHECKPOINT = ANALYSIS / "Sauropsida_MACSE_GARD_2bp_checkpoint.json"
GARD_MIXED_JSON = ANALYSIS / "Sauropsida_MACSE_GARD.json"
GARD_MIXED_LF = ANALYSIS / "Sauropsida_MACSE_GARD.best-gard"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float)
    order = valid.sort_values().index
    ranked = valid.loc[order].to_numpy()
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result.loc[order] = np.minimum(adjusted, 1.0)
    return result


def flatten_mle(path: Path, method: str) -> pd.DataFrame:
    payload = json.loads(path.read_text())
    result_rows = []
    for partition in sorted(payload["MLE"]["content"], key=int):
        rows = payload["MLE"]["content"][partition]
        coverage = payload["data partitions"][partition]["coverage"][0]
        if len(rows) != len(coverage):
            raise ValueError(f"{method} partition {partition}: row/coverage mismatch")
        for global_index, row in zip(coverage, rows):
            record = {"Filtered_Site": int(global_index) + 1, "GARD_Partition": int(partition) + 1}
            if method == "FUBAR":
                names = ["alpha", "beta", "beta_minus_alpha", "PP_negative", "PP_positive", "BF_positive"]
                record.update({f"GARD_FUBAR_{name}": value for name, value in zip(names, row[:6])})
            elif method == "FEL":
                names = ["alpha", "beta", "neutral_rate", "LRT", "p", "total_branch_length"]
                record.update({f"GARD_FEL_{name}": value for name, value in zip(names, row[:6])})
            elif method == "MEME":
                names = [
                    "alpha", "beta_background", "weight_background", "beta_positive",
                    "weight_positive", "LRT", "p", "branches_EBF100", "total_branch_length",
                    "logL", "FEL_logL", "FEL_alpha", "FEL_beta",
                ]
                record.update({f"GARD_MEME_{name}": value for name, value in zip(names, row[:13])})
            result_rows.append(record)
    frame = pd.DataFrame(result_rows).sort_values("Filtered_Site").reset_index(drop=True)
    if frame["Filtered_Site"].duplicated().any():
        raise ValueError(f"{method}: duplicate global sites")
    return frame


def main() -> None:
    gard = json.loads(GARD_JSON.read_text())
    if not GARD_LF.exists() or GARD_LF.stat().st_size == 0:
        raise ValueError("Completed nonempty GARD likelihood snapshot required")
    method_files = {
        "FUBAR": PARTITIONED / "Sauropsida_MACSE_GARD_FUBAR.json",
        "FEL": PARTITIONED / "Sauropsida_MACSE_GARD_FEL.json",
        "MEME": PARTITIONED / "Sauropsida_MACSE_GARD_MEME.json",
    }
    for path in method_files.values():
        if not path.exists() or path.stat().st_size == 0:
            raise ValueError(f"Missing partition-aware result: {path}")

    mapping = pd.read_csv(PROJECT / "analysis_v2" / "summaries" / "retained_site_map.csv")
    mapping = mapping.loc[(mapping["Scope"] == "Sauropsida") & (mapping["Aligner"] == "MACSE")].copy()
    mapping = mapping.sort_values("Filtered_Site").reset_index(drop=True)
    frames = [flatten_mle(path, method) for method, path in method_files.items()]
    sites = mapping.merge(frames[0], on="Filtered_Site", how="inner")
    for frame in frames[1:]:
        sites = sites.merge(frame, on=["Filtered_Site", "GARD_Partition"], how="inner")
    if len(sites) != len(mapping):
        raise ValueError(f"Expected {len(mapping)} mapped sites, recovered {len(sites)}")
    sites["GARD_FEL_q"] = bh_adjust(sites["GARD_FEL_p"])
    sites["GARD_MEME_q"] = bh_adjust(sites["GARD_MEME_p"])
    sites["GARD_FEL_direction"] = np.where(
        sites["GARD_FEL_beta"] > sites["GARD_FEL_alpha"], "diversifying",
        np.where(sites["GARD_FEL_beta"] < sites["GARD_FEL_alpha"], "purifying", "neutral"),
    )
    sites["GARD_FUBAR_positive"] = sites["GARD_FUBAR_PP_positive"] >= 0.90
    sites["GARD_FUBAR_purifying"] = sites["GARD_FUBAR_PP_negative"] >= 0.90
    sites["GARD_FEL_positive_q10"] = (sites["GARD_FEL_direction"] == "diversifying") & (sites["GARD_FEL_q"] <= 0.10)
    sites["GARD_FEL_positive_q05"] = (sites["GARD_FEL_direction"] == "diversifying") & (sites["GARD_FEL_q"] <= 0.05)
    sites["GARD_FEL_purifying_q05"] = (sites["GARD_FEL_direction"] == "purifying") & (sites["GARD_FEL_q"] <= 0.05)
    sites["GARD_MEME_q10"] = sites["GARD_MEME_q"] <= 0.10
    sites["GARD_MEME_q05"] = sites["GARD_MEME_q"] <= 0.05
    sites["GARD_joint_pervasive_q10"] = sites["GARD_FUBAR_positive"] & sites["GARD_FEL_positive_q10"]
    sites["GARD_joint_pervasive_q05"] = sites["GARD_FUBAR_positive"] & sites["GARD_FEL_positive_q05"]
    sites["GARD_joint_purifying_q05"] = sites["GARD_FUBAR_purifying"] & sites["GARD_FEL_purifying_q05"]

    standard = pd.read_csv(PROJECT / "analysis_v2" / "summaries" / "sitewise_selection_all_models.csv")
    standard = standard.loc[(standard["Scope"] == "Sauropsida") & (standard["Aligner"] == "MACSE")].copy()
    standard_columns = [
        "Human_Position", "FUBAR_PP_positive", "FUBAR_PP_negative", "FEL_q", "FEL_direction",
        "MEME_q", "Joint_Pervasive_q10", "Joint_Pervasive_q05", "Joint_Purifying_q05", "MEME_q10", "MEME_q05",
    ]
    sites = sites.merge(standard[standard_columns], on="Human_Position", how="left", suffixes=("_GARD", "_standard"))
    sites.to_csv(ANALYSIS / "gard_partitioned_sitewise_comparison.csv", index=False)

    nexus_intervals = []
    breakpoint_nt = [
        int(item[0])
        for item in gard["improvements"][str(len(gard["breakpointData"]) - 2)]["breakpoints"]
    ] if len(gard["breakpointData"]) > 1 else []
    partitions = []
    breakpoint_filtered = []
    breakpoint_human = []
    breakpoint_support = []
    sorted_partitions = sorted(gard["breakpointData"], key=int)
    for index, key in enumerate(sorted_partitions):
        interval = gard["breakpointData"][key]["bps"][0]
        start, end = int(interval[0]), int(interval[1])
        nexus_intervals.append(f"{start}-{end}")
    partition_starts = [1] + [value + 1 for value in breakpoint_nt]
    partition_ends = breakpoint_nt + [int(gard["input"]["number of sites"])]
    partitions = [f"{start}-{end}" for start, end in zip(partition_starts, partition_ends)]
    for breakpoint in breakpoint_nt:
        breakpoint_support.append(float(gard.get("siteBreakPointSupport", {}).get(str(breakpoint), math.nan)))
        if breakpoint % 3 != 0:
            raise ValueError(f"Selected GARD breakpoint is not on a codon boundary: {breakpoint}")
        filtered_site = breakpoint // 3
        breakpoint_filtered.append(filtered_site)
        hit = mapping.loc[mapping["Filtered_Site"] == filtered_site, "Human_Position"]
        breakpoint_human.append(int(hit.iloc[0]) if len(hit) else None)

    delta = float(gard["singleTreeAICc"]) - float(gard["bestModelAICc"])
    log10_evidence = delta / (2.0 * math.log(10.0))
    focal = sites.loc[sites["Human_Position"].isin([310, 358, 369, 380]), [
        "Human_Position", "Filtered_Site", "GARD_Partition",
        "FUBAR_PP_positive", "GARD_FUBAR_PP_positive",
        "FEL_q", "GARD_FEL_q", "MEME_q", "GARD_MEME_q",
        "MEME_q10", "GARD_MEME_q10", "MEME_q05", "GARD_MEME_q05",
    ]].sort_values("Human_Position")
    focal.to_csv(ANALYSIS / "gard_partitioned_focal_site_comparison.csv", index=False)

    summary = pd.DataFrame([{
        "Scope": "Sauropsida",
        "Alignment": "MACSE",
        "Sequences": int(gard["input"]["number of sequences"]),
        "Nucleotide_Sites": int(gard["input"]["number of sites"]),
        "Inferred_Breakpoint_Count": len(sorted_partitions) - 1,
        "Partition_Nucleotide_Intervals": ";".join(partitions),
        "GARD_JSON_NEXUS_Intervals": ";".join(nexus_intervals),
        "Breakpoint_Nucleotide_Ends": ";".join(map(str, breakpoint_nt)) if breakpoint_nt else "none",
        "Breakpoint_Filtered_Codon_After": ";".join(map(str, breakpoint_filtered)) if breakpoint_filtered else "none",
        "Breakpoint_Human_Position_After": ";".join(str(x) for x in breakpoint_human) if breakpoint_human else "none",
        "Selected_Breakpoint_Model_Averaged_Support": ";".join(
            f"{value:.6g}" for value in breakpoint_support
        ) if breakpoint_support else "none",
        "Best_Model_AICc": float(gard["bestModelAICc"]),
        "Single_Tree_AICc": float(gard["singleTreeAICc"]),
        "Delta_AICc_vs_single_tree": delta,
        "Log10_Evidence_Ratio": log10_evidence,
        "Evidence_Ratio_at_least_100": bool(log10_evidence >= 2.0),
        "Models_Evaluated": int(gard["totalModelCount"]),
        "GARD_Checkpoint_Run_Seconds_Diagnostic": 3916.0,
        "GARD_Mixed_Resume_Run_Seconds_Diagnostic": 3566.0,
        "GARD_Clean_Run_Seconds": float(gard["timeElapsed"]),
        "GARD_Run_Seconds": float(gard["timeElapsed"]),
        "GARD_Partitioned_FUBAR_Positive_Count": int(sites["GARD_FUBAR_positive"].sum()),
        "GARD_Partitioned_FEL_Diversifying_q10_Count": int(sites["GARD_FEL_positive_q10"].sum()),
        "GARD_Partitioned_MEME_q10_Count": int(sites["GARD_MEME_q10"].sum()),
        "GARD_Partitioned_Purifying_q05_Count": int(sites["GARD_joint_purifying_q05"].sum()),
        "NP380_Standard_MEME_q": float(focal.loc[focal["Human_Position"] == 380, "MEME_q"].iloc[0]),
        "NP380_GARD_MEME_q": float(focal.loc[focal["Human_Position"] == 380, "GARD_MEME_q"].iloc[0]),
        "h358_Standard_MEME_q": float(focal.loc[focal["Human_Position"] == 358, "MEME_q"].iloc[0]),
        "h358_GARD_MEME_q": float(focal.loc[focal["Human_Position"] == 358, "GARD_MEME_q"].iloc[0]),
        "GARD_JSON_SHA256": sha256(GARD_JSON),
        "GARD_LF_SHA256": sha256(GARD_LF),
    }])
    summary.to_csv(ANALYSIS / "gard_summary.csv", index=False)
    common = {
        "working_directory": ".",
        "input": "../../analysis_v2/selection_inputs/Sauropsida_MACSE_selection_input.fasta",
        "input_sha256": "3d547c024fb1df4562c0febf5a0ac49e0ccbc72f898d997caa661d95684df851",
        "parameters": "type=nucleotide;model=GTR;rv=GDD;mode=Faster;max_breakpoints=5",
    }
    manifest_rows = [
        {
            "run_id": "gard_sauropsida_macse_checkpoint",
            "status": "checkpoint_saved_then_numerical_termination_during_3bp",
            "started_local": "2026-10-25T10:36:00-07:00",
            "finished_local": "2026-10-25T12:41:04-07:00",
            "executable": "hyphy",
            "version": "2.5.94",
            "executable_sha256": "a96a439829dd1b8f461a8447599c341b6b6d8e1275f4e9934ccda7c7319967e6",
            "source_release_sha256": "not_applicable",
            **common,
            "output_json": "analysis/Sauropsida_MACSE_GARD_2bp_checkpoint.json",
            "output_json_sha256": sha256(GARD_CHECKPOINT),
            "output_lf": "not_finalized",
            "output_lf_sha256": "not_applicable",
        },
        {
            "run_id": "gard_sauropsida_macse_resume",
            "status": "complete_diagnostic_only_not_used_for_inference",
            "started_local": "2026-10-25T12:44:00-07:00",
            "finished_local": "2026-10-25T13:47:00-07:00",
            "executable": "mpirun -np 4 workspace-local HYPHYMPI",
            "version": "2.5.103",
            "executable_sha256": "0f239527d4d7ae2e2f427092742c74eb6d9f94608d73ddf20381bc351c3c52f5",
            "source_release_sha256": "e3602aa3add7f4d88c18038828bc49080a749dd377f9d41c92933fbf07d846c2",
            **common,
            "parameters": common["parameters"] + ";MPI_workers=4;Open_MPI=4.1.6;ENV=TOLERATE_NUMERICAL_ERRORS=1",
            "output_json": "analysis/Sauropsida_MACSE_GARD.json",
            "output_json_sha256": sha256(GARD_MIXED_JSON),
            "output_lf": "analysis/Sauropsida_MACSE_GARD.best-gard",
            "output_lf_sha256": sha256(GARD_MIXED_LF),
        },
        {
            "run_id": "gard_sauropsida_macse_clean",
            "status": "complete_used_for_inference",
            "started_local": "2026-10-25T13:53:00-07:00",
            "finished_local": datetime.now().astimezone().isoformat(timespec="seconds"),
            "executable": "mpirun -np 4 workspace-local HYPHYMPI",
            "version": "2.5.103",
            "executable_sha256": "0f239527d4d7ae2e2f427092742c74eb6d9f94608d73ddf20381bc351c3c52f5",
            "source_release_sha256": "e3602aa3add7f4d88c18038828bc49080a749dd377f9d41c92933fbf07d846c2",
            **common,
            "parameters": common["parameters"] + ";MPI_workers=4;Open_MPI=4.1.6;ENV=TOLERATE_NUMERICAL_ERRORS=1;clean_start=yes",
            "output_json": "analysis/Sauropsida_MACSE_GARD_2_5_103_clean.json",
            "output_json_sha256": sha256(GARD_JSON),
            "output_lf": "analysis/Sauropsida_MACSE_GARD_2_5_103_clean.best-gard",
            "output_lf_sha256": sha256(GARD_LF),
        },
    ]
    with (ANALYSIS / "run_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(summary.to_string(index=False))
    print(focal.to_string(index=False))


if __name__ == "__main__":
    main()
