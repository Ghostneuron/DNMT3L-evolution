from __future__ import annotations

import hashlib
import math
import random
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import BiopythonWarning
from Bio.Data import CodonTable
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds


warnings.simplefilter("ignore", BiopythonWarning)

PROJECT = Path(__file__).resolve().parents[1]
ALIGNMENT_DIR = PROJECT / "analysis_v2" / "selection_inputs"
OUT = PROJECT / "analysis_v2" / "saturation_diagnostics"
SEED = 20261006
MAX_PAIRS = 3000
MIN_VALID_CODONS = 100
TRANSITIONS = {frozenset(("A", "G")), frozenset(("C", "T"))}
STOPS = set(CodonTable.unambiguous_dna_by_id[1].stop_codons)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks).upper()))
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    return records


def valid_pair_codons(first: str, second: str) -> tuple[list[str], list[str]]:
    first_codons: list[str] = []
    second_codons: list[str] = []
    for start in range(0, len(first), 3):
        codon_a = first[start : start + 3]
        codon_b = second[start : start + 3]
        if (
            len(codon_a) == 3
            and len(codon_b) == 3
            and set(codon_a + codon_b) <= set("ACGT")
            and codon_a not in STOPS
            and codon_b not in STOPS
        ):
            first_codons.append(codon_a)
            second_codons.append(codon_b)
    return first_codons, second_codons


def substitution_proportions(
    first_codons: list[str], second_codons: list[str], position: int
) -> tuple[int, float, float, float, float]:
    transitions = 0
    transversions = 0
    valid = 0
    for codon_a, codon_b in zip(first_codons, second_codons):
        base_a = codon_a[position]
        base_b = codon_b[position]
        valid += 1
        if base_a == base_b:
            continue
        if frozenset((base_a, base_b)) in TRANSITIONS:
            transitions += 1
        else:
            transversions += 1
    if not valid:
        return 0, math.nan, math.nan, math.nan, math.nan
    transition_p = transitions / valid
    transversion_p = transversions / valid
    p_distance = transition_p + transversion_p
    first_term = 1 - 2 * transition_p - transversion_p
    second_term = 1 - 2 * transversion_p
    k80 = (
        -0.5 * math.log(first_term) - 0.25 * math.log(second_term)
        if first_term > 0 and second_term > 0
        else math.nan
    )
    return valid, transition_p, transversion_p, p_distance, k80


def sampled_pairs(sequence_count: int) -> list[tuple[int, int]]:
    pairs = list(combinations(range(sequence_count), 2))
    if len(pairs) <= MAX_PAIRS:
        return pairs
    random.Random(SEED + sequence_count).shuffle(pairs)
    return sorted(pairs[:MAX_PAIRS])


def finite_percentile(values: pd.Series, percentile: float) -> float:
    finite = values[np.isfinite(values)]
    return float(np.percentile(finite, percentile)) if len(finite) else math.nan


def diagnostic_label(third_k80_failure: float, ds_ge_1: float, ds_ge_2: float) -> str:
    if third_k80_failure > 10 or (math.isfinite(ds_ge_2) and ds_ge_2 > 25):
        return "substantial multiple-hit/saturation risk"
    if third_k80_failure > 1 or (math.isfinite(ds_ge_1) and ds_ge_1 > 25):
        return "moderate multiple-hit/saturation risk"
    return "limited multiple-hit signal in sampled pairs"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pair_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    paths = sorted(
        path
        for path in ALIGNMENT_DIR.glob("*_selection_input.fasta")
        if not path.name.startswith("._")
    )
    if len(paths) != 12:
        raise ValueError(f"Expected 12 selection alignments, found {len(paths)}")

    for path in paths:
        dataset = path.name.removesuffix("_selection_input.fasta")
        scope, aligner = dataset.rsplit("_", 1)
        records = read_fasta(path)
        all_pair_count = len(records) * (len(records) - 1) // 2
        pairs = sampled_pairs(len(records))
        current_rows: list[dict[str, object]] = []
        for first_index, second_index in pairs:
            first_name, first_sequence = records[first_index]
            second_name, second_sequence = records[second_index]
            first_codons, second_codons = valid_pair_codons(first_sequence, second_sequence)
            row: dict[str, object] = {
                "Scope": scope,
                "Aligner": aligner,
                "Dataset": dataset,
                "First": first_name,
                "Second": second_name,
                "Valid_Codons": len(first_codons),
                "NG86_dN": math.nan,
                "NG86_dS": math.nan,
                "NG86_Status": "too_few_valid_codons",
            }
            for position in range(3):
                valid, transition, transversion, distance, k80 = substitution_proportions(
                    first_codons, second_codons, position
                )
                prefix = f"Pos{position + 1}"
                row[f"{prefix}_Valid_N"] = valid
                row[f"{prefix}_Transition_P"] = transition
                row[f"{prefix}_Transversion_P"] = transversion
                row[f"{prefix}_PDistance"] = distance
                row[f"{prefix}_K80"] = k80
            if len(first_codons) >= MIN_VALID_CODONS:
                try:
                    estimate = cal_dn_ds(
                        CodonSeq("".join(first_codons)),
                        CodonSeq("".join(second_codons)),
                        method="NG86",
                    )
                    if estimate is None:
                        raise ArithmeticError("NG86 returned no estimate")
                    d_n, d_s = estimate
                    row["NG86_dN"] = float(d_n)
                    row["NG86_dS"] = float(d_s)
                    row["NG86_Status"] = (
                        "finite"
                        if math.isfinite(float(d_n)) and math.isfinite(float(d_s))
                        else "nonfinite"
                    )
                except (ArithmeticError, TypeError, ValueError, ZeroDivisionError) as error:
                    row["NG86_Status"] = f"failed:{type(error).__name__}"
            current_rows.append(row)
            pair_rows.append(row)

        table = pd.DataFrame(current_rows)
        finite_ds = table["NG86_dS"].where(np.isfinite(table["NG86_dS"]))
        finite_count = int(finite_ds.notna().sum())
        third_k80_failure = float(100 * table["Pos3_K80"].isna().mean())
        ds_ge_1 = float(100 * (finite_ds >= 1).sum() / finite_count) if finite_count else math.nan
        ds_ge_2 = float(100 * (finite_ds >= 2).sum() / finite_count) if finite_count else math.nan
        correlation_table = table.loc[
            np.isfinite(table["NG86_dS"]) & np.isfinite(table["Pos3_PDistance"]),
            ["NG86_dS", "Pos3_PDistance"],
        ]
        spearman = (
            float(correlation_table.corr(method="spearman").iloc[0, 1])
            if len(correlation_table) >= 3
            else math.nan
        )
        summary_rows.append(
            {
                "Scope": scope,
                "Aligner": aligner,
                "Dataset": dataset,
                "Alignment_SHA256": sha256(path),
                "Sequence_Count": len(records),
                "All_Pair_Count": all_pair_count,
                "Sampled_Pair_Count": len(table),
                "Minimum_Valid_Codons": int(table["Valid_Codons"].min()),
                "Median_Valid_Codons": float(table["Valid_Codons"].median()),
                "NG86_Finite_N": finite_count,
                "NG86_Failure_or_Nonfinite_Pct": 100 * (1 - finite_count / len(table)),
                "NG86_dS_Median": finite_percentile(finite_ds, 50),
                "NG86_dS_P90": finite_percentile(finite_ds, 90),
                "NG86_dS_ge_1_Pct": ds_ge_1,
                "NG86_dS_ge_2_Pct": ds_ge_2,
                "Pos1_PDistance_Median": float(table["Pos1_PDistance"].median()),
                "Pos2_PDistance_Median": float(table["Pos2_PDistance"].median()),
                "Pos3_PDistance_Median": float(table["Pos3_PDistance"].median()),
                "Pos3_PDistance_P90": finite_percentile(table["Pos3_PDistance"], 90),
                "Pos3_K80_Undefined_Pct": third_k80_failure,
                "Spearman_NG86_dS_vs_Pos3_PDistance": spearman,
                "Diagnostic_Result": diagnostic_label(third_k80_failure, ds_ge_1, ds_ge_2),
            }
        )
        print(f"{dataset}: {summary_rows[-1]['Diagnostic_Result']}", flush=True)

    pd.DataFrame(pair_rows).to_csv(OUT / "pairwise_codon_saturation_diagnostics.csv", index=False)
    summary = pd.DataFrame(summary_rows).sort_values(["Scope", "Aligner"])
    summary.to_csv(OUT / "codon_saturation_summary.csv", index=False)
    print(summary[["Dataset", "NG86_dS_Median", "NG86_dS_P90", "NG86_dS_ge_1_Pct", "NG86_dS_ge_2_Pct", "Pos3_K80_Undefined_Pct", "Diagnostic_Result"]].to_string(index=False))


if __name__ == "__main__":
    main()
