#!/usr/bin/env python3
from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.stats import mannwhitneyu


PROJECT = Path(__file__).resolve().parents[1]
ALIGNMENT_ROOT = PROJECT / "analysis/codon_alignments_337/Placentalia"
RESULTS = PROJECT / "results"
HUMAN = "Homo_sapiens"
ALIGNMENTS = {
    "MACSE": ALIGNMENT_ROOT / "MACSE/Placentalia_MACSE_NT.fasta",
    "MAFFT": ALIGNMENT_ROOT / "MAFFT/Placentalia_MAFFT_NT.fasta",
    "PRANK": ALIGNMENT_ROOT / "PRANK/Placentalia_PRANK.best.fas",
    "MUSCLE5": ALIGNMENT_ROOT / "MUSCLE5/Placentalia_MUSCLE5_NT.fasta",
    "COBALT": ALIGNMENT_ROOT / "COBALT/Placentalia_COBALT_NT.fasta",
}


def read_fasta(path: Path) -> dict[str, str]:
    records = {
        record.id: str(record.seq).upper().replace("U", "T")
        for record in SeqIO.parse(path, "fasta")
    }
    if len(records) != 220:
        raise ValueError(f"Expected 220 records in {path}, found {len(records)}")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1 or next(iter(lengths)) % 3:
        raise ValueError(f"Nonrectangular or noncodon alignment: {path}")
    return records


def resolved_codon(codon: str) -> bool:
    return (
        len(codon) == 3
        and not (set(codon) - set("ACGT"))
        and str(Seq(codon).translate()) != "*"
    )


def human_coordinate_columns(sequence: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    position = 0
    for column in range(len(sequence) // 3):
        codon = sequence[column * 3 : (column + 1) * 3]
        if codon == "---":
            continue
        if not resolved_codon(codon):
            if set(codon) <= set("ACGT") and str(Seq(codon).translate()) == "*":
                continue
            raise ValueError(f"Unresolved human codon at alignment column {column + 1}")
        position += 1
        mapping[position] = column
    if position != 387:
        raise ValueError(f"Mapped {position}, not 387, human residues")
    return mapping


def domain(position: int) -> str:
    if 41 <= position <= 173:
        return "ADD"
    if 178 <= position <= 380:
        return "C_terminal_MTase_like"
    if position <= 40:
        return "N_terminal"
    if position <= 177:
        return "linker"
    return "C_terminal_tail"


def structural_region(position: int) -> str:
    if 341 <= position <= 355:
        return "switching_helix"
    if 356 <= position <= 360:
        return "adjacent_loop"
    return "other"


def entropy(counts: Counter[str]) -> float:
    total = sum(counts.values())
    if not total:
        return float("nan")
    return -sum((count / total) * math.log2(count / total) for count in counts.values())


def cliffs_delta(first: np.ndarray, second: np.ndarray) -> float:
    statistic = mannwhitneyu(first, second, alternative="two-sided").statistic
    return float((2 * statistic) / (len(first) * len(second)) - 1)


def bh(values: list[float]) -> list[float]:
    array = np.asarray(values, dtype=float)
    order = np.argsort(array)
    ranked = array[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output.tolist()


def main() -> None:
    mask = pd.read_csv(
        RESULTS / "placental_five_aligner_consensus_mask_by_site.tsv", sep="\t"
    )
    retained = set(
        mask.loc[
            mask["retained_by_80pct_all_method_gate"].eq(True), "human_position"
        ].astype(int)
    )
    if len(retained) != 367:
        raise ValueError(f"Expected 367 shared positions, found {len(retained)}")

    site_rows: list[dict[str, object]] = []
    for method, path in ALIGNMENTS.items():
        records = read_fasta(path)
        mapping = human_coordinate_columns(records[HUMAN])
        for position in sorted(retained):
            column = mapping[position]
            counts: Counter[str] = Counter()
            for sequence in records.values():
                codon = sequence[column * 3 : (column + 1) * 3]
                if resolved_codon(codon):
                    counts[str(Seq(codon).translate())] += 1
            resolved = sum(counts.values())
            site_rows.append(
                {
                    "method": method,
                    "human_position": position,
                    "domain": domain(position),
                    "structural_region": structural_region(position),
                    "resolved_records": resolved,
                    "occupancy": resolved / len(records),
                    "amino_acid_entropy": entropy(counts),
                    "major_residue": counts.most_common(1)[0][0] if counts else "",
                    "major_residue_frequency": (
                        counts.most_common(1)[0][1] / resolved if resolved else float("nan")
                    ),
                    "residue_counts": ";".join(
                        f"{residue}:{count}" for residue, count in sorted(counts.items())
                    ),
                }
            )

    sites = pd.DataFrame(site_rows)
    sites.to_csv(
        RESULTS / "placental_five_aligner_constraint_by_method_site.tsv",
        sep="\t",
        index=False,
    )

    consensus = (
        sites.groupby(["human_position", "domain", "structural_region"], as_index=False)
        .agg(
            median_entropy=("amino_acid_entropy", "median"),
            minimum_occupancy=("occupancy", "min"),
            minimum_major_residue_frequency=("major_residue_frequency", "min"),
            methods=("method", "nunique"),
        )
    )
    consensus.to_csv(
        RESULTS / "placental_five_aligner_constraint_by_site.tsv",
        sep="\t",
        index=False,
    )

    comparison_rows: list[dict[str, object]] = []
    for label, frame in [("FIVE_ALIGNER_MEDIAN", consensus)]:
        add = frame.loc[frame["domain"].eq("ADD"), "median_entropy"].to_numpy()
        cterminal = frame.loc[
            frame["domain"].eq("C_terminal_MTase_like"), "median_entropy"
        ].to_numpy()
        test = mannwhitneyu(add, cterminal, alternative="two-sided")
        comparison_rows.append(
            {
                "analysis": label,
                "comparison": "ADD_vs_C_terminal_MTase_like",
                "first_n": len(add),
                "second_n": len(cterminal),
                "first_median_entropy": float(np.median(add)),
                "second_median_entropy": float(np.median(cterminal)),
                "cliffs_delta": cliffs_delta(add, cterminal),
                "mann_whitney_u": float(test.statistic),
                "p": float(test.pvalue),
                "q_within_analysis": float(test.pvalue),
            }
        )

        background = frame.loc[
            frame["domain"].eq("C_terminal_MTase_like")
            & frame["structural_region"].eq("other"),
            "median_entropy",
        ].to_numpy()
        local_rows: list[dict[str, object]] = []
        for region in ("switching_helix", "adjacent_loop"):
            focal = frame.loc[
                frame["domain"].eq("C_terminal_MTase_like")
                & frame["structural_region"].eq(region),
                "median_entropy",
            ].to_numpy()
            test = mannwhitneyu(focal, background, alternative="two-sided")
            local_rows.append(
                {
                    "analysis": label,
                    "comparison": f"{region}_vs_other_C_terminal",
                    "first_n": len(focal),
                    "second_n": len(background),
                    "first_median_entropy": float(np.median(focal)),
                    "second_median_entropy": float(np.median(background)),
                    "cliffs_delta": cliffs_delta(focal, background),
                    "mann_whitney_u": float(test.statistic),
                    "p": float(test.pvalue),
                }
            )
        local_q = bh([row["p"] for row in local_rows])
        for row, q_value in zip(local_rows, local_q):
            row["q_within_analysis"] = q_value
            comparison_rows.append(row)

    for method, frame in sites.groupby("method"):
        add = frame.loc[frame["domain"].eq("ADD"), "amino_acid_entropy"].to_numpy()
        cterminal = frame.loc[
            frame["domain"].eq("C_terminal_MTase_like"), "amino_acid_entropy"
        ].to_numpy()
        test = mannwhitneyu(add, cterminal, alternative="two-sided")
        comparison_rows.append(
            {
                "analysis": method,
                "comparison": "ADD_vs_C_terminal_MTase_like",
                "first_n": len(add),
                "second_n": len(cterminal),
                "first_median_entropy": float(np.median(add)),
                "second_median_entropy": float(np.median(cterminal)),
                "cliffs_delta": cliffs_delta(add, cterminal),
                "mann_whitney_u": float(test.statistic),
                "p": float(test.pvalue),
                "q_within_analysis": float(test.pvalue),
            }
        )

    comparisons = pd.DataFrame(comparison_rows)
    comparisons.to_csv(
        RESULTS / "placental_five_aligner_constraint_tests.tsv",
        sep="\t",
        index=False,
    )
    print(comparisons.to_string(index=False))


if __name__ == "__main__":
    main()
