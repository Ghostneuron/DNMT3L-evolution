from __future__ import annotations

from collections import Counter
from pathlib import Path

import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices


PROJECT = Path(__file__).resolve().parents[1]
MATRIX = (
    PROJECT
    / "source_package"
    / "analysis"
    / "orthology_rebuild_20260826"
    / "orthology_audit"
    / "DNMT3L_orthology_evidence_matrix.csv"
)
OUT = PROJECT / "analysis_v2" / "summaries"
START = 341
END = 360
STRUCTURE_INFORMED = {349: "UniProt/PDB Q348", 352: "UniProt/PDB Q351"}
MAMMALIAN_CLADES = {
    "Afrotheria",
    "Euarchontoglires",
    "Laurasiatheria",
    "Marsupialia",
    "Xenarthra",
}


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    aligner.target_end_gap_score = -2.0
    aligner.query_end_gap_score = -2.0
    return aligner


def map_positions(reference: str, query: str, aligner: PairwiseAligner) -> tuple[dict[int, tuple[int, str]], float]:
    alignment = aligner.align(reference, query)[0]
    reference_blocks, query_blocks = alignment.aligned
    mapped: dict[int, tuple[int, str]] = {}
    for (ref_start, ref_end), (query_start, query_end) in zip(reference_blocks, query_blocks):
        for offset, query_aa in enumerate(query[query_start:query_end]):
            ref_position = int(ref_start + offset + 1)
            if START <= ref_position <= END:
                mapped[ref_position] = (int(query_start + offset + 1), query_aa)
    return mapped, float(alignment.score)


def count_string(values: pd.Series) -> str:
    counts = Counter(values.dropna().astype(str))
    return ";".join(f"{residue}:{counts[residue]}" for residue in sorted(counts))


def summarize_groups(mapped: pd.DataFrame, group_column: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (group_name, position), group in mapped.groupby([group_column, "Human_Position"], sort=True):
        resolved = group["Query_Residue"].dropna()
        counts = Counter(resolved.astype(str))
        major_residue, major_count = counts.most_common(1)[0] if counts else (None, 0)
        rows.append(
            {
                group_column: group_name,
                "Human_Position": int(position),
                "Human_Residue": group["Human_Residue"].iloc[0],
                "Structural_Annotation": group["Structural_Annotation"].iloc[0],
                "Records": len(group),
                "Resolved_Records": len(resolved),
                "Occupancy": len(resolved) / len(group),
                "Major_Residue": major_residue,
                "Major_Residue_Frequency": major_count / len(resolved) if len(resolved) else None,
                "Human_State_Count": counts.get(group["Human_Residue"].iloc[0], 0),
                "Human_State_Frequency": counts.get(group["Human_Residue"].iloc[0], 0) / len(resolved) if len(resolved) else None,
                "Amino_Acid_State_Count": len(counts),
                "Amino_Acid_Counts": count_string(resolved),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    records = pd.read_csv(MATRIX, dtype={"gene_id": str})
    ready = records.loc[records["Selection_Ready"].astype(bool)].copy()
    if len(ready) != 269:
        raise ValueError(f"Expected 269 selection-ready records, found {len(ready)}")
    human = ready.loc[ready["species"] == "Homo sapiens"]
    if len(human) != 1:
        raise ValueError(f"Expected one selection-ready human record, found {len(human)}")
    reference = str(human.iloc[0]["protein_seq"]).upper()
    if len(reference) != 387 or reference[348] != "Q" or reference[351] != "Q":
        raise ValueError("Unexpected NP_037501.2 coordinate anchor")

    aligner = make_aligner()
    rows: list[dict[str, object]] = []
    for record in ready.to_dict(orient="records"):
        mapped, score = map_positions(reference, str(record["protein_seq"]).upper(), aligner)
        for position in range(START, END + 1):
            query_position, residue = mapped.get(position, (None, None))
            rows.append(
                {
                    "safe_id": record["safe_id"],
                    "species": record["species"],
                    "gene_id": record["gene_id"],
                    "Clade": record["Clade"],
                    "Human_Position": position,
                    "Human_Residue": reference[position - 1],
                    "Query_Position": query_position,
                    "Query_Residue": residue,
                    "Pairwise_Alignment_Score": score,
                    "Structural_Annotation": STRUCTURE_INFORMED.get(position, "switching_helix" if position <= 355 else "adjacent_loop"),
                }
            )
    mapped = pd.DataFrame(rows)
    mapped.to_csv(OUT / "structure_informed_residue_states_by_record.csv", index=False)

    summary = summarize_groups(mapped, "Clade")
    summary.to_csv(OUT / "structure_informed_residue_summary_by_clade.csv", index=False)
    focal = summary.loc[summary["Human_Position"].isin(STRUCTURE_INFORMED)].copy()
    focal.to_csv(OUT / "acidic_patch_counterpart_residue_summary.csv", index=False)

    mapped["Lineage_Group"] = mapped["Clade"].where(
        ~mapped["Clade"].isin(MAMMALIAN_CLADES),
        "Mammalia",
    )
    lineage_summary = summarize_groups(mapped, "Lineage_Group")
    lineage_summary.to_csv(OUT / "structure_informed_residue_summary_by_lineage.csv", index=False)
    lineage_focal = lineage_summary.loc[lineage_summary["Human_Position"].isin(STRUCTURE_INFORMED)].copy()
    lineage_focal.to_csv(OUT / "acidic_patch_counterpart_residue_summary_by_lineage.csv", index=False)

    print(f"Mapped {len(ready)} records across NP positions {START}-{END}")
    print(focal.to_string(index=False))
    print("\nLineage aggregates")
    print(lineage_focal.to_string(index=False))


if __name__ == "__main__":
    main()
