from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


PROJECT = Path(__file__).resolve().parents[1]
MATRIX = (
    PROJECT
    / "source_package"
    / "analysis"
    / "orthology_rebuild_20260826"
    / "orthology_audit"
    / "DNMT3L_orthology_evidence_matrix.csv"
)
STRUCTURE = PROJECT / "analysis_v2" / "structure" / "9MPP.cif"
OUTPUT = PROJECT / "analysis_v2" / "summaries" / "NP_037501_2_to_9MPP_Q9UJW3_crosswalk.csv"


def structure_sequence() -> str:
    payload = MMCIF2Dict(str(STRUCTURE))
    strands = payload["_entity_poly.pdbx_strand_id"]
    entities = payload["_entity_poly.entity_id"]
    sequences = payload["_entity_poly.pdbx_seq_one_letter_code_can"]
    accessions = dict(
        zip(payload["_struct_ref.entity_id"], payload["_struct_ref.pdbx_db_accession"])
    )
    matches = [
        index
        for index, (entity, strand) in enumerate(zip(entities, strands))
        if "N" in strand.split(",") and accessions.get(entity) == "Q9UJW3"
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one 9MPP chain-N/Q9UJW3 entity, found {matches}")
    return "".join(sequences[matches[0]].split()).upper()


def main() -> None:
    matrix = pd.read_csv(MATRIX, dtype={"gene_id": str})
    human = matrix.loc[(matrix["species"] == "Homo sapiens") & matrix["Selection_Ready"].astype(bool)]
    if len(human) != 1:
        raise ValueError(f"Expected one selection-ready human record, found {len(human)}")
    np_sequence = str(human.iloc[0]["protein_seq"]).upper()
    q9ujw3_sequence = structure_sequence()
    if len(np_sequence) != 387 or len(q9ujw3_sequence) != 386:
        raise ValueError(
            f"Unexpected reference lengths: NP_037501.2={len(np_sequence)}, Q9UJW3={len(q9ujw3_sequence)}"
        )

    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    alignment = aligner.align(np_sequence, q9ujw3_sequence)[0]

    mapping: dict[int, int] = {}
    np_blocks, q9ujw3_blocks = alignment.aligned
    for (np_start, np_end), (q_start, q_end) in zip(np_blocks, q9ujw3_blocks):
        if np_end - np_start != q_end - q_start:
            raise ValueError("Unequal aligned block lengths in reference crosswalk")
        for offset in range(int(np_end - np_start)):
            mapping[int(np_start + offset + 1)] = int(q_start + offset + 1)

    # NP332-333 is an SS run aligned to one Q9UJW3 serine. Either gap
    # placement has the same score; use the established rightmost convention.
    if (
        mapping.get(332) is None
        and mapping.get(333) == 332
        and np_sequence[331:333] == "SS"
        and q9ujw3_sequence[331] == "S"
    ):
        mapping[332] = 332
        del mapping[333]

    rows: list[dict[str, object]] = []
    for np_position, np_residue in enumerate(np_sequence, start=1):
        q_position = mapping.get(np_position)
        q_residue = q9ujw3_sequence[q_position - 1] if q_position is not None else ""
        rows.append(
            {
                "NP_037501_2_Position": np_position,
                "NP_037501_2_Residue": np_residue,
                "Q9UJW3_9MPP_Entity_Position": q_position if q_position is not None else "",
                "Q9UJW3_9MPP_Entity_Residue": q_residue,
                "Relationship": (
                    "NP_specific_insertion"
                    if q_position is None
                    else "exact"
                    if np_residue == q_residue
                    else "substitution"
                ),
            }
        )

    by_np = {int(row["NP_037501_2_Position"]): row for row in rows}
    expected = {
        333: ("", "NP_specific_insertion"),
        349: (348, "exact"),
        352: (351, "exact"),
        358: (357, "exact"),
        380: (379, "exact"),
    }
    for np_position, (q_position, relationship) in expected.items():
        row = by_np[np_position]
        if row["Q9UJW3_9MPP_Entity_Position"] != q_position or row["Relationship"] != relationship:
            raise ValueError(f"Unexpected coordinate mapping at NP{np_position}: {row}")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUTPUT.relative_to(PROJECT)} ({len(rows)} NP positions)")
    print("Verified NP333 insertion and NP349->348, NP352->351, NP358->357, NP380->379 mappings")


if __name__ == "__main__":
    main()
