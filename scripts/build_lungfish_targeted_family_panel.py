#!/usr/bin/env python3
"""Add the lungfish translated segment to the established targeted DNMT3 panel.

The internal stop is represented as an unknown residue in the full-segment input.
A second input ends immediately before the stop as a sensitivity analysis.
"""
from __future__ import annotations

from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
LOCUS_DIR = PROJECT / "analysis/locus_history_337"
BASE_PANEL = LOCUS_DIR / "targeted_partial_candidate_panel.fasta"
LUNGFISH = LOCUS_DIR / "lungfish_candidate_hit_peptide.fasta"
FULL_OUTPUT = LOCUS_DIR / "lungfish_targeted_candidate_panel_full_masked_stop.fasta"
PRESTOP_OUTPUT = LOCUS_DIR / "lungfish_targeted_candidate_panel_prestop.fasta"
FULL_CANDIDATE_OUTPUT = LOCUS_DIR / "lungfish_candidate_full_masked_stop.fasta"
PRESTOP_CANDIDATE_OUTPUT = LOCUS_DIR / "lungfish_candidate_prestop.fasta"
LUNGFISH_ID = "LUNGFISH_candidate_translated_segment"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    sequence: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(sequence).upper()))
            name = line[1:].split()[0]
            sequence = []
        else:
            if name is None:
                raise ValueError(f"Sequence before header in {path}")
            sequence.append(line)
    if name is not None:
        records.append((name, "".join(sequence).upper()))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    panel = read_fasta(BASE_PANEL)
    if len(panel) != 13:
        raise ValueError(f"Expected the established 13-sequence panel, found {len(panel)}")
    candidate_records = read_fasta(LUNGFISH)
    if len(candidate_records) != 1 or candidate_records[0][0] != LUNGFISH_ID:
        raise ValueError("Unexpected lungfish candidate FASTA record")
    candidate = candidate_records[0][1]
    if len(candidate) != 244 or candidate.count("*") != 1:
        raise ValueError(
            f"Expected a 244-aa translation with one internal stop, got length={len(candidate)} "
            f"and stops={candidate.count('*')}"
        )
    stop_index = candidate.index("*")
    if stop_index != 229:
        raise ValueError(f"Expected the stop at one-based amino acid 230, got {stop_index + 1}")

    existing_ids = {name for name, _ in panel}
    if LUNGFISH_ID in existing_ids:
        raise ValueError("Lungfish candidate is already present in the base panel")

    full_candidate = (LUNGFISH_ID, candidate.replace("*", "X"))
    prestop_candidate = (LUNGFISH_ID, candidate[:stop_index])
    write_fasta(FULL_OUTPUT, panel + [full_candidate])
    write_fasta(PRESTOP_OUTPUT, panel + [prestop_candidate])
    write_fasta(FULL_CANDIDATE_OUTPUT, [full_candidate])
    write_fasta(PRESTOP_CANDIDATE_OUTPUT, [prestop_candidate])
    print(f"Wrote {FULL_OUTPUT.relative_to(PROJECT)} with masked stop (X), 244 residues")
    print(f"Wrote {PRESTOP_OUTPUT.relative_to(PROJECT)} with residues 1-{stop_index}")


if __name__ == "__main__":
    main()
