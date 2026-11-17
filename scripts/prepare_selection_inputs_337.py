#!/usr/bin/env python3
from __future__ import annotations

import hashlib
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


PROJECT = Path(__file__).resolve().parents[1]
ALIGNMENT_ROOT = PROJECT / "analysis/codon_alignments_337/Placentalia"
OUT = PROJECT / "analysis/selection_reanalysis_337/inputs"
RESULTS = PROJECT / "results"
HUMAN = "Homo_sapiens"
MIN_OCCUPANCY = 0.8
ALIGNMENTS = {
    "MACSE": ALIGNMENT_ROOT / "MACSE/Placentalia_MACSE_NT.fasta",
    "MAFFT": ALIGNMENT_ROOT / "MAFFT/Placentalia_MAFFT_NT.fasta",
    "PRANK": ALIGNMENT_ROOT / "PRANK/Placentalia_PRANK.best.fas",
    "MUSCLE5": ALIGNMENT_ROOT / "MUSCLE5/Placentalia_MUSCLE5_NT.fasta",
    "COBALT": ALIGNMENT_ROOT / "COBALT/Placentalia_COBALT_NT.fasta",
    "CONSENSUS": (
        ALIGNMENT_ROOT
        / "five_aligner_consensus_filtered/Placentalia_MACSE_five_aligner_consensus_filtered.fasta"
    ),
}
CONSENSUS_MASK = RESULTS / "placental_five_aligner_consensus_mask_by_site.tsv"


def read_fasta(path: Path) -> dict[str, str]:
    records = {record.id: str(record.seq).upper().replace("U", "T") for record in SeqIO.parse(path, "fasta")}
    if len(records) != 220:
        raise ValueError(f"Expected 220 records in {path}, found {len(records)}")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1 or next(iter(lengths)) % 3:
        raise ValueError(f"Nonrectangular or noncodon alignment: {path}")
    return records


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolved_codon(codon: str) -> bool:
    return len(codon) == 3 and not (set(codon) - set("ACGT")) and str(Seq(codon).translate()) != "*"


def human_coordinate_columns(sequence: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    position = 0
    for column in range(len(sequence) // 3):
        codon = sequence[column * 3 : (column + 1) * 3]
        if codon == "---":
            continue
        if set(codon) <= set("ACGT") and str(Seq(codon).translate()) == "*":
            continue
        if not resolved_codon(codon):
            raise ValueError(f"Unresolved human codon at source alignment column {column + 1}")
        position += 1
        mapping[position] = column
    if position != 387:
        raise ValueError(f"Human coordinate mapping returned {position}, not 387, residues")
    return mapping


def write_fasta(records: dict[str, str], path: Path) -> None:
    with path.open("w") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def collapse_exact_duplicates(records: dict[str, str]) -> tuple[dict[str, str], list[dict[str, object]]]:
    by_sequence: dict[str, list[str]] = defaultdict(list)
    for name, sequence in records.items():
        by_sequence[sequence].append(name)
    collapsed: dict[str, str] = {}
    manifest: list[dict[str, object]] = []
    for sequence, names in sorted(by_sequence.items(), key=lambda item: sorted(item[1])[0]):
        names = sorted(names)
        representative = HUMAN if HUMAN in names else names[0]
        collapsed[representative] = sequence
        for name in names:
            manifest.append(
                {
                    "representative": representative,
                    "member": name,
                    "group_size": len(names),
                    "removed_as_duplicate": name != representative,
                }
            )
    return dict(sorted(collapsed.items())), manifest


def raw_method_filter(method: str, records: dict[str, str]) -> tuple[dict[str, str], list[dict[str, object]]]:
    mapping = human_coordinate_columns(records[HUMAN])
    retained: list[tuple[int, int, float]] = []
    for position, column in mapping.items():
        values = [sequence[column * 3 : (column + 1) * 3] for sequence in records.values()]
        occupancy = sum(resolved_codon(value) for value in values) / len(values)
        if occupancy >= MIN_OCCUPANCY:
            retained.append((position, column, occupancy))
    filtered = {
        name: "".join(
            (
                sequence[column * 3 : (column + 1) * 3]
                if resolved_codon(sequence[column * 3 : (column + 1) * 3])
                else "NNN"
            )
            for _, column, _ in retained
        )
        for name, sequence in records.items()
    }
    crosswalk = [
        {
            "analysis_id": method,
            "human_position": position,
            "source_alignment_codon_column": column + 1,
            "filtered_codon_site": site,
            "pre_filter_occupancy": occupancy,
            "post_filter_occupancy": occupancy,
        }
        for site, (position, column, occupancy) in enumerate(retained, start=1)
    ]
    return filtered, crosswalk


def consensus_filter(records: dict[str, str]) -> tuple[dict[str, str], list[dict[str, object]]]:
    mask = pd.read_csv(CONSENSUS_MASK, sep="\t")
    retained_original_positions = list(
        mask.loc[mask["retained_by_80pct_all_method_gate"].eq(True), "human_position"].astype(int)
    )
    if len(retained_original_positions) * 3 != len(records[HUMAN]):
        raise ValueError("Consensus mask and consensus alignment dimensions disagree")
    keep: list[tuple[int, int, float]] = []
    for source_index, human_position in enumerate(retained_original_positions):
        values = [sequence[source_index * 3 : (source_index + 1) * 3] for sequence in records.values()]
        occupancy = sum(resolved_codon(value) for value in values) / len(values)
        if occupancy >= MIN_OCCUPANCY:
            keep.append((human_position, source_index, occupancy))
    filtered = {
        name: "".join(sequence[index * 3 : (index + 1) * 3] for _, index, _ in keep)
        for name, sequence in records.items()
    }
    crosswalk = [
        {
            "analysis_id": "CONSENSUS",
            "human_position": position,
            "source_alignment_codon_column": index + 1,
            "filtered_codon_site": site,
            "pre_filter_occupancy": occupancy,
            "post_filter_occupancy": occupancy,
        }
        for site, (position, index, occupancy) in enumerate(keep, start=1)
    ]
    return filtered, crosswalk


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, object]] = []
    all_crosswalks: list[dict[str, object]] = []
    all_duplicates: list[dict[str, object]] = []
    for method, source_path in ALIGNMENTS.items():
        records = read_fasta(source_path)
        if method == "CONSENSUS":
            filtered, crosswalk = consensus_filter(records)
        else:
            filtered, crosswalk = raw_method_filter(method, records)
        collapsed, duplicate_manifest = collapse_exact_duplicates(filtered)
        filtered_path = OUT / f"Placentalia_{method}_human_coordinate_occupancy80.fasta"
        deduplicated_path = OUT / f"Placentalia_{method}_human_coordinate_occupancy80_deduplicated.fasta"
        write_fasta(dict(sorted(filtered.items())), filtered_path)
        write_fasta(collapsed, deduplicated_path)
        for row in crosswalk:
            all_crosswalks.append(row)
        for row in duplicate_manifest:
            all_duplicates.append({"analysis_id": method, **row})
        site_count = len(crosswalk)
        h358 = [row for row in crosswalk if row["human_position"] == 358]
        summary_rows.append(
            {
                "analysis_id": method,
                "primary_or_sensitivity": "PRIMARY" if method == "CONSENSUS" else "ALIGNER_SENSITIVITY",
                "source_alignment": str(source_path.relative_to(PROJECT)),
                "source_sha256": sha256(source_path),
                "taxa_before_deduplication": len(filtered),
                "taxa_after_deduplication": len(collapsed),
                "exact_duplicate_records_removed": len(filtered) - len(collapsed),
                "retained_human_positions": site_count,
                "h358_retained": bool(h358),
                "h358_filtered_site": h358[0]["filtered_codon_site"] if h358 else "",
                "h358_occupancy": h358[0]["post_filter_occupancy"] if h358 else "",
                "filtered_alignment_sha256": sha256(filtered_path),
                "deduplicated_alignment_sha256": sha256(deduplicated_path),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(RESULTS / "selection_337_input_summary.tsv", sep="\t", index=False)
    pd.DataFrame(all_crosswalks).to_csv(
        RESULTS / "selection_337_human_coordinate_crosswalk.tsv", sep="\t", index=False
    )
    pd.DataFrame(all_duplicates).to_csv(
        RESULTS / "selection_337_exact_duplicate_manifest.tsv", sep="\t", index=False
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
