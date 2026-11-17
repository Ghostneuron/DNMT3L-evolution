#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "analysis/codon_alignments_337/Placentalia"
INPUT_CDS = PROJECT / "analysis/codon_alignments_337/inputs/Placentalia_cds.fasta"
SITE_CROSSWALK = PROJECT / "results/functional_reference_coordinate_crosswalk.tsv"
RESULTS = PROJECT / "results"
ALIGNMENTS = {
    "MACSE": ROOT / "MACSE/Placentalia_MACSE_NT.fasta",
    "MAFFT": ROOT / "MAFFT/Placentalia_MAFFT_NT.fasta",
    "PRANK": ROOT / "PRANK/Placentalia_PRANK.best.fas",
    "MUSCLE5": ROOT / "MUSCLE5/Placentalia_MUSCLE5_NT.fasta",
    "COBALT": ROOT / "COBALT/Placentalia_COBALT_NT.fasta",
}
HUMAN = "Homo_sapiens"


def read_records(path: Path) -> dict[str, str]:
    records = {record.id: str(record.seq).upper().replace("U", "T") for record in SeqIO.parse(path, "fasta")}
    if not records:
        raise ValueError(f"No records in {path}")
    if len(records) != sum(1 for _ in SeqIO.parse(path, "fasta")):
        raise ValueError(f"Duplicate identifiers in {path}")
    return records


def resolved_residue(codon: str) -> str | None:
    if len(codon) != 3 or set(codon) - set("ACGT"):
        return None
    residue = str(Seq(codon).translate())
    return None if residue == "*" else residue


def allowed_source_contents(sequence: str) -> set[str]:
    allowed = {sequence}
    if len(sequence) >= 3 and str(Seq(sequence[-3:]).translate()) == "*":
        allowed.add(sequence[:-3])
    return allowed


def audit_alignment(
    method: str, path: Path, expected_cds: dict[str, str]
) -> tuple[dict[str, object], dict[int, int], dict[int, dict[str, str | None]]]:
    records = read_records(path)
    if set(records) != set(expected_cds):
        missing = sorted(set(expected_cds) - set(records))
        extra = sorted(set(records) - set(expected_cds))
        raise ValueError(f"{method} identifier mismatch; missing={missing[:5]}, extra={extra[:5]}")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        raise ValueError(f"{method} is not rectangular")
    nt_length = next(iter(lengths))
    if nt_length % 3:
        raise ValueError(f"{method} length is not divisible by three")

    exact_content = 0
    partial_gap_codons = 0
    internal_stop_codons = 0
    terminal_stop_codons = 0
    frameshift_symbols = 0
    records_with_frameshift_symbols = 0
    for name, sequence in records.items():
        observed_codons = [
            sequence[start : start + 3]
            for start in range(0, nt_length, 3)
            if sequence[start : start + 3] != "---"
        ]
        raw_observed = "".join(observed_codons)
        record_frameshifts = raw_observed.count("!")
        frameshift_symbols += record_frameshifts
        records_with_frameshift_symbols += int(record_frameshifts > 0)
        observed = raw_observed.replace("!", "")
        if observed in allowed_source_contents(expected_cds[name]):
            exact_content += 1
        for observed_index, codon in enumerate(observed_codons):
            if codon != "---" and "-" in codon:
                partial_gap_codons += 1
            residue = resolved_residue(codon)
            if residue is None and set(codon) <= set("ACGT") and codon != "---":
                if str(Seq(codon).translate()) == "*":
                    if observed_index == len(observed_codons) - 1:
                        terminal_stop_codons += 1
                    else:
                        internal_stop_codons += 1
    if exact_content != len(records):
        raise ValueError(f"{method} preserved exact CDS content for {exact_content}/{len(records)} records")
    if partial_gap_codons or internal_stop_codons:
        raise ValueError(
            f"{method} has {partial_gap_codons} partial-gap and "
            f"{internal_stop_codons} internal-stop codons"
        )

    human = records[HUMAN]
    coordinate_to_column: dict[int, int] = {}
    position = 0
    for column in range(nt_length // 3):
        codon = human[column * 3 : (column + 1) * 3]
        if codon == "---":
            continue
        if set(codon) <= set("ACGT") and str(Seq(codon).translate()) == "*":
            continue
        if resolved_residue(codon) is None:
            raise ValueError(f"{method} human sequence has unresolved codon at column {column + 1}")
        position += 1
        coordinate_to_column[position] = column
    if position != 387:
        raise ValueError(f"{method} maps {position}, not 387, human residues")

    calls: dict[int, dict[str, str | None]] = {}
    for human_position, column in coordinate_to_column.items():
        calls[human_position] = {
            name: resolved_residue(sequence[column * 3 : (column + 1) * 3])
            for name, sequence in records.items()
        }
    h358_column = coordinate_to_column[358]
    h358_residue = calls[358][HUMAN]
    summary = {
        "method": method,
        "sequence_count": len(records),
        "alignment_codon_columns": nt_length // 3,
        "human_mapped_residues": position,
        "human_h358_alignment_codon_column": h358_column + 1,
        "human_h358_residue": h358_residue,
        "exact_ungapped_CDS_content_records": exact_content,
        "partial_gap_codons": partial_gap_codons,
        "internal_stop_codons": internal_stop_codons,
        "terminal_stop_codons": terminal_stop_codons,
        "frameshift_symbols": frameshift_symbols,
        "records_with_frameshift_symbols": records_with_frameshift_symbols,
        "status": "PASS_RAW_REQUIRES_MASKING" if frameshift_symbols else "PASS",
    }
    return summary, coordinate_to_column, calls


def main() -> None:
    expected_cds = read_records(INPUT_CDS)
    summaries: list[dict[str, object]] = []
    columns_by_method: dict[str, dict[int, int]] = {}
    calls_by_method: dict[str, dict[int, dict[str, str | None]]] = {}
    for method, path in ALIGNMENTS.items():
        summary, columns, calls = audit_alignment(method, path, expected_cds)
        summaries.append(summary)
        columns_by_method[method] = columns
        calls_by_method[method] = calls

    pd.DataFrame(summaries).to_csv(
        RESULTS / "placental_five_aligner_summary.tsv", sep="\t", index=False
    )
    methods = list(ALIGNMENTS)
    coordinate_rows: list[dict[str, object]] = []
    for position in range(1, 388):
        row: dict[str, object] = {"human_NP_037501_2_position": position}
        for method in methods:
            calls = calls_by_method[method][position]
            row[f"{method}_alignment_codon_column"] = columns_by_method[method][position] + 1
            row[f"{method}_occupancy"] = sum(value is not None for value in calls.values()) / len(calls)
        resolved_all = 0
        concordant_all = 0
        discordant_species: list[str] = []
        for species in sorted(expected_cds):
            values = [calls_by_method[method][position][species] for method in methods]
            if all(value is not None for value in values):
                resolved_all += 1
                if len(set(values)) == 1:
                    concordant_all += 1
                else:
                    discordant_species.append(species)
        row["species_resolved_in_all_five"] = resolved_all
        row["species_concordant_in_all_five"] = concordant_all
        row["all_five_concordance_fraction"] = (
            concordant_all / resolved_all if resolved_all else ""
        )
        row["discordant_species"] = ";".join(discordant_species)
        row["minimum_method_occupancy"] = min(float(row[f"{method}_occupancy"]) for method in methods)
        row["passes_80pct_occupancy_in_all_five"] = row["minimum_method_occupancy"] >= 0.8
        coordinate_rows.append(row)
    coordinate_audit = pd.DataFrame(coordinate_rows)
    coordinate_audit.to_csv(
        RESULTS / "placental_five_aligner_human_coordinate_audit.tsv", sep="\t", index=False
    )

    sites = pd.read_csv(SITE_CROSSWALK, sep="\t").loc[lambda frame: frame["status"].eq("PASS")]
    site_rows: list[dict[str, object]] = []
    for site in sites.itertuples(index=False):
        position = int(site.mapped_human_position)
        coordinate = coordinate_audit.loc[
            coordinate_audit["human_NP_037501_2_position"].eq(position)
        ].iloc[0]
        site_rows.append(
            {
                "feature_id": site.feature_id,
                "human_position": position,
                "human_residue": site.mapped_human_residue,
                "minimum_method_occupancy": coordinate["minimum_method_occupancy"],
                "species_resolved_in_all_five": coordinate["species_resolved_in_all_five"],
                "species_concordant_in_all_five": coordinate["species_concordant_in_all_five"],
                "all_five_concordance_fraction": coordinate["all_five_concordance_fraction"],
                "discordant_species": coordinate["discordant_species"],
                "passes_80pct_occupancy_in_all_five": coordinate["passes_80pct_occupancy_in_all_five"],
            }
        )
    pd.DataFrame(site_rows).to_csv(
        RESULTS / "placental_functional_sites_five_aligner_audit.tsv", sep="\t", index=False
    )
    print(pd.DataFrame(summaries).to_string(index=False))
    print("\nFunctional-site five-aligner audit")
    print(pd.DataFrame(site_rows).to_string(index=False))


if __name__ == "__main__":
    main()
