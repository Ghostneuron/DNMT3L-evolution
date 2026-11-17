#!/usr/bin/env python3
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "analysis/codon_alignments_337/Placentalia"
OUTPUT = ROOT / "five_aligner_consensus_filtered"
RESULTS = PROJECT / "results"
ALIGNMENTS = {
    "MACSE": ROOT / "MACSE/Placentalia_MACSE_NT.fasta",
    "MAFFT": ROOT / "MAFFT/Placentalia_MAFFT_NT.fasta",
    "PRANK": ROOT / "PRANK/Placentalia_PRANK.best.fas",
    "MUSCLE5": ROOT / "MUSCLE5/Placentalia_MUSCLE5_NT.fasta",
    "COBALT": ROOT / "COBALT/Placentalia_COBALT_NT.fasta",
}
HUMAN = "Homo_sapiens"
MIN_OCCUPANCY = 0.8


def read_fasta(path: Path) -> dict[str, str]:
    return {record.id: str(record.seq).upper().replace("U", "T") for record in SeqIO.parse(path, "fasta")}


def resolved_codon(codon: str) -> bool:
    return len(codon) == 3 and not (set(codon) - set("ACGT")) and str(Seq(codon).translate()) != "*"


def human_coordinate_columns(records: dict[str, str]) -> dict[int, int]:
    sequence = records[HUMAN]
    mapping: dict[int, int] = {}
    position = 0
    for column in range(len(sequence) // 3):
        codon = sequence[column * 3 : (column + 1) * 3]
        if codon == "---":
            continue
        if set(codon) <= set("ACGT") and str(Seq(codon).translate()) == "*":
            continue
        if not resolved_codon(codon):
            raise ValueError(f"Unresolved human codon at alignment column {column + 1}")
        position += 1
        mapping[position] = column
    if position != 387:
        raise ValueError(f"Mapped {position}, not 387, human residues")
    return mapping


def write_fasta(records: dict[str, str], path: Path) -> None:
    with path.open("w") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    records_by_method = {method: read_fasta(path) for method, path in ALIGNMENTS.items()}
    taxa_sets = {method: set(records) for method, records in records_by_method.items()}
    if len({frozenset(taxa) for taxa in taxa_sets.values()}) != 1:
        raise ValueError("Five aligners do not contain the same taxon set")
    taxa = sorted(next(iter(taxa_sets.values())))
    columns_by_method = {
        method: human_coordinate_columns(records)
        for method, records in records_by_method.items()
    }

    site_rows: list[dict[str, object]] = []
    retained_positions: list[int] = []
    codons: dict[str, dict[int, dict[str, str]]] = {
        method: {} for method in ALIGNMENTS
    }
    for method, records in records_by_method.items():
        for position, column in columns_by_method[method].items():
            codons[method][position] = {
                taxon: sequence[column * 3 : (column + 1) * 3]
                for taxon, sequence in records.items()
            }

    for position in range(1, 388):
        occupancies = {
            method: sum(resolved_codon(value) for value in codons[method][position].values()) / len(taxa)
            for method in ALIGNMENTS
        }
        keep = min(occupancies.values()) >= MIN_OCCUPANCY
        if keep:
            retained_positions.append(position)
        resolved_all = 0
        amino_acid_concordant = 0
        exact_codon_concordant = 0
        discordant_taxa: list[str] = []
        for taxon in taxa:
            values = [codons[method][position][taxon] for method in ALIGNMENTS]
            if all(resolved_codon(value) for value in values):
                resolved_all += 1
                residues = [str(Seq(value).translate()) for value in values]
                amino_acid_concordant += int(len(set(residues)) == 1)
                exact_codon_concordant += int(len(set(values)) == 1)
                if len(set(values)) != 1:
                    discordant_taxa.append(taxon)
            else:
                discordant_taxa.append(taxon)
        site_rows.append(
            {
                "human_position": position,
                **{f"{method}_occupancy": occupancy for method, occupancy in occupancies.items()},
                "minimum_method_occupancy": min(occupancies.values()),
                "retained_by_80pct_all_method_gate": keep,
                "taxa_resolved_in_all_five": resolved_all,
                "taxa_amino_acid_concordant_in_all_five": amino_acid_concordant,
                "taxa_exact_codon_concordant_in_all_five": exact_codon_concordant,
                "taxa_masked_by_strict_codon_consensus": len(taxa) - exact_codon_concordant,
                "masked_taxa": ";".join(discordant_taxa),
            }
        )

    filtered_by_method: dict[str, dict[str, str]] = {method: {} for method in ALIGNMENTS}
    masked_calls = 0
    total_retained_calls = len(retained_positions) * len(taxa)
    for taxon in taxa:
        output_codons = {method: [] for method in ALIGNMENTS}
        for position in retained_positions:
            values = [codons[method][position][taxon] for method in ALIGNMENTS]
            concordant = all(resolved_codon(value) for value in values) and len(set(values)) == 1
            if not concordant:
                masked_calls += 1
            for method, value in zip(ALIGNMENTS, values):
                output_codons[method].append(value if concordant else "NNN")
        for method in ALIGNMENTS:
            filtered_by_method[method][taxon] = "".join(output_codons[method])

    filtered_paths = []
    for method, records in filtered_by_method.items():
        path = OUTPUT / f"Placentalia_{method}_five_aligner_consensus_filtered.fasta"
        write_fasta(records, path)
        filtered_paths.append(path)
    filtered_hashes = {sha256(path) for path in filtered_paths}

    site_table = pd.DataFrame(site_rows)
    site_table.to_csv(RESULTS / "placental_five_aligner_consensus_mask_by_site.tsv", sep="\t", index=False)
    summary = pd.DataFrame(
        [
            ("input_human_positions", 387),
            ("retained_human_positions", len(retained_positions)),
            ("removed_low_occupancy_positions", 387 - len(retained_positions)),
            ("taxa", len(taxa)),
            ("retained_taxon_position_calls", total_retained_calls),
            ("taxon_position_calls_masked_in_all_methods", masked_calls),
            ("masked_fraction", masked_calls / total_retained_calls),
            ("h358_retained", 358 in retained_positions),
            ("filtered_alignments_byte_identical", len(filtered_hashes) == 1),
            ("consensus_alignment_sha256", next(iter(filtered_hashes))),
        ],
        columns=["metric", "value"],
    )
    summary.to_csv(RESULTS / "placental_five_aligner_consensus_mask_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print("\nPrespecified functional sites")
    functional = pd.read_csv(RESULTS / "functional_reference_coordinate_crosswalk.tsv", sep="\t")
    positions = set(functional.loc[functional["status"].eq("PASS"), "mapped_human_position"].astype(int))
    print(site_table.loc[site_table["human_position"].isin(positions)].to_string(index=False))


if __name__ == "__main__":
    main()
