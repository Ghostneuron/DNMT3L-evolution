#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

from Bio.Seq import Seq


PROJECT = Path(__file__).resolve().parents[1]
HITS = PROJECT / "results" / "marker_centered_diamond_hits_full_family.tsv"
MANIFEST = PROJECT / "metadata" / "marker_centered_interval_manifest.tsv"
INTERVALS = PROJECT / "data" / "interim" / "marker_centered_intervals.fasta"
FAMILY_FASTA = (
    PROJECT
    / "analysis"
    / "gene_family_tree_337"
    / "DNMT3_family_337_candidates_plus_paralog_references.fasta"
)
OUTDIR = PROJECT / "analysis" / "locus_history_337"
DIAMOND_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qframe",
)
WINDOW_FLANK = 100_000


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = ""
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            else:
                records[current].append(line)
    return {name: "".join(parts) for name, parts in records.items()}


def write_fasta(path: Path, name: str, sequence: str, width: int = 80) -> None:
    lines = [f">{name}"]
    lines.extend(sequence[index : index + width] for index in range(0, len(sequence), width))
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with HITS.open() as handle:
        data_lines = (line for line in handle if line.strip() and not line.startswith("#"))
        hits = list(csv.DictReader(data_lines, fieldnames=DIAMOND_FIELDS, delimiter="\t"))
    lungfish_hits = [row for row in hits if row["qseqid"].startswith("LUNGFISH|")]
    if not lungfish_hits:
        raise ValueError("No lungfish translated hit was found")
    best = max(lungfish_hits, key=lambda row: float(row["bitscore"]))
    intervals = read_fasta(INTERVALS)
    interval_sequence = intervals[best["qseqid"]]
    with MANIFEST.open(newline="") as handle:
        manifest = {
            row["record_id"]: row for row in csv.DictReader(handle, delimiter="\t")
        }
    source = manifest[best["qseqid"]]

    qstart = int(best["qstart"])
    qend = int(best["qend"])
    left = min(qstart, qend) - 1
    right = max(qstart, qend)
    hit_nt = interval_sequence[left:right]
    if int(best["qframe"]) < 0:
        hit_nt = str(Seq(hit_nt).reverse_complement())
    peptide = str(Seq(hit_nt).translate())

    window_left = max(0, left - WINDOW_FLANK)
    window_right = min(len(interval_sequence), right + WINDOW_FLANK)
    window_nt = interval_sequence[window_left:window_right]
    genomic_start = int(source["fetch_start"]) + window_left
    genomic_end = int(source["fetch_start"]) + window_right - 1
    hit_genomic_start = int(source["fetch_start"]) + left
    hit_genomic_end = int(source["fetch_start"]) + right - 1

    family = read_fasta(FAMILY_FASTA)
    subject = best["sseqid"]
    if subject not in family:
        raise ValueError(f"Subject {subject} not found in family FASTA")

    write_fasta(
        OUTDIR / "lungfish_candidate_hit_nt.fasta",
        f"LUNGFISH_candidate_hit|{source['genomic_accession']}:{hit_genomic_start}-{hit_genomic_end}",
        hit_nt,
    )
    write_fasta(
        OUTDIR / "lungfish_candidate_hit_peptide.fasta",
        "LUNGFISH_candidate_translated_segment",
        peptide,
    )
    write_fasta(
        OUTDIR / "lungfish_candidate_target_window.fasta",
        f"LUNGFISH_target_window|{source['genomic_accession']}:{genomic_start}-{genomic_end}",
        window_nt,
    )
    write_fasta(
        OUTDIR / "geotrypetes_dnmt3l_query.fasta",
        subject,
        family[subject],
    )
    summary = {
        "panel_id": "LUNGFISH",
        "species": source["species"],
        "genomic_accession": source["genomic_accession"],
        "hit_genomic_start": hit_genomic_start,
        "hit_genomic_end": hit_genomic_end,
        "query_frame": best["qframe"],
        "translated_length_aa": len(peptide),
        "top_subject": subject,
        "percent_identity": best["pident"],
        "alignment_length_aa": best["length"],
        "evalue": best["evalue"],
        "bitscore": best["bitscore"],
        "stop_count_in_direct_translation": peptide.count("*"),
        "stop_position_aa": peptide.find("*") + 1 if "*" in peptide else "",
        "stop_codon": (
            hit_nt[(peptide.find("*")) * 3 : (peptide.find("*")) * 3 + 3]
            if "*" in peptide
            else ""
        ),
        "interpretation": (
            "single_in_frame_stop_unresolved_without_transcript_or_read_validation"
            if peptide.count("*") == 1
            else "translated_candidate_requires_manual_review"
        ),
        "target_window_start": genomic_start,
        "target_window_end": genomic_end,
    }
    with (OUTDIR / "lungfish_candidate_extraction_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary), delimiter="\t")
        writer.writeheader()
        writer.writerow(summary)
    print(summary)


if __name__ == "__main__":
    main()
