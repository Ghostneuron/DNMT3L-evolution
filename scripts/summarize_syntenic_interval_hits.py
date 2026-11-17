#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
HITS = PROJECT / "results" / "syntenic_interval_diamond_hits_expanded.tsv"
MANIFEST = PROJECT / "metadata" / "syntenic_interval_manifest.tsv"
OUTPUT = PROJECT / "results" / "syntenic_interval_dnmt3_summary_expanded.tsv"
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


def family(subject: str) -> str:
    return subject.split("__", 1)[0]


def main() -> None:
    with HITS.open(newline="") as handle:
        data_lines = (line for line in handle if line.strip() and not line.startswith("#"))
        rows = list(csv.DictReader(data_lines, fieldnames=DIAMOND_FIELDS, delimiter="\t"))
    with MANIFEST.open(newline="") as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
    manifest_by_panel = {row["panel_id"]: row for row in manifest_rows}
    by_panel: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        panel_id = row["qseqid"].split("|", 1)[0]
        row["family"] = family(row["sseqid"])
        by_panel[panel_id].append(row)

    summaries: list[dict[str, object]] = []
    for panel_id in sorted(manifest_by_panel):
        panel_rows = by_panel.get(panel_id, [])
        manifest = manifest_by_panel[panel_id]
        if not panel_rows:
            summaries.append(
                {
                    "panel_id": panel_id,
                    "species": manifest["species"],
                    "genomic_accession": manifest["genomic_accession"],
                    "fetch_start": manifest["fetch_start"],
                    "fetch_end": manifest["fetch_end"],
                    "total_hsps": 0,
                    "best_family": "",
                    "best_subject": "",
                    "best_bitscore": "",
                    "best_evalue": "",
                    "best_alignment_length_aa": "",
                    "best_query_start": "",
                    "best_query_end": "",
                    "dnmt3l_hsp_count": 0,
                    "best_dnmt3l_subject": "",
                    "best_dnmt3l_bitscore": "",
                    "best_dnmt3l_evalue": "",
                    "best_dnmt3l_alignment_length_aa": "",
                    "best_non_dnmt3l_family": "",
                    "best_non_dnmt3l_bitscore": "",
                    "dnmt3l_minus_competitor_bitscore": "",
                    "screen_interpretation": "no_DNMT3_family_match_in_representative_db_at_Evalue_1e-10",
                    "locus_status": "UNRESOLVED",
                }
            )
            continue
        best_by_family: dict[str, dict[str, str]] = {}
        for row in panel_rows:
            current = best_by_family.get(row["family"])
            if current is None or float(row["bitscore"]) > float(current["bitscore"]):
                best_by_family[row["family"]] = row
        best_overall = max(panel_rows, key=lambda item: float(item["bitscore"]))
        dnmt3l_hits = [row for row in panel_rows if row["family"] == "DNMT3L"]
        dnmt3l_best = max(dnmt3l_hits, key=lambda item: float(item["bitscore"])) if dnmt3l_hits else None
        best_competing = max(
            (row for row in panel_rows if row["family"] != "DNMT3L"),
            key=lambda item: float(item["bitscore"]),
            default=None,
        )
        dnmt3l_best_score = float(dnmt3l_best["bitscore"]) if dnmt3l_best else None
        competitor_score = float(best_competing["bitscore"]) if best_competing else None
        summaries.append(
            {
                "panel_id": panel_id,
                "species": manifest["species"],
                "genomic_accession": manifest["genomic_accession"],
                "fetch_start": manifest["fetch_start"],
                "fetch_end": manifest["fetch_end"],
                "total_hsps": len(panel_rows),
                "best_family": best_overall["family"],
                "best_subject": best_overall["sseqid"],
                "best_bitscore": best_overall["bitscore"],
                "best_evalue": best_overall["evalue"],
                "best_alignment_length_aa": best_overall["length"],
                "best_query_start": best_overall["qstart"],
                "best_query_end": best_overall["qend"],
                "dnmt3l_hsp_count": len(dnmt3l_hits),
                "best_dnmt3l_subject": dnmt3l_best["sseqid"] if dnmt3l_best else "",
                "best_dnmt3l_bitscore": dnmt3l_best["bitscore"] if dnmt3l_best else "",
                "best_dnmt3l_evalue": dnmt3l_best["evalue"] if dnmt3l_best else "",
                "best_dnmt3l_alignment_length_aa": dnmt3l_best["length"] if dnmt3l_best else "",
                "best_non_dnmt3l_family": best_competing["family"] if best_competing else "",
                "best_non_dnmt3l_bitscore": best_competing["bitscore"] if best_competing else "",
                "dnmt3l_minus_competitor_bitscore": (
                    dnmt3l_best_score - competitor_score
                    if dnmt3l_best_score is not None and competitor_score is not None
                    else ""
                ),
                "screen_interpretation": (
                    "DNMT3L_like_translated_match_detected"
                    if dnmt3l_best is not None and best_overall["family"] == "DNMT3L"
                    else "DNMT3L_match_not_best_family"
                    if dnmt3l_best is not None
                    else "no_DNMT3L_match_in_representative_db_at_Evalue_1e-10"
                ),
                "locus_status": "UNRESOLVED",
            }
        )

    temporary = OUTPUT.with_suffix(OUTPUT.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)
    temporary.replace(OUTPUT)
    print(f"Wrote {OUTPUT.relative_to(PROJECT)} ({len(summaries)} interval summaries)")
    for row in summaries:
        print(
            f"{row['panel_id']}: best={row['best_family'] or 'none'} {row['best_bitscore'] or ''}; "
            f"DNMT3L={row['best_dnmt3l_bitscore'] or 'none'}; {row['screen_interpretation']}"
        )


if __name__ == "__main__":
    main()

