#!/usr/bin/env python3
from __future__ import annotations

import csv
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SCREEN = PROJECT / "results" / "conserved_neighborhood_screen.tsv"
MARKER_HITS = PROJECT / "results" / "conserved_neighborhood_marker_hits.tsv"
OUTDIR = PROJECT / "data" / "raw" / "marker_centered_intervals"
COMBINED = PROJECT / "data" / "interim" / "marker_centered_intervals.fasta"
MANIFEST = PROJECT / "metadata" / "marker_centered_interval_manifest.tsv"
FLANK_BP = 250_000
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch_fasta(accession: str, start: int, end: int, destination: Path) -> str:
    if destination.exists():
        return destination.read_text()
    params = urllib.parse.urlencode(
        {
            "db": "nuccore",
            "id": accession,
            "seq_start": str(start),
            "seq_stop": str(end),
            "rettype": "fasta",
            "retmode": "text",
            "tool": "dnmt3l_marker_centered_fetch",
        }
    )
    request = urllib.request.Request(
        f"{BASE}?{params}",
        headers={"User-Agent": "dnmt3l_marker_centered_fetch/1.0"},
    )
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                text = response.read().decode("utf-8")
            if not text.startswith(">"):
                raise ValueError(f"Unexpected efetch response for {accession}: {text[:120]!r}")
            destination.parent.mkdir(parents=True, exist_ok=True)
            temporary = destination.with_suffix(destination.suffix + ".tmp")
            temporary.write_text(text)
            temporary.replace(destination)
            time.sleep(0.36)
            return text
        except Exception as error:
            last_error = error
            time.sleep(1.0 + attempt)
    raise RuntimeError(f"Failed to fetch {accession}:{start}-{end}") from last_error


def fasta_sequence(text: str) -> str:
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


def main() -> None:
    with SCREEN.open(newline="") as handle:
        screen_rows = list(csv.DictReader(handle, delimiter="\t"))
    unresolved = {
        row["panel_id"]: row
        for row in screen_rows
        if row["neighborhood_screen_status"] == "INSUFFICIENT_MARKERS"
    }
    with MARKER_HITS.open(newline="") as handle:
        hit_rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["panel_id"] in unresolved
        ]
    by_panel_accession: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in hit_rows:
        by_panel_accession[(row["panel_id"], row["genomic_accession"])].append(row)
    if not by_panel_accession:
        raise ValueError("No marker locations were available for insufficient-marker lineages")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    COMBINED.parent.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, object]] = []
    combined_records: list[str] = []
    for (panel_id, accession), rows in sorted(by_panel_accession.items()):
        marker_start = min(int(row["start"]) for row in rows)
        marker_end = max(int(row["end"]) for row in rows)
        start = max(1, marker_start - FLANK_BP)
        end = marker_end + FLANK_BP
        destination = OUTDIR / f"{panel_id}__{accession}__{start}_{end}.fasta"
        text = fetch_fasta(accession, start, end, destination)
        sequence = fasta_sequence(text).upper()
        expected = end - start + 1
        if len(sequence) != expected:
            raise ValueError(
                f"Length mismatch for {panel_id} {accession}: expected {expected}, recovered {len(sequence)}"
            )
        markers = sorted({row["queried_marker"] for row in rows})
        record_id = f"{panel_id}|MARKER_CENTERED|{accession}|{start}-{end}"
        combined_records.append(f">{record_id}\n{sequence}\n")
        manifest_rows.append(
            {
                "panel_id": panel_id,
                "major_lineage": unresolved[panel_id]["major_lineage"],
                "species": unresolved[panel_id]["species"],
                "genomic_accession": accession,
                "markers": ";".join(markers),
                "marker_interval_start": marker_start,
                "marker_interval_end": marker_end,
                "fetch_start": start,
                "fetch_end": end,
                "fetched_length": len(sequence),
                "record_id": record_id,
                "local_fasta": str(destination.relative_to(PROJECT)),
                "retrieved_at": "2026-11-04",
                "evidence_tier": "marker_centered_noncolocalized",
            }
        )
        print(
            f"{panel_id}: {accession}:{start}-{end} ({len(sequence):,} bp); "
            f"markers={';'.join(markers)}"
        )

    temporary_combined = COMBINED.with_suffix(COMBINED.suffix + ".tmp")
    temporary_combined.write_text("".join(combined_records))
    temporary_combined.replace(COMBINED)
    temporary_manifest = MANIFEST.with_suffix(MANIFEST.suffix + ".tmp")
    with temporary_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    temporary_manifest.replace(MANIFEST)
    print(f"Wrote {COMBINED.relative_to(PROJECT)}")
    print(f"Wrote {MANIFEST.relative_to(PROJECT)}")


if __name__ == "__main__":
    main()
