#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
PANEL = PROJECT / "metadata" / "lineage_panel.tsv"
CACHE = PROJECT / "data" / "raw" / "ncbi_eutils_neighborhood_screen"
HITS_OUT = PROJECT / "results" / "conserved_neighborhood_marker_hits.tsv"
SUMMARY_OUT = PROJECT / "results" / "conserved_neighborhood_screen.tsv"
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "dnmt3l_neighborhood_screen"
MARKERS = ("ICOSLG", "GATD3", "PWP2", "TRAPPC10")
MAX_SPAN_BP = 5_000_000


def get_json(endpoint: str, params: dict[str, str], cache_path: Path) -> dict:
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    query = dict(params)
    query.update({"retmode": "json", "tool": TOOL})
    url = f"{BASE}/{endpoint}?{urllib.parse.urlencode(query)}"
    request = urllib.request.Request(url, headers={"User-Agent": f"{TOOL}/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = json.loads(response.read().decode("utf-8"))
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = cache_path.with_suffix(cache_path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(cache_path)
    time.sleep(0.36)
    return payload


def gene_summaries(species: str, marker: str, panel_id: str) -> list[dict]:
    term = f'{marker}[Gene Name] AND "{species}"[Organism]'
    search = get_json(
        "esearch.fcgi",
        {"db": "gene", "term": term, "retmax": "20"},
        CACHE / f"{panel_id}_{marker}_esearch.json",
    )["esearchresult"]
    ids = list(search.get("idlist", []))
    if not ids:
        return []
    result = get_json(
        "esummary.fcgi",
        {"db": "gene", "id": ",".join(ids)},
        CACHE / f"{panel_id}_{marker}_esummary.json",
    )["result"]
    return [result[uid] for uid in result.get("uids", [])]


def main() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    HITS_OUT.parent.mkdir(parents=True, exist_ok=True)
    with PANEL.open(newline="") as handle:
        panel = list(csv.DictReader(handle, delimiter="\t"))

    hits: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    for row in panel:
        panel_id = row["panel_id"]
        species = row["species"]
        by_accession: dict[str, list[dict[str, object]]] = defaultdict(list)
        for marker in MARKERS:
            records = gene_summaries(species, marker, panel_id)
            for record in records:
                for location in record.get("genomicinfo", []):
                    accession = str(location.get("chraccver", ""))
                    if not accession:
                        continue
                    start = min(int(location["chrstart"]), int(location["chrstop"])) + 1
                    end = max(int(location["chrstart"]), int(location["chrstop"])) + 1
                    hit = {
                        "panel_id": panel_id,
                        "major_lineage": row["major_lineage"],
                        "species": species,
                        "queried_marker": marker,
                        "gene_id": record.get("uid", ""),
                        "reported_symbol": record.get("name", ""),
                        "description": record.get("description", ""),
                        "genomic_accession": accession,
                        "chromosome": location.get("chrloc", ""),
                        "start": start,
                        "end": end,
                        "orientation": "minus" if int(location["chrstart"]) > int(location["chrstop"]) else "plus",
                    }
                    hits.append(hit)
                    by_accession[accession].append(hit)

        candidates: list[tuple[int, int, str, list[dict[str, object]]]] = []
        for accession, accession_hits in by_accession.items():
            ordered_hits = sorted(accession_hits, key=lambda hit: int(hit["start"]))
            for left in range(len(ordered_hits)):
                for right in range(left, len(ordered_hits)):
                    window_hits = ordered_hits[left : right + 1]
                    start = min(int(hit["start"]) for hit in window_hits)
                    end = max(int(hit["end"]) for hit in window_hits)
                    span = end - start + 1
                    if span > MAX_SPAN_BP:
                        break
                    unique_markers = {str(hit["queried_marker"]) for hit in window_hits}
                    if len(unique_markers) >= 2:
                        candidates.append((len(unique_markers), -span, accession, window_hits))
        candidates.sort(reverse=True)

        if candidates:
            marker_count, negative_span, accession, best_hits = candidates[0]
            interval_start = min(int(hit["start"]) for hit in best_hits)
            interval_end = max(int(hit["end"]) for hit in best_hits)
            unique_markers = sorted({str(hit["queried_marker"]) for hit in best_hits})
            status = "COLOCALIZED_MARKERS"
        else:
            accession = ""
            interval_start = ""
            interval_end = ""
            unique_markers = []
            marker_count = 0
            negative_span = 0
            status = "INSUFFICIENT_MARKERS"

        summaries.append(
            {
                "panel_id": panel_id,
                "major_lineage": row["major_lineage"],
                "species": species,
                "queried_markers": ";".join(MARKERS),
                "markers_recovered_anywhere": ";".join(
                    sorted({str(hit["queried_marker"]) for hit in hits if hit["panel_id"] == panel_id})
                ),
                "best_genomic_accession": accession,
                "colocalized_markers": ";".join(unique_markers),
                "colocalized_marker_count": marker_count,
                "marker_interval_start": interval_start,
                "marker_interval_end": interval_end,
                "marker_span_bp": -negative_span if candidates else "",
                "neighborhood_screen_status": status,
                "locus_interpretation": (
                    "candidate conserved neighborhood; DNMT3L presence/absence unresolved"
                    if status == "COLOCALIZED_MARKERS"
                    else "marker screen insufficient; DNMT3L presence/absence unresolved"
                ),
                "retrieved_at": "2026-11-02",
            }
        )
        print(
            f"{panel_id}: {status}; colocalized={';'.join(unique_markers) or 'none'}; "
            f"accession={accession or 'none'}"
        )

    with HITS_OUT.with_suffix(HITS_OUT.suffix + ".tmp").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(hits[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(hits)
    HITS_OUT.with_suffix(HITS_OUT.suffix + ".tmp").replace(HITS_OUT)

    with SUMMARY_OUT.with_suffix(SUMMARY_OUT.suffix + ".tmp").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)
    SUMMARY_OUT.with_suffix(SUMMARY_OUT.suffix + ".tmp").replace(SUMMARY_OUT)
    print(f"Wrote {HITS_OUT.relative_to(PROJECT)} ({len(hits)} marker-location rows)")
    print(f"Wrote {SUMMARY_OUT.relative_to(PROJECT)} ({len(summaries)} species summaries)")


if __name__ == "__main__":
    main()

