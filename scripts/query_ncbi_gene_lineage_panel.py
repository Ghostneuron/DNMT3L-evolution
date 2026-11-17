#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
PANEL = PROJECT / "metadata" / "lineage_panel.tsv"
CACHE = PROJECT / "data" / "raw" / "ncbi_eutils_annotation_screen"
OUTPUT = PROJECT / "results" / "vertebrate_dnmt3l_annotation_screen.tsv"
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "dnmt3l_gene_history_redesign"


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


def esearch(db: str, term: str, cache_path: Path, retmax: int = 20) -> tuple[int, list[str]]:
    payload = get_json(
        "esearch.fcgi",
        {"db": db, "term": term, "retmax": str(retmax)},
        cache_path,
    )
    result = payload["esearchresult"]
    return int(result["count"]), list(result.get("idlist", []))


def esummary(db: str, ids: list[str], cache_path: Path) -> list[dict]:
    if not ids:
        return []
    payload = get_json(
        "esummary.fcgi",
        {"db": db, "id": ",".join(ids)},
        cache_path,
    )
    result = payload["result"]
    return [result[uid] for uid in result.get("uids", [])]


def main() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with PANEL.open(newline="") as handle:
        panel_rows = list(csv.DictReader(handle, delimiter="\t"))

    output_rows: list[dict[str, object]] = []
    for row in panel_rows:
        panel_id = row["panel_id"]
        species = row["species"]
        gene_term = f'DNMT3L[Gene Name] AND "{species}"[Organism]'
        protein_term = (
            f'(DNMT3L[Protein Name] OR "DNA (cytosine-5)-methyltransferase 3-like"[Protein Name]) '
            f'AND "{species}"[Organism]'
        )
        gene_count, gene_ids = esearch(
            "gene", gene_term, CACHE / f"{panel_id}_gene_esearch.json"
        )
        protein_count, protein_ids = esearch(
            "protein", protein_term, CACHE / f"{panel_id}_protein_esearch.json"
        )
        gene_summaries = esummary(
            "gene", gene_ids, CACHE / f"{panel_id}_gene_esummary.json"
        )
        protein_summaries = esummary(
            "protein", protein_ids, CACHE / f"{panel_id}_protein_esummary.json"
        )

        gene_symbols = sorted({str(item.get("name", "")) for item in gene_summaries if item.get("name")})
        gene_descriptions = sorted(
            {str(item.get("description", "")) for item in gene_summaries if item.get("description")}
        )
        protein_accessions = sorted(
            {
                str(item.get("accessionversion", item.get("caption", "")))
                for item in protein_summaries
                if item.get("accessionversion") or item.get("caption")
            }
        )
        status = "ANNOTATED_CANDIDATE_PRESENT" if gene_count or protein_count else "NO_ANNOTATED_CANDIDATE"
        output_rows.append(
            {
                "panel_id": panel_id,
                "major_lineage": row["major_lineage"],
                "species": species,
                "gene_query": gene_term,
                "gene_count": gene_count,
                "gene_ids": ";".join(gene_ids),
                "gene_symbols": ";".join(gene_symbols),
                "gene_descriptions": ";".join(gene_descriptions),
                "protein_query": protein_term,
                "protein_count": protein_count,
                "protein_ids": ";".join(protein_ids),
                "protein_accessions": ";".join(protein_accessions),
                "annotation_screen_status": status,
                "locus_status": "UNRESOLVED",
                "absence_interpretation": row["absence_interpretation"],
                "retrieved_at": "2026-11-02",
            }
        )
        print(f"{panel_id}: genes={gene_count}, proteins={protein_count}, status={status}")

    fieldnames = list(output_rows[0])
    temporary = OUTPUT.with_suffix(OUTPUT.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    temporary.replace(OUTPUT)
    print(f"Wrote {OUTPUT.relative_to(PROJECT)} ({len(output_rows)} lineage-panel species)")


if __name__ == "__main__":
    main()

