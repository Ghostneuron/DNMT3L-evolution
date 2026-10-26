from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd


REPRESENTATIVES = (
    "Homo_sapiens",
    "Chelonia_mydas",
    "Alligator_sinensis",
    "Crocodylus_porosus",
    "Pogona_vitticeps",
    "Lacerta_agilis",
    "Ahaetulla_prasina",
    "Eublepharis_macularius",
)


def request_json(url: str, cache_path: Path) -> dict:
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "DNMT3L-revision-synteny/20261025"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = response.read()
    cache_path.write_bytes(payload)
    time.sleep(0.4)
    return json.loads(payload)


def gene_ids_for_region(species: str, chromosome: str, start: int, end: int, cache_path: Path) -> list[str]:
    terms = [f'"{species}"[Organism]']
    if chromosome and chromosome != "Unplaced Scaffold":
        terms.append(f'"{chromosome}"[Chromosome]')
    terms.append(f"{start}:{end}[ChrPos]")
    query = " AND ".join(terms)
    params = urllib.parse.urlencode(
        {"db": "gene", "term": query, "retmode": "json", "retmax": 2000}
    )
    payload = request_json(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + params,
        cache_path,
    )
    return [str(value) for value in payload["esearchresult"].get("idlist", [])]


def summaries_for_ids(ids: list[str], cache_dir: Path, prefix: str) -> list[dict]:
    records: list[dict] = []
    for batch_index, start in enumerate(range(0, len(ids), 100), start=1):
        batch = ids[start : start + 100]
        params = urllib.parse.urlencode(
            {"db": "gene", "id": ",".join(batch), "retmode": "json"}
        )
        payload = request_json(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + params,
            cache_dir / f"{prefix}_esummary_{batch_index:03d}.json",
        )
        result = payload.get("result", {})
        for uid in result.get("uids", []):
            if uid in result:
                records.append(result[uid])
    return records


def location_on_accession(record: dict, accession: str) -> dict | None:
    for location in record.get("genomicinfo", []):
        if str(location.get("chraccver", "")) == accession:
            return location
    return None


def normalize_symbols(record: dict) -> list[str]:
    tokens = [str(record.get("name", ""))]
    tokens.extend(part.strip() for part in str(record.get("otheraliases", "")).split(","))
    result = []
    for token in tokens:
        if not token or token.startswith("LOC"):
            continue
        if re.fullmatch(r"[A-Za-z0-9_.-]+", token):
            result.append(token.upper())
    return sorted(set(result))


def signed_distance(gene_start: int, gene_end: int, dnmt_start: int, dnmt_end: int) -> int:
    if gene_end < dnmt_start:
        return gene_end - dnmt_start
    if gene_start > dnmt_end:
        return gene_start - dnmt_end
    return 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=1_000_000)
    args = parser.parse_args()
    project = args.project_root.resolve()
    outdir = args.outdir.resolve()
    cache_dir = outdir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    matrix_path = (
        project
        / "source_package"
        / "analysis"
        / "orthology_rebuild_20260826"
        / "orthology_audit"
        / "DNMT3L_orthology_evidence_matrix.csv"
    )
    matrix = pd.read_csv(matrix_path)
    selected = matrix.loc[matrix["safe_id"].isin(REPRESENTATIVES)].copy()
    missing = sorted(set(REPRESENTATIVES) - set(selected["safe_id"]))
    if missing:
        raise ValueError(f"Missing representatives: {missing}")

    neighborhood_rows: list[dict[str, object]] = []
    search_rows: list[dict[str, object]] = []
    for row in selected.sort_values("safe_id").itertuples(index=False):
        begin = int(min(row.NCBI_Genomic_Begin, row.NCBI_Genomic_End))
        end = int(max(row.NCBI_Genomic_Begin, row.NCBI_Genomic_End))
        region_start = max(1, begin - args.flank)
        region_end = end + args.flank
        prefix = row.safe_id
        ids = gene_ids_for_region(
            str(row.species),
            str(row.NCBI_Chromosome_or_Sequence),
            region_start,
            region_end,
            cache_dir / f"{prefix}_esearch.json",
        )
        target_gene_id = str(int(row.gene_id))
        if target_gene_id not in ids:
            ids.append(target_gene_id)
        ids = sorted(set(ids), key=int)
        summaries = summaries_for_ids(ids, cache_dir, f"{prefix}_with_target")
        kept = 0
        for record in summaries:
            location = location_on_accession(record, str(row.NCBI_Genomic_Accession))
            if location is None:
                continue
            gene_start = int(min(location["chrstart"], location["chrstop"])) + 1
            gene_end = int(max(location["chrstart"], location["chrstop"])) + 1
            if gene_end < region_start or gene_start > region_end:
                continue
            kept += 1
            gene_id = str(record.get("uid", ""))
            neighborhood_rows.append(
                {
                    "safe_id": row.safe_id,
                    "species": row.species,
                    "assembly": row.NCBI_Assembly_Accession,
                    "genomic_accession": row.NCBI_Genomic_Accession,
                    "chromosome_or_sequence": row.NCBI_Chromosome_or_Sequence,
                    "DNMT3L_gene_id": str(int(row.gene_id)),
                    "DNMT3L_begin": begin,
                    "DNMT3L_end": end,
                    "neighbor_gene_id": gene_id,
                    "neighbor_name": record.get("name", ""),
                    "neighbor_description": record.get("description", ""),
                    "neighbor_aliases": record.get("otheraliases", ""),
                    "normalized_symbols": ";".join(normalize_symbols(record)),
                    "neighbor_begin": gene_start,
                    "neighbor_end": gene_end,
                    "neighbor_orientation": "minus" if int(location["chrstart"]) > int(location["chrstop"]) else "plus",
                    "distance_from_DNMT3L_bp": signed_distance(gene_start, gene_end, begin, end),
                    "is_target_DNMT3L": gene_id == str(int(row.gene_id)),
                }
            )
        search_rows.append(
            {
                "safe_id": row.safe_id,
                "species": row.species,
                "assembly": row.NCBI_Assembly_Accession,
                "genomic_accession": row.NCBI_Genomic_Accession,
                "chromosome_or_sequence": row.NCBI_Chromosome_or_Sequence,
                "region_start": region_start,
                "region_end": region_end,
                "esearch_gene_ids": len(ids),
                "same_accession_genes_retained": kept,
            }
        )

    neighborhoods = pd.DataFrame(neighborhood_rows).sort_values(
        ["safe_id", "neighbor_begin", "neighbor_end", "neighbor_gene_id"]
    )
    neighborhoods.to_csv(outdir / "representative_sauropsid_gene_neighborhoods.csv", index=False)
    pd.DataFrame(search_rows).to_csv(outdir / "representative_sauropsid_synteny_search_manifest.csv", index=False)

    human = neighborhoods.loc[neighborhoods["safe_id"] == "Homo_sapiens"]
    human_symbols: set[str] = set()
    for value in human.loc[~human["is_target_DNMT3L"], "normalized_symbols"]:
        human_symbols.update(symbol for symbol in str(value).split(";") if symbol)
    human_symbols.discard("DNMT3L")

    summary_rows = []
    for safe_id, group in neighborhoods.groupby("safe_id", sort=True):
        target = group.loc[group["is_target_DNMT3L"]]
        symbols: set[str] = set()
        for value in group.loc[~group["is_target_DNMT3L"], "normalized_symbols"]:
            symbols.update(symbol for symbol in str(value).split(";") if symbol)
        shared = sorted(symbols & human_symbols)
        nearest = (
            group.loc[~group["is_target_DNMT3L"]]
            .assign(abs_distance=lambda x: x["distance_from_DNMT3L_bp"].abs())
            .sort_values(["abs_distance", "neighbor_gene_id"])
            .head(6)
        )
        summary_rows.append(
            {
                "safe_id": safe_id,
                "species": group.iloc[0]["species"],
                "DNMT3L_target_recovered": not target.empty,
                "genes_on_same_accession_within_flank": len(group),
                "shared_named_markers_with_human": ";".join(shared),
                "shared_named_marker_count": len(shared),
                "nearest_six_genes": ";".join(
                    f"{row.neighbor_name}:{int(row.distance_from_DNMT3L_bp)}"
                    for row in nearest.itertuples(index=False)
                ),
                "neighborhood_support": (
                    "two_or_more_shared_named_markers"
                    if safe_id != "Homo_sapiens" and len(shared) >= 2
                    else "one_shared_named_marker"
                    if safe_id != "Homo_sapiens" and len(shared) == 1
                    else "no_named_marker_support"
                    if safe_id != "Homo_sapiens"
                    else "human_reference"
                ),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(outdir / "representative_sauropsid_synteny_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
