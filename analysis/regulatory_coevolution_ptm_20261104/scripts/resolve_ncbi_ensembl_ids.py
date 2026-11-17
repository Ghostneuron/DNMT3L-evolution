#!/usr/bin/env python3
"""Resolve NCBI Gene and Ensembl identifiers for the nine germline-atlas species."""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.error
import urllib.request
from pathlib import Path


SPECIES = [
    ("Homo sapiens", 9606),
    ("Pan troglodytes", 9598),
    ("Gorilla gorilla", 9595),
    ("Macaca mulatta", 9544),
    ("Callithrix jacchus", 9483),
    ("Mus musculus", 10090),
    ("Monodelphis domestica", 13616),
    ("Ornithorhynchus anatinus", 9258),
    ("Gallus gallus", 9031),
]


def fetch(symbol: str, taxid: int):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene/symbol/{symbol}/taxon/{taxid}"
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "DNMT3L-comparative-audit/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return url, json.load(response), "PASS"
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return url, {}, "NO_SYMBOL_RECORD"
        return url, {}, f"HTTP_{exc.code}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for species, taxid in SPECIES:
        for symbol in ("DNMT3L", "DNMT3A"):
            url, payload, query_status = fetch(symbol, taxid)
            reports = payload.get("reports", [])
            genes = [report.get("gene", {}) for report in reports]
            if query_status == "PASS" and not genes:
                query_status = "NO_RECORD"
            exact = [gene for gene in genes if gene.get("symbol", "").upper() == symbol]
            gene = exact[0] if exact else (genes[0] if len(genes) == 1 else {})
            rows.append({
                "species": species,
                "taxid": taxid,
                "query_symbol": symbol,
                "query_status": query_status,
                "returned_symbol": gene.get("symbol", ""),
                "ncbi_gene_id": gene.get("gene_id", ""),
                "ensembl_gene_ids": ";".join(gene.get("ensembl_gene_ids", [])),
                "description": gene.get("description", ""),
                "record_taxname": gene.get("taxname", ""),
                "api_url": url,
            })
            time.sleep(0.15)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".partial")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
