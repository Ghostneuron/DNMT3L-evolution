from __future__ import annotations

import csv
import json
import zipfile
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
ARCHIVE = PROJECT.parents[3] / "ncbi_dataset_DNMT3L_orthologs.zip"
MATRIX = PROJECT / "results" / "DNMT3L_337_orthology_evidence_matrix.tsv"
ALIGNMENT = (
    PROJECT
    / "analysis/selection_reanalysis_337/inputs/Placentalia_CONSENSUS_human_coordinate_occupancy80_deduplicated.fasta"
)
ROOT = PROJECT / "analysis" / "selection_reanalysis_337" / "taxonomy"
MAP_OUT = PROJECT / "metadata" / "selection_337_taxid_map.tsv"
IDS_OUT = ROOT / "selection_337_taxids.txt"


def fasta_names(path: Path) -> list[str]:
    return [line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")]


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(ARCHIVE) as archive:
        lines = archive.read("ncbi_dataset/data/data_report.jsonl").decode().splitlines()
    tax_by_gene: dict[str, tuple[str, str]] = {}
    for line in lines:
        record = json.loads(line)
        tax_by_gene[str(record["geneId"])] = (str(record["taxId"]), record["taxname"])

    with MATRIX.open(newline="") as handle:
        matrix = {row["safe_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    rows: list[dict[str, str]] = []
    for safe_id in fasta_names(ALIGNMENT):
        gene_id = matrix[safe_id]["gene_id"]
        taxid, taxname = tax_by_gene[gene_id]
        rows.append(
            {
                "safe_id": safe_id,
                "gene_id": gene_id,
                "taxid": taxid,
                "taxname": taxname,
            }
        )
    if len(rows) != len({row["safe_id"] for row in rows}):
        raise ValueError("Duplicate safe IDs")
    if len(rows) != len({row["taxid"] for row in rows}):
        raise ValueError("The deduplicated consensus input has repeated taxonomy IDs")

    with MAP_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    IDS_OUT.write_text(",".join(row["taxid"] for row in rows) + "\n")
    print(f"Prepared {len(rows)} unique taxonomy IDs")
    print(MAP_OUT.relative_to(PROJECT))
    print(IDS_OUT.relative_to(PROJECT))


if __name__ == "__main__":
    main()
