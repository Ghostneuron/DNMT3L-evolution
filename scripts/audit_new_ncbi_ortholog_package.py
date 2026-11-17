#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


PROJECT = Path(__file__).resolve().parents[1]
PACKAGE = PROJECT / "data" / "raw" / "ncbi_dataset_DNMT3L_orthologs_20261102" / "ncbi_dataset" / "data"
OLD_MATRIX = (
    PROJECT
    / "../../../source_package/analysis/orthology_rebuild_20260826/orthology_audit/DNMT3L_orthology_evidence_matrix.csv"
).resolve()
RESULTS = PROJECT / "results"
GENE_ID = re.compile(r"\[GeneID=(\d+)\]")
ORGANISM = re.compile(r"\[organism=([^\]]+)\]")


def metadata(description: str) -> tuple[str, str]:
    gene_match = GENE_ID.search(description)
    organism_match = ORGANISM.search(description)
    if not gene_match or not organism_match:
        raise ValueError(f"Could not parse FASTA metadata: {description}")
    return gene_match.group(1), organism_match.group(1)


def read_fasta(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for record in SeqIO.parse(path, "fasta"):
        gene_id, species = metadata(record.description)
        rows.append(
            {
                "accession": record.id,
                "gene_id": gene_id,
                "species": species,
                "sequence": str(record.seq).upper(),
                "description": record.description,
            }
        )
    return rows


def translate_cds(sequence: str) -> str | None:
    if len(sequence) % 3:
        return None
    protein = str(Seq(sequence).translate())
    if protein.endswith("*"):
        protein = protein[:-1]
    if "*" in protein:
        return None
    return protein


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    table = pd.read_csv(PACKAGE / "gene_data_table.tsv", sep="\t", dtype=str).fillna("")
    # NCBI's paginated web export repeats the header every 25 records.
    table = table.loc[table["NCBI GeneID"].str.fullmatch(r"\d+")].copy()
    old = pd.read_csv(OLD_MATRIX, dtype={"gene_id": str}).fillna("")
    proteins = read_fasta(PACKAGE / "protein.faa")
    cdss = read_fasta(PACKAGE / "cds.fna")

    protein_by_accession = {str(row["accession"]): row for row in proteins}
    cds_by_accession = {str(row["accession"]): row for row in cdss}
    proteins_by_gene: dict[str, list[dict[str, object]]] = defaultdict(list)
    cdss_by_gene: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in proteins:
        proteins_by_gene[str(row["gene_id"])].append(row)
    for row in cdss:
        cdss_by_gene[str(row["gene_id"])].append(row)

    unique_table = table.drop_duplicates(subset=["NCBI GeneID", "Taxonomic Name"]).copy()
    new_species = set(unique_table["Taxonomic Name"])
    old_species = set(old["species"])
    membership_rows: list[dict[str, object]] = []
    table_by_species = {
        species: group for species, group in unique_table.groupby("Taxonomic Name", sort=True)
    }
    old_by_species = {species: group for species, group in old.groupby("species", sort=True)}
    for species in sorted(new_species | old_species):
        new_group = table_by_species.get(species)
        old_group = old_by_species.get(species)
        if species in new_species and species in old_species:
            status = "SHARED"
        elif species in new_species:
            status = "NEW_ONLY"
        else:
            status = "OLD_ONLY"
        membership_rows.append(
            {
                "species": species,
                "membership_status": status,
                "new_gene_ids": (
                    ";".join(sorted(set(new_group["NCBI GeneID"]), key=int))
                    if new_group is not None
                    else ""
                ),
                "old_gene_ids": (
                    ";".join(sorted(set(old_group["gene_id"]), key=int))
                    if old_group is not None
                    else ""
                ),
                "new_gene_id_count": new_group["NCBI GeneID"].nunique() if new_group is not None else 0,
                "old_record_count": len(old_group) if old_group is not None else 0,
                "old_selection_ready_count": (
                    int(old_group["Selection_Ready"].astype(bool).sum()) if old_group is not None else 0
                ),
            }
        )
    membership = pd.DataFrame(membership_rows)
    membership.to_csv(RESULTS / "ncbi_337_species_membership_diff.tsv", sep="\t", index=False)

    recovery_rows: list[dict[str, object]] = []
    for record in old.to_dict(orient="records"):
        protein = protein_by_accession.get(record["protein_accession"])
        cds = cds_by_accession.get(record["cds_accession"])
        recovery_rows.append(
            {
                "safe_id": record["safe_id"],
                "species": record["species"],
                "gene_id": record["gene_id"],
                "old_protein_accession": record["protein_accession"],
                "old_protein_recovered": protein is not None,
                "old_protein_sequence_identical": (
                    protein is not None and str(protein["sequence"]) == record["protein_seq"]
                ),
                "old_cds_accession": record["cds_accession"],
                "old_cds_recovered": cds is not None,
                "old_cds_sequence_identical": (
                    cds is not None and str(cds["sequence"]) == record["cds_seq"]
                ),
                "old_selection_ready": bool(record["Selection_Ready"]),
            }
        )
    recovery = pd.DataFrame(recovery_rows)
    recovery.to_csv(RESULTS / "ncbi_337_old_accession_recovery.tsv", sep="\t", index=False)

    pairing_rows: list[dict[str, object]] = []
    for gene_id, gene_group in unique_table.groupby("NCBI GeneID", sort=True):
        species_values = sorted(set(gene_group["Taxonomic Name"]))
        gene_proteins = proteins_by_gene.get(str(gene_id), [])
        gene_cdss = cdss_by_gene.get(str(gene_id), [])
        protein_sequences = defaultdict(list)
        for protein in gene_proteins:
            protein_sequences[str(protein["sequence"])].append(str(protein["accession"]))
        exact_pairs: list[str] = []
        valid_translation_count = 0
        for cds in gene_cdss:
            translated = translate_cds(str(cds["sequence"]))
            if translated is None:
                continue
            valid_translation_count += 1
            for protein_accession in protein_sequences.get(translated, []):
                exact_pairs.append(f"{protein_accession}|{cds['accession']}")
        pairing_rows.append(
            {
                "gene_id": gene_id,
                "species": ";".join(species_values),
                "protein_records": len(gene_proteins),
                "unique_protein_sequences": len(protein_sequences),
                "cds_records": len(gene_cdss),
                "valid_complete_translations": valid_translation_count,
                "exact_protein_cds_pairs": len(exact_pairs),
                "exact_pair_examples": ";".join(exact_pairs[:5]),
                "pairing_status": "PASS" if exact_pairs else "NO_EXACT_PAIR",
            }
        )
    pairing = pd.DataFrame(pairing_rows)
    pairing.to_csv(RESULTS / "ncbi_337_protein_cds_pairing_audit.tsv", sep="\t", index=False)

    membership_counts = Counter(membership["membership_status"])
    summary_rows = [
        ("archive_sha256", "f4b9e71a64cb4287fc13cf7182988ffef49464164f9a0890dde511db1fec1f12"),
        ("gene_table_rows", len(table)),
        ("unique_gene_ids", unique_table["NCBI GeneID"].nunique()),
        ("unique_species", len(new_species)),
        ("protein_records", len(proteins)),
        ("cds_records", len(cdss)),
        ("shared_species", membership_counts["SHARED"]),
        ("new_only_species", membership_counts["NEW_ONLY"]),
        ("old_only_species", membership_counts["OLD_ONLY"]),
        ("old_protein_accessions_recovered", int(recovery["old_protein_recovered"].sum())),
        ("old_protein_sequences_identical", int(recovery["old_protein_sequence_identical"].sum())),
        ("old_cds_accessions_recovered", int(recovery["old_cds_recovered"].sum())),
        ("old_cds_sequences_identical", int(recovery["old_cds_sequence_identical"].sum())),
        ("genes_with_exact_protein_cds_pair", int(pairing["pairing_status"].eq("PASS").sum())),
        ("genes_without_exact_protein_cds_pair", int(pairing["pairing_status"].ne("PASS").sum())),
    ]
    with (RESULTS / "ncbi_337_package_audit_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(summary_rows)

    print("\n".join(f"{metric}: {value}" for metric, value in summary_rows))
    print("\nNew-only species")
    print(membership.loc[membership["membership_status"].eq("NEW_ONLY")].to_string(index=False))
    print("\nOld-only species")
    print(membership.loc[membership["membership_status"].eq("OLD_ONLY")].to_string(index=False))


if __name__ == "__main__":
    main()

