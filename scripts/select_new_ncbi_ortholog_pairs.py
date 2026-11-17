#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq


PROJECT = Path(__file__).resolve().parents[1]
PACKAGE = PROJECT / "data/raw/ncbi_dataset_DNMT3L_orthologs_20261102/ncbi_dataset/data"
OLD_MATRIX = (
    PROJECT
    / "../../../source_package/analysis/orthology_rebuild_20260826/orthology_audit/"
    "DNMT3L_orthology_evidence_matrix.csv"
).resolve()
MANUAL_CLADES = PROJECT / "metadata/new_gene_clade_assignments.tsv"
RESULTS = PROJECT / "results"
PROCESSED = PROJECT / "data/processed/ncbi_337_representatives"

GENE_ID_PATTERN = re.compile(r"\[GeneID=(\d+)\]")
ORGANISM_PATTERN = re.compile(r"\[organism=([^\]]+)\]")
ADD_RANGE = (41, 173)
CTD_RANGE = (178, 379)
FOCAL_POSITION = 358
MIN_FULL_REFERENCE_COVERAGE = 0.80
MIN_CTD_REFERENCE_COVERAGE = 0.90
MIN_ADD_REFERENCE_COVERAGE = 0.80


def read_fasta(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for record in SeqIO.parse(path, "fasta"):
        gene_match = GENE_ID_PATTERN.search(record.description)
        organism_match = ORGANISM_PATTERN.search(record.description)
        if gene_match is None or organism_match is None:
            raise ValueError(f"Cannot parse GeneID/organism from {record.description}")
        rows.append(
            {
                "accession": record.id,
                "gene_id": gene_match.group(1),
                "species": organism_match.group(1),
                "description": record.description,
                "sequence": str(record.seq).upper().replace("U", "T"),
            }
        )
    return rows


def translation_diagnostics(sequence: str) -> dict[str, object]:
    frame_complete = len(sequence) % 3 == 0
    ambiguous_codons = sum(
        bool(set(sequence[start : start + 3]) - set("ACGT"))
        for start in range(0, len(sequence) - 2, 3)
    )
    translated = str(Seq(sequence).translate()) if frame_complete else ""
    if translated.endswith("*"):
        translated = translated[:-1]
    return {
        "frame_complete": frame_complete,
        "translated": translated,
        "internal_stop_count": translated.count("*"),
        "ambiguous_codon_count": ambiguous_codons,
    }


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    aligner.target_end_gap_score = -2.0
    aligner.query_end_gap_score = -2.0
    return aligner


def reference_metrics(
    reference: str, query: str, aligner: PairwiseAligner
) -> dict[str, object]:
    alignment = aligner.align(reference, query)[0]
    reference_blocks, query_blocks = alignment.aligned
    covered: set[int] = set()
    identical: set[int] = set()
    focal_query_position: int | None = None
    focal_query_residue: str | None = None
    for (ref_start, ref_end), (query_start, query_end) in zip(
        reference_blocks, query_blocks
    ):
        for offset, (ref_aa, query_aa) in enumerate(
            zip(reference[ref_start:ref_end], query[query_start:query_end])
        ):
            ref_position = int(ref_start + offset + 1)
            query_position = int(query_start + offset + 1)
            covered.add(ref_position)
            if ref_aa == query_aa:
                identical.add(ref_position)
            if ref_position == FOCAL_POSITION:
                focal_query_position = query_position
                focal_query_residue = query_aa

    def bounded(bounds: tuple[int, int]) -> tuple[float, float]:
        positions = set(range(bounds[0], bounds[1] + 1))
        observed = positions & covered
        coverage = len(observed) / len(positions)
        identity = len(observed & identical) / len(observed) if observed else 0.0
        return coverage, identity

    add_coverage, add_identity = bounded(ADD_RANGE)
    ctd_coverage, ctd_identity = bounded(CTD_RANGE)
    full_positions = set(range(1, len(reference) + 1))
    full_observed = full_positions & covered
    full_coverage = len(full_observed) / len(reference)
    full_identity = (
        len(full_observed & identical) / len(full_observed) if full_observed else 0.0
    )
    return {
        "Full_Reference_Coverage": full_coverage,
        "Full_Aligned_Identity": full_identity,
        "ADD_Reference_Coverage": add_coverage,
        "ADD_Aligned_Identity": add_identity,
        "CTD_Reference_Coverage": ctd_coverage,
        "CTD_Aligned_Identity": ctd_identity,
        "Pairwise_h358_Query_Position": focal_query_position,
        "Pairwise_h358_Residue": focal_query_residue,
        "Pairwise_Alignment_Score": float(alignment.score),
    }


def is_curated_protein(accession: str) -> bool:
    return accession.startswith(("NP_", "YP_"))


def is_curated_cds(accession: str) -> bool:
    return accession.startswith(("NM_", "NR_"))


def candidate_sort_key(row: dict[str, object]) -> tuple[object, ...]:
    return (
        -int(bool(row["Sequence_QC_Pass"])),
        -int(bool(row["Domain_QC_Pass"])),
        -float(row["Full_Reference_Coverage"]),
        -float(row["CTD_Reference_Coverage"]),
        -float(row["ADD_Reference_Coverage"]),
        -float(row["Pairwise_Alignment_Score"]),
        -float(row["Full_Aligned_Identity"]),
        -int(bool(row["Curated_RefSeq_Protein"])),
        -int(bool(row["Curated_RefSeq_CDS"])),
        -int(bool(row["Protein_Starts_M"])),
        abs(int(row["protein_len"]) - 387),
        str(row["protein_accession"]),
        str(row["cds_accession"]),
    )


def safe_name(species: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", species).strip("_")


def write_fasta(rows: pd.DataFrame, sequence_column: str, path: Path) -> None:
    with path.open("w") as handle:
        for row in rows.itertuples(index=False):
            sequence = str(getattr(row, sequence_column))
            handle.write(f">{row.safe_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)

    table = pd.read_csv(PACKAGE / "gene_data_table.tsv", sep="\t", dtype=str).fillna("")
    table = table.loc[table["NCBI GeneID"].str.fullmatch(r"\d+")].copy()
    table = table.drop_duplicates(subset=["NCBI GeneID"], keep="first")
    if len(table) != 337:
        raise ValueError(f"Expected 337 unique GeneIDs, found {len(table)}")

    reports = [json.loads(line) for line in (PACKAGE / "data_report.jsonl").read_text().splitlines() if line]
    report_by_gene = {str(row["geneId"]): row for row in reports}
    old = pd.read_csv(OLD_MATRIX, dtype={"gene_id": str}).fillna("")
    old_clade_by_gene = dict(zip(old["gene_id"], old["Clade"]))
    manual = pd.read_csv(MANUAL_CLADES, sep="\t", dtype=str).fillna("")
    manual_clade_by_gene = dict(zip(manual["gene_id"], manual["clade"]))

    proteins = read_fasta(PACKAGE / "protein.faa")
    cdss = read_fasta(PACKAGE / "cds.fna")
    proteins_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    cdss_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in proteins:
        proteins_by_gene[row["gene_id"]].append(row)
    for row in cdss:
        cdss_by_gene[row["gene_id"]].append(row)

    human_reference_records = [
        row for row in proteins if row["accession"] == "NP_037501.2"
    ]
    if len(human_reference_records) != 1:
        raise ValueError("NP_037501.2 is absent or duplicated in the new archive")
    human_reference = human_reference_records[0]["sequence"]
    if len(human_reference) != 387 or human_reference[FOCAL_POSITION - 1] != "A":
        raise ValueError("Unexpected NP_037501.2 sequence")

    protein_sequence_counts = Counter(row["sequence"] for row in proteins)
    cds_sequence_counts = Counter(row["sequence"] for row in cdss)
    aligner = make_aligner()
    ranked_rows: list[dict[str, object]] = []
    selected_by_gene: dict[str, dict[str, object]] = {}

    for gene in table.to_dict(orient="records"):
        gene_id = gene["NCBI GeneID"]
        species = gene["Taxonomic Name"]
        exact_candidates: list[dict[str, object]] = []
        for cds in cdss_by_gene.get(gene_id, []):
            diagnostics = translation_diagnostics(cds["sequence"])
            for protein in proteins_by_gene.get(gene_id, []):
                if diagnostics["translated"] != protein["sequence"]:
                    continue
                metrics = reference_metrics(human_reference, protein["sequence"], aligner)
                sequence_qc_pass = (
                    bool(diagnostics["frame_complete"])
                    and int(diagnostics["internal_stop_count"]) == 0
                    and int(diagnostics["ambiguous_codon_count"]) == 0
                    and protein["sequence"].count("X") == 0
                )
                domain_qc_pass = (
                    float(metrics["Full_Reference_Coverage"])
                    >= MIN_FULL_REFERENCE_COVERAGE
                    and float(metrics["CTD_Reference_Coverage"])
                    >= MIN_CTD_REFERENCE_COVERAGE
                    and metrics["Pairwise_h358_Query_Position"] is not None
                )
                exact_candidates.append(
                    {
                        "gene_id": gene_id,
                        "species": species,
                        "protein_accession": protein["accession"],
                        "protein_description": protein["description"],
                        "protein_seq": protein["sequence"],
                        "protein_len": len(protein["sequence"]),
                        "cds_accession": cds["accession"],
                        "cds_description": cds["description"],
                        "cds_seq": cds["sequence"],
                        "cds_len": len(cds["sequence"]),
                        "CDS_Frame_Complete": diagnostics["frame_complete"],
                        "Translation_Exact_Recalculated": True,
                        "Internal_Stop_Count": diagnostics["internal_stop_count"],
                        "Ambiguous_Codon_Count": diagnostics["ambiguous_codon_count"],
                        "Protein_Ambiguous_X_Count": protein["sequence"].count("X"),
                        "Protein_Starts_M": protein["sequence"].startswith("M"),
                        "Curated_RefSeq_Protein": is_curated_protein(protein["accession"]),
                        "Curated_RefSeq_CDS": is_curated_cds(cds["accession"]),
                        "Protein_Exact_Duplicate_Count": protein_sequence_counts[protein["sequence"]],
                        "CDS_Exact_Duplicate_Count": cds_sequence_counts[cds["sequence"]],
                        **metrics,
                        "ADD_Architecture_Complete": (
                            float(metrics["ADD_Reference_Coverage"])
                            >= MIN_ADD_REFERENCE_COVERAGE
                        ),
                        "Sequence_QC_Pass": sequence_qc_pass,
                        "Domain_QC_Pass": domain_qc_pass,
                    }
                )
        exact_candidates.sort(key=candidate_sort_key)
        for rank, candidate in enumerate(exact_candidates, start=1):
            candidate["Pair_Selection_Rank"] = rank
            candidate["Representative_Selected"] = rank == 1
            ranked_rows.append(candidate)
        if exact_candidates:
            selected_by_gene[gene_id] = exact_candidates[0]

    ranked = pd.DataFrame(ranked_rows)
    ranked.to_csv(RESULTS / "ncbi_337_isoform_pair_ranking.tsv", sep="\t", index=False)

    evidence_rows: list[dict[str, object]] = []
    for gene in table.to_dict(orient="records"):
        gene_id = gene["NCBI GeneID"]
        species = gene["Taxonomic Name"]
        selected = selected_by_gene.get(gene_id, {})
        report = report_by_gene.get(gene_id, {})
        annotations = report.get("annotations", [])
        annotation = annotations[0] if annotations else {}
        locations = annotation.get("genomicLocations", [])
        location = locations[0] if locations else {}
        genomic_range = location.get("genomicRange", {})
        clade = old_clade_by_gene.get(gene_id, manual_clade_by_gene.get(gene_id, ""))
        if not clade:
            raise ValueError(f"No clade assignment for {gene_id} {species}")
        hybrid = " x " in species.lower()
        sequence_qc = bool(selected.get("Sequence_QC_Pass", False))
        domain_qc = bool(selected.get("Domain_QC_Pass", False))
        preliminary_ready = bool(selected) and sequence_qc and domain_qc and not hybrid
        evidence_rows.append(
            {
                "safe_id": safe_name(species),
                "species": species,
                "gene_id": gene_id,
                "symbol": gene["Symbol"],
                "description": gene["Description"],
                "clade": clade,
                "protein_record_count": len(proteins_by_gene.get(gene_id, [])),
                "cds_record_count": len(cdss_by_gene.get(gene_id, [])),
                "exact_pair_candidate_count": int(
                    (ranked["gene_id"] == gene_id).sum()
                ) if len(ranked) else 0,
                "protein_accession": selected.get("protein_accession", ""),
                "protein_description": selected.get("protein_description", ""),
                "protein_seq": selected.get("protein_seq", ""),
                "protein_len": selected.get("protein_len", ""),
                "cds_accession": selected.get("cds_accession", ""),
                "cds_description": selected.get("cds_description", ""),
                "cds_seq": selected.get("cds_seq", ""),
                "cds_len": selected.get("cds_len", ""),
                "NCBI_DNMT3L_Ortholog_Group": True,
                "NCBI_Assembly_Accession": annotation.get("assemblyAccession", ""),
                "NCBI_Annotation_Name": annotation.get("annotationName", ""),
                "NCBI_Annotation_Release_Date": annotation.get("annotationReleaseDate", ""),
                "NCBI_Genomic_Accession": location.get("genomicAccessionVersion", ""),
                "NCBI_Chromosome_or_Sequence": location.get("sequenceName", ""),
                "NCBI_Genomic_Begin": genomic_range.get("begin", ""),
                "NCBI_Genomic_End": genomic_range.get("end", ""),
                "NCBI_Genomic_Orientation": genomic_range.get("orientation", ""),
                "Hybrid_Taxon": hybrid,
                **{
                    key: selected.get(key, "")
                    for key in [
                        "CDS_Frame_Complete",
                        "Translation_Exact_Recalculated",
                        "Internal_Stop_Count",
                        "Ambiguous_Codon_Count",
                        "Protein_Ambiguous_X_Count",
                        "Protein_Starts_M",
                        "Curated_RefSeq_Protein",
                        "Curated_RefSeq_CDS",
                        "Protein_Exact_Duplicate_Count",
                        "CDS_Exact_Duplicate_Count",
                        "Full_Reference_Coverage",
                        "Full_Aligned_Identity",
                        "ADD_Reference_Coverage",
                        "ADD_Aligned_Identity",
                        "CTD_Reference_Coverage",
                        "CTD_Aligned_Identity",
                        "Pairwise_h358_Query_Position",
                        "Pairwise_h358_Residue",
                        "Pairwise_Alignment_Score",
                        "ADD_Architecture_Complete",
                        "Sequence_QC_Pass",
                        "Domain_QC_Pass",
                    ]
                },
                "Family_Tree_Orthology_Status": "PENDING",
                "Preliminary_Alignment_Ready": preliminary_ready,
            }
        )

    evidence = pd.DataFrame(evidence_rows)
    evidence.to_csv(
        RESULTS / "DNMT3L_337_preliminary_orthology_matrix.tsv", sep="\t", index=False
    )
    ready = evidence.loc[evidence["Preliminary_Alignment_Ready"]].copy()
    write_fasta(ready, "protein_seq", PROCESSED / "DNMT3L_337_preliminary_ready_protein.fasta")
    write_fasta(ready, "cds_seq", PROCESSED / "DNMT3L_337_preliminary_ready_cds.fasta")

    subset_definitions = {
        "Mammalia": {"Afrotheria", "Euarchontoglires", "Laurasiatheria", "Xenarthra", "Marsupialia"},
        "Placentalia": {"Afrotheria", "Euarchontoglires", "Laurasiatheria", "Xenarthra"},
        "Euarchontoglires": {"Euarchontoglires"},
        "Laurasiatheria": {"Laurasiatheria"},
        "Sauropsida": {"Sauropsida"},
        "Amphibia": {"Amphibia"},
    }
    subset_rows: list[dict[str, object]] = []
    for subset, clades in subset_definitions.items():
        selected = ready.loc[ready["clade"].isin(clades)].copy()
        write_fasta(selected, "protein_seq", PROCESSED / f"DNMT3L_337_{subset}_protein.fasta")
        write_fasta(selected, "cds_seq", PROCESSED / f"DNMT3L_337_{subset}_cds.fasta")
        subset_rows.append(
            {
                "subset": subset,
                "included_clades": ";".join(sorted(clades)),
                "n": len(selected),
            }
        )
    pd.DataFrame(subset_rows).to_csv(
        RESULTS / "ncbi_337_preliminary_ready_subsets.tsv", sep="\t", index=False
    )

    summary_rows = [
        ("NCBI_ortholog_group_GeneIDs", len(evidence)),
        ("GeneIDs_with_exact_pair", len(selected_by_gene)),
        ("GeneIDs_without_exact_pair", int((evidence["exact_pair_candidate_count"] == 0).sum())),
        ("Sequence_QC_pass", int(evidence["Sequence_QC_Pass"].eq(True).sum())),
        ("Domain_QC_pass", int(evidence["Domain_QC_Pass"].eq(True).sum())),
        ("Hybrid_taxa", int(evidence["Hybrid_Taxon"].sum())),
        ("Preliminary_alignment_ready", int(evidence["Preliminary_Alignment_Ready"].sum())),
        ("Family_tree_validation", "PENDING"),
    ]
    pd.DataFrame(summary_rows, columns=["metric", "value"]).to_csv(
        RESULTS / "ncbi_337_representative_selection_summary.tsv", sep="\t", index=False
    )
    print(pd.DataFrame(summary_rows, columns=["metric", "value"]).to_string(index=False))
    print("\nPreliminary-ready by clade")
    print(ready.groupby("clade").size().rename("n").to_string())
    print("\nNot preliminary-ready")
    print(
        evidence.loc[~evidence["Preliminary_Alignment_Ready"],
                     ["species", "gene_id", "protein_accession", "Sequence_QC_Pass", "Domain_QC_Pass"]]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
