#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

from Bio import AlignIO, Phylo


PROJECT = Path(__file__).resolve().parents[1]
PANEL = PROJECT / "metadata/lineage_panel.tsv"
ANNOTATION = PROJECT / "results/vertebrate_dnmt3l_annotation_screen.tsv"
NEIGHBORHOOD = PROJECT / "results/conserved_neighborhood_screen.tsv"
ORTHOLOGY = PROJECT / "results/DNMT3L_337_orthology_evidence_matrix.tsv"
COMPACT_SUMMARY = PROJECT / "results/syntenic_interval_dnmt3_summary_expanded.tsv"
COMPACT_FULL_HITS = PROJECT / "results/syntenic_interval_diamond_hits_full_family_expanded.tsv"
MARKER_FULL_HITS = PROJECT / "results/marker_centered_diamond_hits_full_family.tsv"
LUNGFISH_CANDIDATE = PROJECT / "analysis/locus_history_337/lungfish_candidate_extraction_summary.tsv"
FAMILY_ALIGNMENT = PROJECT / "analysis/gene_family_tree_337/DNMT3_family_337_MAFFT.fasta"
TARGETED_TREE = (
    PROJECT
    / "analysis/locus_history_337/targeted_partial_candidate_panel_human_columns_LG_F_R4.treefile"
)
OUTPUT = PROJECT / "results/vertebrate_dnmt3l_locus_matrix.tsv"
SUMMARY = PROJECT / "results/vertebrate_dnmt3l_locus_summary.tsv"
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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_diamond(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        lines = (line for line in handle if line.strip() and not line.startswith("#"))
        return list(csv.DictReader(lines, fieldnames=DIAMOND_FIELDS, delimiter="\t"))


def truth(value: str) -> bool:
    return str(value).strip().lower() == "true"


def best_by_panel(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["qseqid"].split("|", 1)[0]].append(row)
    return {
        panel_id: max(panel_rows, key=lambda row: float(row["bitscore"]))
        for panel_id, panel_rows in grouped.items()
    }


def subject_core_overlap(
    hit: dict[str, str] | None, aligned: dict[str, str], human: str
) -> tuple[int | None, float | None, bool, str]:
    if hit is None:
        return None, None, False, "NO_HIT"
    if hit["sseqid"] not in aligned:
        return None, None, False, "UNAVAILABLE_SUBJECT_NOT_IN_ALIGNMENT"
    subject = aligned[hit["sseqid"]]
    start = min(int(hit["sstart"]), int(hit["send"]))
    end = max(int(hit["sstart"]), int(hit["send"]))
    subject_position = 0
    mapped = 0
    for index, residue in enumerate(subject):
        if residue != "-":
            subject_position += 1
        if start <= subject_position <= end and residue != "-" and human[index] != "-":
            mapped += 1
    alignment_length = max(1, int(hit["length"]))
    fraction = mapped / alignment_length
    compatible = mapped >= 40 and fraction >= 0.30
    return mapped, fraction, compatible, "PASS" if compatible else "FAIL"


def targeted_distances() -> dict[str, dict[str, object]]:
    tree = Phylo.read(TARGETED_TREE, "newick")
    tips = [tip.name for tip in tree.get_terminals()]
    dnmt3l = [tip for tip in tips if tip.startswith("DNMT3L__")]
    paralogs = [tip for tip in tips if not tip.startswith("DNMT3L__")]
    results: dict[str, dict[str, object]] = {}
    for tip in dnmt3l:
        nearest_l_distance, nearest_l = min(
            (tree.distance(tip, other), other) for other in dnmt3l if other != tip
        )
        nearest_p_distance, nearest_p = min(
            (tree.distance(tip, other), other) for other in paralogs
        )
        results[tip.removeprefix("DNMT3L__")] = {
            "nearest_dnmt3l": nearest_l,
            "nearest_dnmt3l_distance": nearest_l_distance,
            "nearest_paralog": nearest_p,
            "nearest_paralog_distance": nearest_p_distance,
            "paralog_to_dnmt3l_distance_ratio": nearest_p_distance / nearest_l_distance,
        }
    return results


def main() -> None:
    panel = read_tsv(PANEL)
    annotation = {row["panel_id"]: row for row in read_tsv(ANNOTATION)}
    neighborhood = {row["panel_id"]: row for row in read_tsv(NEIGHBORHOOD)}
    compact = {row["panel_id"]: row for row in read_tsv(COMPACT_SUMMARY)}
    orthology = {row["species"]: row for row in read_tsv(ORTHOLOGY)}
    compact_hits = best_by_panel(read_diamond(COMPACT_FULL_HITS))
    marker_hits = best_by_panel(read_diamond(MARKER_FULL_HITS))
    lungfish_candidate = {
        row["panel_id"]: row for row in read_tsv(LUNGFISH_CANDIDATE)
    }
    targeted = targeted_distances()

    family_alignment = AlignIO.read(FAMILY_ALIGNMENT, "fasta")
    aligned = {record.id: str(record.seq) for record in family_alignment}
    human_aligned = aligned["DNMT3L__Homo_sapiens"]

    output_rows: list[dict[str, object]] = []
    for panel_row in panel:
        panel_id = panel_row["panel_id"]
        species = panel_row["species"]
        ann = annotation[panel_id]
        hood = neighborhood[panel_id]
        orth = orthology.get(species)
        compact_row = compact.get(panel_id)
        full_hit = compact_hits.get(panel_id) or marker_hits.get(panel_id)
        mapped_count, mapped_fraction, core_compatible, core_mapping_status = subject_core_overlap(
            full_hit, aligned, human_aligned
        )
        target = targeted.get(orth["safe_id"] if orth else "", {})
        targeted_ratio = float(target.get("paralog_to_dnmt3l_distance_ratio", 0.0))
        marker_candidate = lungfish_candidate.get(panel_id, {})

        if orth and truth(orth["Selection_Ready"]) and orth["Family_Tree_Orthology_Status"] == "PASS":
            evidence_grade = "A_CONFIRMED_PRESENT"
            presence_call = "PRESENT"
            locus_status = "CONFIRMED_DNMT3L_LOCUS"
            rationale = "exact protein-CDS pair, sequence/domain QC, and broad family-tree validation"
        elif (
            orth
            and truth(orth["NCBI_DNMT3L_Ortholog_Group"])
            and truth(orth["Translation_Exact_Recalculated"])
            and targeted_ratio > 1.0
        ):
            syntenic_core = bool(
                compact_row
                and compact_row["screen_interpretation"] == "DNMT3L_like_translated_match_detected"
                and core_compatible
            )
            evidence_grade = "B1_SUPPORTED_PARTIAL_SYNTENIC" if syntenic_core else "B2_SUPPORTED_PARTIAL_NO_CORE_HIT"
            presence_call = "PRESENT"
            locus_status = "SUPPORTED_PARTIAL_DNMT3L_LOCUS"
            if syntenic_core:
                suffix = " with a conserved-neighborhood core hit"
            elif (
                compact_row
                and compact_row["screen_interpretation"] == "DNMT3L_like_translated_match_detected"
                and core_mapping_status == "UNAVAILABLE_SUBJECT_NOT_IN_ALIGNMENT"
            ):
                suffix = " with a conserved-neighborhood translated hit whose core overlap was not evaluable"
            elif compact_row and compact_row["screen_interpretation"] == "DNMT3L_like_translated_match_detected":
                suffix = " with a conserved-neighborhood translated hit below the core-overlap safeguard"
            else:
                suffix = (
                    " without a DNMT3-family hit in the primary representative-database "
                    "compact-interval search at E-value 1e-10"
                )
            rationale = "exact translated ortholog-group record and targeted-tree proximity to DNMT3L" + suffix
        elif hood["neighborhood_screen_status"] == "COLOCALIZED_MARKERS":
            evidence_grade = "C_NOT_DETECTED_CONSERVED_NEIGHBORHOOD"
            presence_call = "NOT_DETECTED"
            locus_status = "NO_DNMT3L_CORE_DETECTED_IN_SEARCHED_INTERVAL"
            rationale = (
                "no annotation and no DNMT3L-core match in the primary representative-database "
                "compact-interval search at E-value 1e-10"
            )
        else:
            evidence_grade = "D_UNRESOLVED_MARKER_CENTERED"
            presence_call = "UNRESOLVED"
            locus_status = "MARKER_CENTERED_SEARCH_ONLY"
            rationale = "no compact multi-marker neighborhood; marker-centered evidence cannot establish absence"
            if marker_candidate:
                rationale += (
                    "; translated candidate is "
                    f"{marker_candidate['translated_length_aa']} aa and "
                    f"{marker_candidate['percent_identity']}% identical to the top DNMT3L reference across "
                    f"{marker_candidate['alignment_length_aa']} aligned aa, with one in-frame "
                    f"{marker_candidate['stop_codon']} stop at aa {marker_candidate['stop_position_aa']}; "
                    "pseudogene versus assembly or gene-model error is unresolved"
                )

        representative_core_hit = bool(
            compact_row
            and compact_row["screen_interpretation"] == "DNMT3L_like_translated_match_detected"
            and core_compatible
        )
        if representative_core_hit:
            hit_interpretation = "REPRESENTATIVE_DNMT3L_CORE_HIT"
        elif full_hit and full_hit["sseqid"].startswith("DNMT3L__") and core_mapping_status == "UNAVAILABLE_SUBJECT_NOT_IN_ALIGNMENT":
            hit_interpretation = "DNMT3L_LABELLED_HIT_CORE_MAPPING_UNAVAILABLE"
        elif full_hit and full_hit["sseqid"].startswith("DNMT3L__") and not core_compatible:
            hit_interpretation = "DNMT3L_LABELLED_HIT_BELOW_CORE_THRESHOLD"
        elif full_hit and core_compatible:
            hit_interpretation = "DNMT3L_CORE_COMPATIBLE_HIT"
        elif full_hit:
            hit_interpretation = "NON_DNMT3L_OR_UNMAPPABLE_HIT"
        else:
            hit_interpretation = "NO_FULL_FAMILY_HIT_AT_THRESHOLD"

        output_rows.append(
            {
                "panel_id": panel_id,
                "major_lineage": panel_row["major_lineage"],
                "species": species,
                "presence_call": presence_call,
                "evidence_grade": evidence_grade,
                "locus_status": locus_status,
                "loss_inference_allowed": False,
                "annotation_screen_status": ann["annotation_screen_status"],
                "annotation_gene_ids": ann["gene_ids"],
                "annotation_protein_accessions": ann["protein_accessions"],
                "orthology_matrix_record": bool(orth),
                "protein_accession": orth["protein_accession"] if orth else "",
                "cds_accession": orth["cds_accession"] if orth else "",
                "assembly_accession": orth["NCBI_Assembly_Accession"] if orth else "",
                "genomic_accession": orth["NCBI_Genomic_Accession"] if orth else hood["best_genomic_accession"],
                "translation_exact": orth["Translation_Exact_Recalculated"] if orth else "",
                "sequence_qc_pass": orth["Sequence_QC_Pass"] if orth else "",
                "domain_qc_pass": orth["Domain_QC_Pass"] if orth else "",
                "broad_family_tree_status": orth["Family_Tree_Orthology_Status"] if orth else "",
                "targeted_nearest_dnmt3l": target.get("nearest_dnmt3l", ""),
                "targeted_nearest_paralog": target.get("nearest_paralog", ""),
                "targeted_paralog_to_dnmt3l_distance_ratio": (
                    f"{targeted_ratio:.6f}" if targeted_ratio else ""
                ),
                "neighborhood_screen_status": hood["neighborhood_screen_status"],
                "colocalized_markers": hood["colocalized_markers"],
                "translated_search_interpretation": (
                    compact_row["screen_interpretation"] if compact_row else "marker_centered_only"
                ),
                "full_family_top_subject": full_hit["sseqid"] if full_hit else "",
                "full_family_top_evalue": full_hit["evalue"] if full_hit else "",
                "full_family_top_bitscore": full_hit["bitscore"] if full_hit else "",
                "top_hit_human_coordinate_mappable_residues": "" if mapped_count is None else mapped_count,
                "top_hit_human_coordinate_mappable_fraction": "" if mapped_fraction is None else f"{mapped_fraction:.6f}",
                "top_hit_core_mapping_status": core_mapping_status,
                "top_hit_core_rule_pass": core_compatible,
                "top_hit_core_interpretation": hit_interpretation,
                "marker_candidate_translated_length_aa": marker_candidate.get("translated_length_aa", ""),
                "marker_candidate_percent_identity": marker_candidate.get("percent_identity", ""),
                "marker_candidate_alignment_length_aa": marker_candidate.get("alignment_length_aa", ""),
                "marker_candidate_stop_count": marker_candidate.get("stop_count_in_direct_translation", ""),
                "marker_candidate_stop_position_aa": marker_candidate.get("stop_position_aa", ""),
                "marker_candidate_stop_codon": marker_candidate.get("stop_codon", ""),
                "marker_candidate_interpretation": marker_candidate.get("interpretation", ""),
                "evidence_rationale": rationale,
                "absence_caveat": (
                    "representative assembled interval only; not evidence of gene loss"
                    if presence_call != "PRESENT"
                    else ""
                ),
            }
        )

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    by_lineage: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in output_rows:
        by_lineage[str(row["major_lineage"])].append(row)
    summary_rows: list[dict[str, object]] = []
    for lineage, rows in sorted(by_lineage.items()):
        counts = Counter(str(row["presence_call"]) for row in rows)
        if counts["PRESENT"]:
            interpretation = "DNMT3L_PRESENT_IN_AT_LEAST_ONE_REPRESENTATIVE"
        elif counts["NOT_DETECTED"] >= 2 and counts["UNRESOLVED"] == 0:
            interpretation = "CANDIDATE_MULTISPECIES_LINEAGE_ABSENCE_NOT_PROVEN_LOSS"
        elif counts["NOT_DETECTED"] and counts["UNRESOLVED"] == 0:
            interpretation = "NOT_DETECTED_IN_ONE_CONSERVED_NEIGHBORHOOD"
        else:
            interpretation = "UNRESOLVED"
        summary_rows.append(
            {
                "major_lineage": lineage,
                "n_species": len(rows),
                "n_present": counts["PRESENT"],
                "n_not_detected": counts["NOT_DETECTED"],
                "n_unresolved": counts["UNRESOLVED"],
                "lineage_interpretation": interpretation,
                "species": ";".join(str(row["species"]) for row in rows),
                "evidence_grades": ";".join(str(row["evidence_grade"]) for row in rows),
            }
        )
    with SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote {OUTPUT.relative_to(PROJECT)} ({len(output_rows)} species)")
    print(f"Wrote {SUMMARY.relative_to(PROJECT)} ({len(summary_rows)} lineages)")
    print(Counter(str(row["evidence_grade"]) for row in output_rows))


if __name__ == "__main__":
    main()
