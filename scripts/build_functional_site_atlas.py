#!/usr/bin/env python3
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices


PROJECT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((PROJECT / "config" / "design.json").read_text())
MATRIX = (PROJECT / CONFIG["primary_input_matrix"]).resolve()
REFERENCE_SEQUENCE_ARCHIVE = (PROJECT / CONFIG["reference_sequence_archive"]).resolve()
SITE_MANIFEST = PROJECT / "metadata" / "functional_site_manifest.tsv"
RESULTS = PROJECT / "results"
HUMAN_ACCESSION = CONFIG["human_reference_accession"]
MOUSE_ACCESSION = CONFIG["mouse_reference_accession"]
MIN_LOCAL_IDENTITY = float(CONFIG["minimum_local_identity_11aa"])
MIN_REGION_COVERAGE = float(CONFIG["minimum_site_occupancy"])
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    aligner.end_insertion_score = -2.0
    aligner.end_deletion_score = -2.0
    return aligner


def bidirectional_map(
    reference: str,
    query: str,
    aligner: PairwiseAligner,
) -> tuple[dict[int, tuple[int, str]], dict[int, tuple[int, str]], float]:
    alignment = aligner.align(reference, query)[0]
    reference_to_query: dict[int, tuple[int, str]] = {}
    query_to_reference: dict[int, tuple[int, str]] = {}
    reference_blocks, query_blocks = alignment.aligned
    for (ref_start, ref_end), (query_start, query_end) in zip(reference_blocks, query_blocks):
        block_length = int(ref_end - ref_start)
        if block_length != int(query_end - query_start):
            raise ValueError("Aligned block lengths differ")
        for offset in range(block_length):
            ref_position = int(ref_start + offset + 1)
            query_position = int(query_start + offset + 1)
            reference_to_query[ref_position] = (query_position, query[query_position - 1])
            query_to_reference[query_position] = (ref_position, reference[ref_position - 1])
    return reference_to_query, query_to_reference, float(alignment.score)


def local_context(
    reference: str,
    query: str,
    mapping: dict[int, tuple[int, str]],
    position: int,
    radius: int = 5,
) -> tuple[str, str, int, int, float | None]:
    start = max(1, position - radius)
    end = min(len(reference), position + radius)
    ref_window = reference[start - 1 : end]
    query_chars: list[str] = []
    compared = 0
    identical = 0
    for ref_position in range(start, end + 1):
        mapped = mapping.get(ref_position)
        if mapped is None:
            query_chars.append("-")
            continue
        query_aa = mapped[1]
        query_chars.append(query_aa)
        if query_aa in STANDARD_AA:
            compared += 1
            identical += int(query_aa == reference[ref_position - 1])
    identity = identical / compared if compared else None
    return ref_window, "".join(query_chars), compared, identical, identity


def broad_lineage(clade: str) -> str:
    if clade in {"Afrotheria", "Euarchontoglires", "Laurasiatheria", "Xenarthra"}:
        return "Placentalia"
    if clade == "Marsupialia":
        return "Marsupialia"
    if clade == "Sauropsida":
        return "Sauropsida"
    if clade == "Amphibia":
        return "Amphibia"
    return clade


def count_string(values: pd.Series) -> str:
    counts = Counter(values.dropna().astype(str))
    return ";".join(f"{residue}:{counts[residue]}" for residue in sorted(counts))


def validate_reference_rows(records: pd.DataFrame) -> tuple[str, str]:
    sequences = {
        str(row.protein_accession): str(row.protein_seq).upper()
        for row in records.itertuples(index=False)
        if str(row.protein_accession)
    }
    for record in SeqIO.parse(REFERENCE_SEQUENCE_ARCHIVE, "fasta"):
        sequences.setdefault(record.id, str(record.seq).upper())
    missing = [accession for accession in [HUMAN_ACCESSION, MOUSE_ACCESSION] if accession not in sequences]
    if missing:
        raise ValueError(f"Reference accessions absent from matrix and archive: {missing}")
    return sequences[HUMAN_ACCESSION], sequences[MOUSE_ACCESSION]


def build_reference_crosswalk(
    sites: pd.DataFrame,
    human_sequence: str,
    mouse_sequence: str,
    aligner: PairwiseAligner,
) -> pd.DataFrame:
    _, mouse_to_human, score = bidirectional_map(human_sequence, mouse_sequence, aligner)
    rows: list[dict[str, object]] = []
    for site in sites.loc[sites["feature_type"].eq("site")].to_dict(orient="records"):
        source_position = int(site["source_position"])
        declared_human_position = int(site["human_np_position"])
        if site["source_accession"] == HUMAN_ACCESSION:
            mapped_human_position = source_position
            mapped_human_residue = human_sequence[source_position - 1]
            source_sequence = human_sequence
            mapping_method = "identity"
        elif site["source_accession"] == MOUSE_ACCESSION:
            mapped = mouse_to_human.get(source_position)
            if mapped is None:
                mapped_human_position = None
                mapped_human_residue = None
            else:
                mapped_human_position, mapped_human_residue = mapped
            source_sequence = mouse_sequence
            mapping_method = "global_BLOSUM62"
        else:
            raise ValueError(f"Unsupported source accession: {site['source_accession']}")

        observed_source_residue = source_sequence[source_position - 1]
        status = "PASS"
        notes: list[str] = []
        if observed_source_residue != site["source_residue"]:
            status = "FAIL"
            notes.append("source residue mismatch")
        if mapped_human_position != declared_human_position:
            status = "FAIL"
            notes.append("human position mismatch")
        if mapped_human_residue != site["human_np_residue"]:
            status = "FAIL"
            notes.append("human residue mismatch")
        rows.append(
            {
                "feature_id": site["feature_id"],
                "source_accession": site["source_accession"],
                "source_position": source_position,
                "declared_source_residue": site["source_residue"],
                "observed_source_residue": observed_source_residue,
                "declared_human_position": declared_human_position,
                "mapped_human_position": mapped_human_position,
                "declared_human_residue": site["human_np_residue"],
                "mapped_human_residue": mapped_human_residue,
                "mapping_method": mapping_method,
                "human_mouse_alignment_score": score if mapping_method != "identity" else None,
                "status": status,
                "notes": "; ".join(notes),
            }
        )
    crosswalk = pd.DataFrame(rows)
    failed = crosswalk.loc[~crosswalk["status"].eq("PASS")]
    if not failed.empty:
        raise ValueError("Reference coordinate validation failed:\n" + failed.to_string(index=False))
    return crosswalk


def build_site_atlas(
    ready: pd.DataFrame,
    sites: pd.DataFrame,
    human_sequence: str,
    aligner: PairwiseAligner,
) -> pd.DataFrame:
    site_rows = sites.loc[sites["feature_type"].eq("site")].copy()
    output: list[dict[str, object]] = []
    for record in ready.to_dict(orient="records"):
        query_sequence = str(record["protein_seq"]).upper()
        mapping, _, score = bidirectional_map(human_sequence, query_sequence, aligner)
        for site in site_rows.to_dict(orient="records"):
            human_position = int(site["human_np_position"])
            mapped = mapping.get(human_position)
            query_position = mapped[0] if mapped else None
            query_residue = mapped[1] if mapped else None
            ref_window, query_window, compared, identical, identity = local_context(
                human_sequence, query_sequence, mapping, human_position
            )
            if mapped is None or query_residue not in STANDARD_AA:
                mapping_status = "UNMAPPED"
            elif compared < 7:
                mapping_status = "LOW_CONTEXT_COVERAGE"
            elif identity is None or identity < MIN_LOCAL_IDENTITY:
                mapping_status = "LOW_CONTEXT_IDENTITY"
            else:
                mapping_status = "PASS"
            output.append(
                {
                    "feature_id": site["feature_id"],
                    "atlas_role": site["atlas_role"],
                    "evidence_class": site["evidence_class"],
                    "safe_id": record["safe_id"],
                    "species": record["species"],
                    "gene_id": record["gene_id"],
                    "protein_accession": record["protein_accession"],
                    "clade": record["Clade"],
                    "broad_lineage": broad_lineage(str(record["Clade"])),
                    "human_np_position": human_position,
                    "human_np_residue": human_sequence[human_position - 1],
                    "query_position": query_position,
                    "query_residue": query_residue,
                    "pairwise_alignment_score": score,
                    "reference_window_11aa": ref_window,
                    "query_window_11aa": query_window,
                    "local_compared_positions": compared,
                    "local_identical_positions": identical,
                    "local_identity_11aa": identity,
                    "mapping_status": mapping_status,
                }
            )
    return pd.DataFrame(output)


def summarize_site_atlas(atlas: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_columns = ["feature_id", "atlas_role", "human_np_position", "human_np_residue", "broad_lineage"]
    for key, group in atlas.groupby(group_columns, sort=True, dropna=False):
        feature_id, role, position, human_residue, lineage = key
        resolved = group.loc[group["query_residue"].isin(STANDARD_AA), "query_residue"]
        high_confidence = group.loc[
            group["mapping_status"].eq("PASS") & group["query_residue"].isin(STANDARD_AA),
            "query_residue",
        ]
        counts = Counter(resolved.astype(str))
        pass_counts = Counter(high_confidence.astype(str))
        major, major_n = counts.most_common(1)[0] if counts else (None, 0)
        pass_major, pass_major_n = pass_counts.most_common(1)[0] if pass_counts else (None, 0)
        rows.append(
            {
                "feature_id": feature_id,
                "atlas_role": role,
                "human_np_position": int(position),
                "human_np_residue": human_residue,
                "broad_lineage": lineage,
                "records": len(group),
                "resolved_records": len(resolved),
                "resolved_occupancy": len(resolved) / len(group),
                "high_confidence_records": len(high_confidence),
                "high_confidence_occupancy": len(high_confidence) / len(group),
                "major_residue_all_mapped": major,
                "major_residue_frequency_all_mapped": major_n / len(resolved) if len(resolved) else None,
                "human_state_count_all_mapped": counts.get(human_residue, 0),
                "human_state_frequency_all_mapped": counts.get(human_residue, 0) / len(resolved) if len(resolved) else None,
                "residue_counts_all_mapped": count_string(resolved),
                "major_residue_high_confidence": pass_major,
                "major_residue_frequency_high_confidence": pass_major_n / len(high_confidence) if len(high_confidence) else None,
                "residue_counts_high_confidence": count_string(high_confidence),
            }
        )
    return pd.DataFrame(rows)


def build_region_metrics(
    ready: pd.DataFrame,
    sites: pd.DataFrame,
    human_sequence: str,
    aligner: PairwiseAligner,
) -> pd.DataFrame:
    region_rows = sites.loc[sites["feature_type"].eq("region")]
    output: list[dict[str, object]] = []
    for region in region_rows.to_dict(orient="records"):
        start_text, end_text = str(region["human_np_position"]).split("-", 1)
        start, end = int(start_text), int(end_text)
        reference_region = human_sequence[start - 1 : end]
        for record in ready.to_dict(orient="records"):
            query_sequence = str(record["protein_seq"]).upper()
            mapping, _, score = bidirectional_map(human_sequence, query_sequence, aligner)
            mapped_residues: list[str] = []
            identity_count = 0
            resolved_count = 0
            positive = 0
            negative = 0
            polar_or_charged = 0
            for position in range(start, end + 1):
                mapped = mapping.get(position)
                aa = mapped[1] if mapped else "-"
                mapped_residues.append(aa)
                if aa in STANDARD_AA:
                    resolved_count += 1
                    identity_count += int(aa == human_sequence[position - 1])
                    positive += int(aa in {"K", "R"})
                    negative += int(aa in {"D", "E"})
                    polar_or_charged += int(aa in {"D", "E", "H", "K", "N", "Q", "R", "S", "T"})
            region_length = end - start + 1
            coverage = resolved_count / region_length
            output.append(
                {
                    "feature_id": region["feature_id"],
                    "atlas_role": region["atlas_role"],
                    "safe_id": record["safe_id"],
                    "species": record["species"],
                    "protein_accession": record["protein_accession"],
                    "clade": record["Clade"],
                    "broad_lineage": broad_lineage(str(record["Clade"])),
                    "human_np_start": start,
                    "human_np_end": end,
                    "reference_region": reference_region,
                    "mapped_region": "".join(mapped_residues),
                    "resolved_positions": resolved_count,
                    "region_length": region_length,
                    "coverage": coverage,
                    "identity_fraction_resolved": identity_count / resolved_count if resolved_count else None,
                    "positive_residue_count": positive,
                    "negative_residue_count": negative,
                    "net_basic_minus_acidic": positive - negative,
                    "polar_or_charged_fraction": polar_or_charged / resolved_count if resolved_count else None,
                    "pairwise_alignment_score": score,
                    "coverage_status": "PASS" if coverage >= MIN_REGION_COVERAGE else "LOW_COVERAGE",
                }
            )
    return pd.DataFrame(output)


def summarize_regions(regions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (feature_id, lineage), group in regions.groupby(["feature_id", "broad_lineage"], sort=True):
        passing = group.loc[group["coverage_status"].eq("PASS")]
        rows.append(
            {
                "feature_id": feature_id,
                "broad_lineage": lineage,
                "records": len(group),
                "coverage_pass_records": len(passing),
                "coverage_pass_fraction": len(passing) / len(group),
                "median_identity_fraction_resolved": passing["identity_fraction_resolved"].median(),
                "median_net_basic_minus_acidic": passing["net_basic_minus_acidic"].median(),
                "median_polar_or_charged_fraction": passing["polar_or_charged_fraction"].median(),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    separator = "\t" if MATRIX.suffix == ".tsv" else ","
    records = pd.read_csv(MATRIX, sep=separator, dtype={"gene_id": str})
    if "Clade" not in records and "clade" in records:
        records = records.rename(columns={"clade": "Clade"})
    ready = records.loc[records["Selection_Ready"].astype(bool)].copy()
    ready = ready.loc[~ready["Hybrid_Taxon"].astype(bool)].copy()
    expected = int(CONFIG["expected_selection_ready_records"])
    if len(ready) != expected:
        raise ValueError(f"Expected {expected} nonhybrid selection-ready records, found {len(ready)}")

    sites = pd.read_csv(SITE_MANIFEST, sep="\t", dtype=str, keep_default_na=False)
    human_sequence, mouse_sequence = validate_reference_rows(records)
    aligner = make_aligner()

    crosswalk = build_reference_crosswalk(sites, human_sequence, mouse_sequence, aligner)
    crosswalk.to_csv(RESULTS / "functional_reference_coordinate_crosswalk.tsv", sep="\t", index=False)

    atlas = build_site_atlas(ready, sites, human_sequence, aligner)
    atlas.to_csv(RESULTS / "functional_site_states_by_record.tsv", sep="\t", index=False)
    summary = summarize_site_atlas(atlas)
    summary.to_csv(RESULTS / "functional_site_summary_by_lineage.tsv", sep="\t", index=False)

    regions = build_region_metrics(ready, sites, human_sequence, aligner)
    regions.to_csv(RESULTS / "functional_region_metrics_by_record.tsv", sep="\t", index=False)
    region_summary = summarize_regions(regions)
    region_summary.to_csv(RESULTS / "functional_region_summary_by_lineage.tsv", sep="\t", index=False)

    print(f"Validated {len(crosswalk)} site coordinate definitions")
    print(f"Mapped {atlas['feature_id'].nunique()} site features across {ready['safe_id'].nunique()} records")
    print(f"Summarized {regions['feature_id'].nunique()} region feature across {len(regions)} record-region rows")
    print("\nPrimary-site lineage summary")
    print(
        summary.loc[summary["atlas_role"].eq("primary")]
        .loc[:, [
            "feature_id",
            "broad_lineage",
            "records",
            "high_confidence_records",
            "major_residue_high_confidence",
            "residue_counts_high_confidence",
        ]]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()

