#!/usr/bin/env python3
"""Build conservative non-mammalian DNMT3L expression-support tables."""

from __future__ import annotations

import csv
import gzip
import hashlib
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LOCUS = ROOT / "results" / "vertebrate_dnmt3l_locus_matrix.tsv"
ANOLIS = ROOT / "data" / "raw" / "expression_processed" / "GSE97367_anolis.7sp.fpkm.txt.gz"
AXO = ROOT / "data" / "raw" / "expression_processed" / "Axo.rsem.isoforms.TMM.fpkm.wAnnot.matrix.gz"
OUT = ROOT / "results" / "nonmammalian_expression_support.tsv"
DETAIL = ROOT / "results" / "nonmammalian_expression_support_by_context.tsv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_loci() -> dict[str, dict[str, str]]:
    with LOCUS.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {
        row["panel_id"]: row
        for row in rows
        if row["major_lineage"] not in {"Marsupialia", "Placentalia"}
        and row["evidence_grade"].startswith(("A_", "B1_", "B2_"))
    }


def read_anolis() -> tuple[list[dict[str, object]], dict[str, object]]:
    stable_id = "ENSACAG00000028767"
    with gzip.open(ANOLIS, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", quotechar='"')
        header = next(reader)
        row = next(r for r in reader if r[0] == stable_id)
    contexts = header[:10]
    values = [float(x) for x in row[1:11]]
    detail = [
        {
            "panel_id": "ANOLE",
            "dataset_id": "ANOLIS_MULTI_TISSUE_2017",
            "context": context,
            "context_unit": "sex_by_tissue_group_mean",
            "processed_measure": "FPKM",
            "processed_value": f"{value:.12g}",
            "detected": str(value > 0),
        }
        for context, value in zip(contexts, values)
    ]
    summary = {
        "dataset_id": "ANOLIS_MULTI_TISSUE_2017",
        "repository_accession": "GSE97367",
        "expression_state": "GROUP_SUMMARY_ONLY",
        "mapping_method": "exact stable Ensembl gene identifier crosswalked from NCBI GeneID 103278356",
        "mapped_identifier": stable_id,
        "processed_measure": "FPKM",
        "context_unit": "sex_by_tissue_group_mean",
        "contexts_tested": len(values),
        "contexts_positive": sum(v > 0 for v in values),
        "positive_contexts": ";".join(c for c, v in zip(contexts, values) if v > 0),
        "maximum_processed_value": f"{max(values):.12g}",
        "artifact_sha256": sha256(ANOLIS),
        "evidence_scope": "secondary descriptive transcription evidence; biological-library counts unavailable",
        "caveat": "group-level table fails the frozen sample-level primary gate; values are not compared across species",
    }
    return detail, summary


def read_axolotl() -> tuple[list[dict[str, object]], dict[str, object]]:
    # Two nearly full-length predicted proteins match the current axolotl
    # DNMT3L candidate at 100% identity across residues 6-341. The two shorter
    # isoforms belong to the same Trinity gene and retain an exact C-terminal
    # segment, so locus-level abundance is the sum of all four isoforms.
    isoforms = {f"c1086735_g1_i{i}" for i in range(1, 5)}
    values_by_isoform: dict[str, list[float]] = {}
    with gzip.open(AXO, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        for row in reader:
            transcript_id = row[0].split("^", 1)[0]
            if transcript_id in isoforms:
                values_by_isoform[transcript_id] = [float(x) for x in row[1:]]
    missing = isoforms - values_by_isoform.keys()
    if missing:
        raise RuntimeError(f"Missing expected axolotl isoforms: {sorted(missing)}")
    contexts = header[1:] if header and header[0] == "" else header
    sums = [sum(values_by_isoform[t][i] for t in sorted(isoforms)) for i in range(len(contexts))]
    detail = [
        {
            "panel_id": "AXOLOTL",
            "dataset_id": "AXOLOTL_MULTI_TISSUE_2016",
            "context": context,
            "context_unit": "biological_library",
            "processed_measure": "summed_isoform_TMM_normalized_FPKM",
            "processed_value": f"{value:.12g}",
            "detected": str(value > 0),
        }
        for context, value in zip(contexts, sums)
    ]
    summary = {
        "dataset_id": "AXOLOTL_MULTI_TISSUE_2016",
        "repository_accession": "GSE92429; GFBM000000000",
        "expression_state": "SUPPORTED" if any(v > 0 for v in sums) else "NOT_DETECTED",
        "mapping_method": "100% protein identity over 336 aa to two Trinity predicted proteins; same-gene isoforms summed",
        "mapped_identifier": "c1086735_g1_i1;c1086735_g1_i2;c1086735_g1_i3;c1086735_g1_i4",
        "processed_measure": "summed_isoform_TMM_normalized_FPKM",
        "context_unit": "biological_library",
        "contexts_tested": len(sums),
        "contexts_positive": sum(v > 0 for v in sums),
        "positive_contexts": ";".join(c for c, v in zip(contexts, sums) if v > 0),
        "maximum_processed_value": f"{max(sums):.12g}",
        "artifact_sha256": sha256(AXO),
        "evidence_scope": "sample-level transcription support in the named tissues",
        "caveat": "legacy de novo isoforms; transcription does not establish biochemical function",
    }
    return detail, summary


def main() -> None:
    loci = read_loci()
    anole_detail, anole_summary = read_anolis()
    axo_detail, axo_summary = read_axolotl()
    measured = {"ANOLE": anole_summary, "AXOLOTL": axo_summary}
    rows = []
    for panel_id, locus in loci.items():
        base = {
            "panel_id": panel_id,
            "major_lineage": locus["major_lineage"],
            "species": locus["species"],
            "locus_evidence_grade": locus["evidence_grade"],
        }
        if panel_id in measured:
            base.update(measured[panel_id])
        else:
            base.update(
                {
                    "dataset_id": "PELODISCUS_GONAD_2019;PELODISCUS_TESTIS_2016" if panel_id == "TURTLE" else "NONE_IDENTIFIED",
                    "repository_accession": "PRJNA528401;PRJNA353665;PRJNA353709" if panel_id == "TURTLE" else "",
                    "expression_state": "NO_PROCESSED_ARTIFACT",
                    "mapping_method": "not tested",
                    "mapped_identifier": "",
                    "processed_measure": "",
                    "context_unit": "",
                    "contexts_tested": 0,
                    "contexts_positive": 0,
                    "positive_contexts": "",
                    "maximum_processed_value": "",
                    "artifact_sha256": "",
                    "evidence_scope": "none",
                    "caveat": "raw-only public data were not downloaded for this presence-support endpoint" if panel_id == "TURTLE" else "no suitable processed public matrix was identified in the frozen screen",
                }
            )
        rows.append(base)

    fields = [
        "panel_id", "major_lineage", "species", "locus_evidence_grade", "dataset_id",
        "repository_accession", "expression_state", "mapping_method", "mapped_identifier",
        "processed_measure", "context_unit", "contexts_tested", "contexts_positive",
        "positive_contexts", "maximum_processed_value", "artifact_sha256", "evidence_scope", "caveat",
    ]
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    detail_fields = ["panel_id", "dataset_id", "context", "context_unit", "processed_measure", "processed_value", "detected"]
    with DETAIL.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=detail_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(axo_detail + anole_detail)

    print(f"Wrote {OUT.relative_to(ROOT)} ({len(rows)} eligible non-mammalian loci)")
    print(f"Wrote {DETAIL.relative_to(ROOT)} ({len(axo_detail) + len(anole_detail)} contexts)")


if __name__ == "__main__":
    main()
