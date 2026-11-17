#!/usr/bin/env python3
"""Audit and descriptively summarize target genes in the nine-species H5AD panel."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csr_matrix


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def text(values):
    return np.asarray([value.decode() if isinstance(value, bytes) else str(value) for value in values])


def obs_column(group, name):
    item = group[name]
    if isinstance(item, h5py.Dataset):
        return text(item[:])
    categories = text(item["categories"][:])
    codes = item["codes"][:]
    return np.asarray([categories[code] if code >= 0 else "NA" for code in codes])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad-manifest", type=Path, required=True)
    parser.add_argument("--feature-manifest", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out-audit", type=Path, required=True)
    parser.add_argument("--out-scale", type=Path, required=True)
    parser.add_argument("--out-summary", type=Path, required=True)
    args = parser.parse_args()

    h5ad_rows = read_tsv(args.h5ad_manifest)
    feature_rows = read_tsv(args.feature_manifest)
    features_by_species = {}
    for row in feature_rows:
        features_by_species.setdefault(row["species"], []).append(row)

    audits, scales, summaries = [], [], []
    for source in h5ad_rows:
        species = source["species"]
        path = args.input_dir / source["filename"]
        with h5py.File(path, "r") as handle:
            genes = text(handle["var"]["_index"][:])
            gene_to_column = {gene: index for index, gene in enumerate(genes)}
            requested = features_by_species[species]
            for feature in requested:
                identifier = feature["h5ad_gene_id"]
                audits.append({
                    **feature,
                    "h5ad_file": source["filename"],
                    "feature_present": "yes" if identifier and identifier in gene_to_column else "no",
                    "h5ad_feature_n": len(genes),
                })

            xgroup = handle["X"]
            sample = xgroup["data"][: min(100000, xgroup["data"].shape[0])]
            finite = sample[np.isfinite(sample)]
            scales.append({
                "species": species,
                "h5ad_file": source["filename"],
                "x_dtype": str(xgroup["data"].dtype),
                "sample_nonzero_n": len(sample),
                "sample_noninteger_fraction": float(np.mean(np.abs(finite - np.rint(finite)) > 1e-6)) if len(finite) else "",
                "sample_min": float(np.min(finite)) if len(finite) else "",
                "sample_max": float(np.max(finite)) if len(finite) else "",
                "analysis_scale": "processed_noninteger" if len(finite) and np.any(np.abs(finite - np.rint(finite)) > 1e-6) else "integer_like",
            })

            admitted = [feature for feature in requested if feature["h5ad_gene_id"] in gene_to_column]
            if not admitted:
                continue
            shape = tuple(int(value) for value in xgroup.attrs["shape"])
            matrix = csr_matrix((xgroup["data"][:], xgroup["indices"][:], xgroup["indptr"][:]), shape=shape)
            obs = handle["obs"]
            library = obs_column(obs, "library_id")
            author_type = obs_column(obs, "cellTypeAuthorAnnotation")
            harmonized_type = obs_column(obs, "cellTypeName")
            for library_id in sorted(set(library)):
                library_selector = library == library_id
                for author_cell_type in sorted(set(author_type[library_selector])):
                    selector = library_selector & (author_type == author_cell_type)
                    indices = np.flatnonzero(selector)
                    harmonized = sorted(set(harmonized_type[selector]))
                    for feature in admitted:
                        values = matrix[indices, gene_to_column[feature["h5ad_gene_id"]]].toarray().ravel()
                        nonzero = int(np.count_nonzero(values))
                        summaries.append({
                            "species": species,
                            "library_id": library_id,
                            "author_cell_type": author_cell_type,
                            "harmonized_cell_type": ";".join(harmonized),
                            "gene_label": feature["gene_label"],
                            "stable_gene_id": feature["h5ad_gene_id"],
                            "n_nuclei": len(indices),
                            "n_nonzero": nonzero,
                            "detection_fraction": f"{nonzero / len(indices):.8f}",
                            "matrix_value_sum": float(values.sum()),
                            "matrix_value_mean_per_nucleus": float(values.mean()),
                            "analysis_scale": scales[-1]["analysis_scale"],
                        })
            del matrix
        print(f"{species}: audited", flush=True)

    write_tsv(args.out_audit, audits, list(audits[0]))
    write_tsv(args.out_scale, scales, list(scales[0]))
    write_tsv(args.out_summary, summaries, list(summaries[0]))


if __name__ == "__main__":
    main()
