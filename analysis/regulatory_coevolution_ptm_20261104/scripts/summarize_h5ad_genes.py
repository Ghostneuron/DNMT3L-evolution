#!/usr/bin/env python3
"""Create donor/library-aware gene summaries from a sparse H5AD matrix."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csr_matrix


def text(values):
    return np.asarray([
        value.decode("utf-8", errors="replace") if isinstance(value, bytes) else str(value)
        for value in values
    ])


def obs_column(group: h5py.Group, name: str) -> np.ndarray:
    item = group[name]
    if isinstance(item, h5py.Dataset):
        return text(item[:])
    categories = text(item["categories"][:])
    codes = item["codes"][:]
    return np.asarray([categories[code] if code >= 0 else "NA" for code in codes])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad", type=Path)
    parser.add_argument("--gene", action="append", required=True, help="LABEL=stable_gene_id")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    requested = []
    for token in args.gene:
        label, stable_id = token.split("=", 1)
        requested.append((label, stable_id))

    with h5py.File(args.h5ad, "r") as handle:
        genes = text(handle["var"]["_index"][:])
        gene_to_column = {gene: index for index, gene in enumerate(genes)}
        missing = [stable_id for _, stable_id in requested if stable_id not in gene_to_column]
        if missing:
            raise SystemExit(f"Missing stable IDs: {', '.join(missing)}")

        xgroup = handle["X"]
        shape = tuple(int(value) for value in xgroup.attrs["shape"])
        matrix = csr_matrix(
            (xgroup["data"][:], xgroup["indices"][:], xgroup["indptr"][:]),
            shape=shape,
        )
        obs = handle["obs"]
        library = obs_column(obs, "library_id")
        author_type = obs_column(obs, "cellTypeAuthorAnnotation")
        harmonized_type = obs_column(obs, "cellTypeName")

        rows = []
        for library_id in sorted(set(library)):
            for author_cell_type in sorted(set(author_type[library == library_id])):
                selector = (library == library_id) & (author_type == author_cell_type)
                indices = np.flatnonzero(selector)
                harmonized = sorted(set(harmonized_type[selector]))
                for label, stable_id in requested:
                    values = matrix[indices, gene_to_column[stable_id]].toarray().ravel()
                    rows.append({
                        "library_id": library_id,
                        "author_cell_type": author_cell_type,
                        "harmonized_cell_type": ";".join(harmonized),
                        "gene_label": label,
                        "stable_gene_id": stable_id,
                        "n_nuclei": len(indices),
                        "n_nonzero": int(np.count_nonzero(values)),
                        "detection_fraction": float(np.count_nonzero(values) / len(indices)),
                        "matrix_value_sum": float(values.sum()),
                        "matrix_value_mean_per_nucleus": float(values.mean()),
                    })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".partial")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
