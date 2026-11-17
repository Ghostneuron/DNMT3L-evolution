#!/usr/bin/env python3
"""Inventory an H5AD file without requiring scanpy/anndata."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py


def decode(values):
    out = []
    for value in values:
        if isinstance(value, bytes):
            out.append(value.decode("utf-8", errors="replace"))
        else:
            out.append(str(value))
    return out


def summarize_group(group: h5py.Group) -> dict:
    summary = {}
    for name, item in group.items():
        if isinstance(item, h5py.Dataset):
            summary[name] = {
                "kind": "dataset",
                "shape": list(item.shape),
                "dtype": str(item.dtype),
            }
        else:
            summary[name] = {"kind": "group", "keys": sorted(item.keys())}
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    result = {"path": str(args.h5ad.resolve()), "size_bytes": args.h5ad.stat().st_size}
    with h5py.File(args.h5ad, "r") as handle:
        result["root_keys"] = sorted(handle.keys())
        result["root_attributes"] = {key: str(value) for key, value in handle.attrs.items()}
        for group_name in ("obs", "var", "layers", "obsm", "uns"):
            if group_name in handle:
                result[group_name] = summarize_group(handle[group_name])
        for axis in ("obs", "var"):
            if axis in handle and "_index" in handle[axis]:
                ds = handle[axis]["_index"]
                result[f"{axis}_n"] = int(ds.shape[0])
                result[f"{axis}_index_preview"] = decode(ds[: min(10, ds.shape[0])])
        if "X" in handle:
            item = handle["X"]
            if isinstance(item, h5py.Dataset):
                result["X"] = {"encoding": "dense", "shape": list(item.shape), "dtype": str(item.dtype)}
            else:
                result["X"] = {"encoding": "group", "keys": sorted(item.keys())}
                if "shape" in item:
                    result["X"]["shape"] = [int(x) for x in item["shape"][:]]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".partial")
    temporary.write_text(json.dumps(result, indent=2) + "\n")
    temporary.replace(args.output)


if __name__ == "__main__":
    main()

