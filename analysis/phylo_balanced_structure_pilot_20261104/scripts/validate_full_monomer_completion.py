#!/usr/bin/env python3
"""Validate all 72 balanced-panel monomer predictions and provenance hashes."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def one(path, pattern):
    found = [item for item in path.glob(pattern) if not item.name.startswith("._")]
    if len(found) != 1:
        raise ValueError((path, pattern, len(found)))
    return found[0]


def sha256(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main():
    status = read_tsv(ROOT / "metadata/full_monomer_run_status.tsv")
    if len(status) != 72 or len({(row["target"], row["seed"]) for row in status}) != 72:
        raise ValueError("Expected 72 unique target-seed rows")
    allowed = {"VALIDATED", "VALIDATED_REUSED"}
    if any(row["status"] not in allowed for row in status):
        raise ValueError("Nonvalidated status remains")
    manifest = {row["safe_id"]: int(row["monomer_length"]) for row in read_tsv(ROOT / "metadata/panel_manifest.tsv")}

    for row in status:
        output = ROOT / row["output"]
        cif = one(output, "*.cif")
        confidence = json.loads(one(output, "confidence_*.json").read_text())
        if not confidence["structure_checks_passed"]:
            raise ValueError(f"Structure checks failed: {output}")
        metadata = json.loads(one(output, "run_metadata.json").read_text())
        for name, expected in metadata["output_sha256"].items():
            if sha256(output / name) != expected:
                raise ValueError(f"Hash mismatch: {output / name}")
        structure = MMCIFParser(QUIET=True).get_structure("prediction", cif)[0]
        chain_lengths = [len(list(chain.get_residues())) for chain in structure]
        expected_length = manifest[row["target"]]
        if chain_lengths != [expected_length]:
            raise ValueError((output, chain_lengths, expected_length))
        plddt = np.load(one(output, "plddt_*.npz"))["plddt"]
        pae = np.load(one(output, "pae_*.npz"))["pae"]
        pde = np.load(one(output, "pde_*.npz"))["pde"]
        if plddt.shape != (expected_length,) or pae.shape != (expected_length, expected_length) or pde.shape != (expected_length, expected_length):
            raise ValueError(f"Array dimension mismatch: {output}")
        if not np.isfinite(plddt).all() or not np.isfinite(pae).all() or not np.isfinite(pde).all():
            raise ValueError(f"Nonfinite confidence values: {output}")
    print("VALIDATED: 72/72 balanced-panel monomers; 24 targets x 3 seeds; all hashes and dimensions pass")


if __name__ == "__main__":
    main()
