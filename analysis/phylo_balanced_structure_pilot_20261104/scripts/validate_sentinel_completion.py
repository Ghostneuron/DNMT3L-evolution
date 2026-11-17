#!/usr/bin/env python3
"""Validate that every archived sentinel prediction is complete and parseable."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser


ROOT = Path(__file__).resolve().parents[1]


def tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def files(path, pattern):
    return [item for item in path.glob(pattern) if not item.name.startswith("._")]


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main():
    status = tsv(ROOT / "metadata/sentinel_run_status.tsv")
    assert len(status) == 42, f"Expected 42 status rows, found {len(status)}"
    assert len({(row['mode'], row['target'], row['seed']) for row in status}) == 42
    assert all(row["status"] == "VALIDATED" for row in status)
    manifest = {row["safe_id"]: row for row in tsv(ROOT / "metadata/panel_manifest.tsv")}

    for row in status:
        output = ROOT / row["output"]
        found = {
            "cif": files(output, "*.cif"),
            "confidence": files(output, "confidence_*.json"),
            "pae": files(output, "pae_*.npz"),
            "pde": files(output, "pde_*.npz"),
            "plddt": files(output, "plddt_*.npz"),
            "metadata": files(output, "run_metadata.json"),
        }
        assert all(len(items) == 1 for items in found.values()), (output, {k: len(v) for k, v in found.items()})
        confidence = json.loads(found["confidence"][0].read_text())
        assert confidence["structure_checks_passed"] is True
        metadata = json.loads(found["metadata"][0].read_text())
        for name, expected in metadata["output_sha256"].items():
            assert digest(output / name) == expected, (output, name)

        structure = MMCIFParser(QUIET=True).get_structure("prediction", found["cif"][0])[0]
        chain_lengths = [len(list(chain.get_residues())) for chain in structure]
        panel = manifest[row["target"]]
        expected_chains = [int(panel["monomer_length"])] if row["mode"] == "monomer" else [290, int(panel["ctd_length"])]
        assert chain_lengths == expected_chains, (output, chain_lengths, expected_chains)
        total = sum(expected_chains)
        plddt = np.load(found["plddt"][0])["plddt"]
        pae = np.load(found["pae"][0])["pae"]
        pde = np.load(found["pde"][0])["pde"]
        assert plddt.shape == (total,), (output, plddt.shape)
        assert pae.shape == (total, total), (output, pae.shape)
        assert pde.shape == (total, total), (output, pde.shape)
        assert np.isfinite(plddt).all() and np.isfinite(pae).all() and np.isfinite(pde).all()

    print("VALIDATED: 42/42 sentinel predictions; 21 monomers and 21 complexes; all hashes and dimensions pass")


if __name__ == "__main__":
    main()
