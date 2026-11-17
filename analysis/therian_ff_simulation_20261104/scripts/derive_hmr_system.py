#!/usr/bin/env python3
"""Derive a 3-amu protein-H HMR system without changing potential energy."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

from openmm import XmlSerializer, unit
from openmm.app import PDBFile


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--hydrogen-mass-amu", type=float, default=3.0)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=False)
    pdb = PDBFile(str(args.source / "minimized_solvated.pdb"))
    system = XmlSerializer.deserialize((args.source / "system.xml").read_text())
    target = args.hydrogen_mass_amu * unit.dalton

    bonded_neighbor: dict[int, int] = {}
    for atom1, atom2 in pdb.topology.bonds():
        if atom1.element.symbol == "H" and atom2.element.symbol != "H":
            bonded_neighbor[atom1.index] = atom2.index
        elif atom2.element.symbol == "H" and atom1.element.symbol != "H":
            bonded_neighbor[atom2.index] = atom1.index

    changed = 0
    transferred = 0.0 * unit.dalton
    for atom in pdb.topology.atoms():
        if atom.element.symbol != "H" or atom.residue.name in {"HOH", "WAT"}:
            continue
        heavy = bonded_neighbor.get(atom.index)
        if heavy is None:
            raise RuntimeError(f"Protein hydrogen {atom.index} lacks a heavy neighbor")
        old_mass = system.getParticleMass(atom.index)
        delta = target - old_mass
        if delta <= 0 * unit.dalton:
            continue
        heavy_mass = system.getParticleMass(heavy)
        if heavy_mass - delta <= 1.0 * unit.dalton:
            raise RuntimeError(f"Invalid HMR mass transfer at atom {atom.index}")
        system.setParticleMass(atom.index, target)
        system.setParticleMass(heavy, heavy_mass - delta)
        changed += 1
        transferred += delta

    for name in ("minimized_solvated.pdb", "minimized_state.xml"):
        shutil.copy2(args.source / name, args.output / name)
    (args.output / "system.xml").write_text(XmlSerializer.serialize(system))
    metadata = {
        "source": str(args.source.resolve()),
        "source_system_sha256": sha256(args.source / "system.xml"),
        "derived_system_sha256": sha256(args.output / "system.xml"),
        "hydrogen_mass_amu": args.hydrogen_mass_amu,
        "protein_hydrogens_repartitioned": changed,
        "total_mass_transferred_amu": transferred.value_in_unit(unit.dalton),
        "potential_energy_model_changed": False,
        "water_hydrogens_repartitioned": False,
    }
    (args.output / "hmr_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
