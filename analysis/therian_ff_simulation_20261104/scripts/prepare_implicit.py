#!/usr/bin/env python3
"""Prepare one nonperiodic GBn2 screening system from a validated dry model."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from openmm import LangevinMiddleIntegrator, Platform, XmlSerializer, unit
from openmm.app import ForceField, HBonds, NoCutoff, PDBFile, Simulation


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def select_cpu(threads: int):
    platform = Platform.getPlatformByName("CPU")
    return platform, {"Threads": str(threads)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-pdb", required=True, type=Path)
    parser.add_argument("--variant", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--kappa-per-nm", type=float, default=1.0)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=False)
    pdb = PDBFile(str(args.dry_pdb))
    forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
        removeCMMotion=True,
        soluteDielectric=1.0,
        solventDielectric=78.5,
        implicitSolventKappa=args.kappa_per_nm / unit.nanometer,
    )
    integrator = LangevinMiddleIntegrator(
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    platform, properties = select_cpu(args.threads)
    simulation = Simulation(pdb.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.positions)
    initial = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy(maxIterations=5000)
    minimized = simulation.context.getState(getPositions=True, getEnergy=True)

    topology_path = args.output / "minimized_implicit.pdb"
    with topology_path.open("w") as handle:
        PDBFile.writeFile(
            pdb.topology, minimized.getPositions(), handle, keepIds=True
        )
    (args.output / "system.xml").write_text(XmlSerializer.serialize(system))
    (args.output / "minimized_state.xml").write_text(
        XmlSerializer.serialize(minimized)
    )
    metadata = {
        "variant": args.variant,
        "source_dry_pdb": str(args.dry_pdb.resolve()),
        "source_dry_pdb_sha256": sha256(args.dry_pdb),
        "topology_sha256": sha256(topology_path),
        "force_field": ["amber14-all.xml", "implicit/gbn2.xml"],
        "nonbonded_method": "NoCutoff",
        "solute_dielectric": 1.0,
        "solvent_dielectric": 78.5,
        "kappa_per_nm": args.kappa_per_nm,
        "atom_count": system.getNumParticles(),
        "initial_potential_energy_kj_mol": initial.value_in_unit(
            unit.kilojoule_per_mole
        ),
        "minimized_potential_energy_kj_mol": minimized.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        ),
        "platform": platform.getName(),
        "platform_properties": properties,
    }
    (args.output / "prepare_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
