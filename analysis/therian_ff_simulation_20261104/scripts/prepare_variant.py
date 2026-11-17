#!/usr/bin/env python3
"""Prepare a minimized explicit-solvent DNMT3A-DNMT3L interface variant."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from pdbfixer import PDBFixer
from openmm import LangevinMiddleIntegrator, Platform, Vec3, XmlSerializer, unit
from openmm.app import (
    ForceField,
    HBonds,
    Modeller,
    PDBFile,
    PME,
    Simulation,
)


VARIANTS = {
    "amniote_YEE": ["PHE-261-TYR", "ASN-291-GLU"],
    "therian_FTD": ["ASN-291-THR", "GLU-293-ASP"],
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def residue_state(topology, chain_id: str, residue_ids: set[str]) -> dict[str, str]:
    states = {}
    for chain in topology.chains():
        if chain.id != chain_id:
            continue
        for residue in chain.residues():
            if residue.id in residue_ids:
                states[residue.id] = residue.name
    return states


def select_platform(requested: str, threads: int):
    available = {
        Platform.getPlatform(index).getName(): Platform.getPlatform(index)
        for index in range(Platform.getNumPlatforms())
    }
    if requested == "auto":
        name = next((x for x in ("OpenCL", "CUDA", "CPU") if x in available), None)
        if name is None:
            name = "Reference"
    else:
        name = requested
    if name not in available:
        raise RuntimeError(f"Requested platform {name} not in {sorted(available)}")
    properties = {}
    if name in {"OpenCL", "CUDA"}:
        properties["Precision"] = "mixed"
    elif name == "CPU":
        properties["Threads"] = str(threads)
    return available[name], properties


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--variant", required=True, choices=sorted(VARIANTS))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--padding-nm", type=float, default=1.0)
    parser.add_argument(
        "--box-nm",
        type=float,
        default=None,
        help="Use a fixed cubic box instead of mutation-dependent padding.",
    )
    parser.add_argument("--ionic-strength-m", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=20261104)
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=False)
    fixer = PDBFixer(filename=str(args.input))

    remove = [
        index
        for index, chain in enumerate(fixer.topology.chains())
        if chain.id not in {"A", "B"}
    ]
    fixer.removeChains(remove)
    fixer.removeHeterogens(False)
    fixer.applyMutations(VARIANTS[args.variant], "B")

    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues):
        chain_index, insertion_index = key
        residue_count = len(list(chains[chain_index].residues()))
        if insertion_index == 0 or insertion_index == residue_count:
            del fixer.missingResidues[key]

    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=args.seed)
    fixer.addMissingHydrogens(7.0)

    dry_states = residue_state(fixer.topology, "B", {"261", "291", "293"})
    expected = {
        "amniote_YEE": {"261": "TYR", "291": "GLU", "293": "GLU"},
        "therian_FTD": {"261": "PHE", "291": "THR", "293": "ASP"},
    }[args.variant]
    if dry_states != expected:
        raise RuntimeError(f"Mutation validation failed: {dry_states} != {expected}")

    dry_path = args.output / "prepared_dry.pdb"
    with dry_path.open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    modeller = Modeller(fixer.topology, fixer.positions)
    solvent_options = {
        "model": "tip3p",
        "ionicStrength": args.ionic_strength_m * unit.molar,
        "neutralize": True,
        "positiveIon": "Na+",
        "negativeIon": "Cl-",
    }
    if args.box_nm is None:
        solvent_options["padding"] = args.padding_nm * unit.nanometer
    else:
        solvent_options["boxSize"] = Vec3(
            args.box_nm, args.box_nm, args.box_nm
        ) * unit.nanometer
    modeller.addSolvent(forcefield, **solvent_options)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        removeCMMotion=True,
    )

    integrator = LangevinMiddleIntegrator(
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    platform, platform_properties = select_platform(args.platform, args.threads)
    simulation = Simulation(
        modeller.topology,
        system,
        integrator,
        platform,
        platform_properties,
    )
    simulation.context.setPositions(modeller.positions)
    initial_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy(maxIterations=2000)
    minimized = simulation.context.getState(getPositions=True, getEnergy=True)
    final_energy = minimized.getPotentialEnergy()

    solvated_path = args.output / "minimized_solvated.pdb"
    with solvated_path.open("w") as handle:
        PDBFile.writeFile(
            modeller.topology,
            minimized.getPositions(),
            handle,
            keepIds=True,
        )
    (args.output / "system.xml").write_text(XmlSerializer.serialize(system))
    (args.output / "minimized_state.xml").write_text(
        XmlSerializer.serialize(minimized)
    )

    metadata = {
        "variant": args.variant,
        "mutations_from_2QRV": VARIANTS[args.variant],
        "validated_dnmt3l_states": dry_states,
        "source": str(args.input.resolve()),
        "source_sha256": sha256(args.input),
        "dry_pdb_sha256": sha256(dry_path),
        "solvated_pdb_sha256": sha256(solvated_path),
        "force_field": ["amber14-all.xml", "amber14/tip3pfb.xml"],
        "water_model": "tip3p coordinates / TIP3P-FB parameters",
        "padding_nm": args.padding_nm,
        "fixed_cubic_box_nm": args.box_nm,
        "ionic_strength_m": args.ionic_strength_m,
        "atom_count": system.getNumParticles(),
        "force_count": system.getNumForces(),
        "initial_potential_energy_kj_mol": initial_energy.value_in_unit(
            unit.kilojoule_per_mole
        ),
        "minimized_potential_energy_kj_mol": final_energy.value_in_unit(
            unit.kilojoule_per_mole
        ),
        "platform": platform.getName(),
        "platform_properties": platform_properties,
    }
    (args.output / "prepare_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

