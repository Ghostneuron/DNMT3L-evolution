#!/usr/bin/env python3
"""Run one GBn2 screening replica from a prepared ancestral endpoint."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from openmm import LangevinMiddleIntegrator, Platform, XmlSerializer, unit
from openmm.app import CheckpointReporter, DCDReporter, PDBFile, Simulation, StateDataReporter


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--equilibration-steps", type=int, default=50_000)
    parser.add_argument("--production-steps", required=True, type=int)
    parser.add_argument("--report-interval", type=int, default=5_000)
    parser.add_argument("--checkpoint-interval", type=int, default=50_000)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    done_path = args.output / "COMPLETE.json"
    if done_path.exists():
        raise RuntimeError(f"Completed output already exists: {done_path}")

    pdb = PDBFile(str(args.prepared / "minimized_implicit.pdb"))
    system = XmlSerializer.deserialize((args.prepared / "system.xml").read_text())
    integrator = LangevinMiddleIntegrator(
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    integrator.setRandomNumberSeed(args.seed)
    platform = Platform.getPlatformByName("CPU")
    properties = {"Threads": str(args.threads)}
    simulation = Simulation(pdb.topology, system, integrator, platform, properties)
    checkpoint = args.output / "checkpoint.chk"
    if checkpoint.exists():
        raise RuntimeError(
            "Checkpoint exists without COMPLETE.json; automatic resume is disabled "
            "to prevent duplicate production steps. Start a new output directory."
        )
    state = XmlSerializer.deserialize((args.prepared / "minimized_state.xml").read_text())
    simulation.context.setState(state)
    simulation.context.setVelocitiesToTemperature(300.0 * unit.kelvin, args.seed)

    total_steps = args.equilibration_steps + args.production_steps
    simulation.reporters.append(
        StateDataReporter(
            str(args.output / "state.tsv"),
            args.report_interval,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            temperature=True,
            speed=True,
            remainingTime=True,
            totalSteps=total_steps,
            separator="\t",
        )
    )
    simulation.reporters.append(
        CheckpointReporter(str(checkpoint), args.checkpoint_interval)
    )

    started = time.time()
    if args.equilibration_steps:
        simulation.step(args.equilibration_steps)
    simulation.reporters.append(
        DCDReporter(str(args.output / "trajectory.dcd"), args.report_interval)
    )
    simulation.step(args.production_steps)
    final_state = simulation.context.getState(
        getPositions=True, getVelocities=True, getEnergy=True
    )
    (args.output / "final_state.xml").write_text(XmlSerializer.serialize(final_state))
    simulation.saveCheckpoint(str(checkpoint))
    completion = {
        "seed": args.seed,
        "equilibration_steps": args.equilibration_steps,
        "production_steps": args.production_steps,
        "report_interval": args.report_interval,
        "time_step_fs": 2.0,
        "production_time_ns": args.production_steps * 0.002 / 1000.0,
        "runtime_seconds": time.time() - started,
        "final_potential_energy_kj_mol": final_state.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        ),
        "platform": "CPU",
        "platform_properties": properties,
    }
    done_path.write_text(json.dumps(completion, indent=2) + "\n")
    print(json.dumps(completion, indent=2))


if __name__ == "__main__":
    main()
