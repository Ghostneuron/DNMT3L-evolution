#!/usr/bin/env python3
"""Run or resume one matched OpenMM replica from a prepared variant."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, XmlSerializer, unit
from openmm.app import (
    CheckpointReporter,
    DCDReporter,
    PDBFile,
    Simulation,
    StateDataReporter,
)


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
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--equilibration-steps", type=int, default=50_000)
    parser.add_argument("--production-steps", type=int, required=True)
    parser.add_argument("--report-interval", type=int, default=5_000)
    parser.add_argument("--checkpoint-interval", type=int, default=50_000)
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--platform", default="CPU")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    done_path = args.output / "COMPLETE.json"
    if done_path.exists():
        raise RuntimeError(f"Completed output already exists: {done_path}")

    pdb = PDBFile(str(args.prepared / "minimized_solvated.pdb"))
    system = XmlSerializer.deserialize((args.prepared / "system.xml").read_text())
    system.addForce(MonteCarloBarostat(1.0 * unit.bar, 300.0 * unit.kelvin, 25))
    integrator = LangevinMiddleIntegrator(
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        args.timestep_fs * unit.femtoseconds,
    )
    integrator.setRandomNumberSeed(args.seed)
    for force in system.getForces():
        if isinstance(force, MonteCarloBarostat):
            force.setRandomNumberSeed(args.seed + 1)
    platform, platform_properties = select_platform(args.platform, args.threads)
    simulation = Simulation(
        pdb.topology,
        system,
        integrator,
        platform,
        platform_properties,
    )

    checkpoint = args.output / "checkpoint.chk"
    resumed = checkpoint.exists()
    if resumed:
        simulation.loadCheckpoint(str(checkpoint))
    else:
        state = XmlSerializer.deserialize(
            (args.prepared / "minimized_state.xml").read_text()
        )
        simulation.context.setState(state)
        simulation.context.setVelocitiesToTemperature(300.0 * unit.kelvin, args.seed)

    simulation.reporters.append(
        StateDataReporter(
            str(args.output / "state.tsv"),
            args.report_interval,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True,
            remainingTime=True,
            totalSteps=args.equilibration_steps + args.production_steps,
            separator="\t",
            append=resumed,
        )
    )
    simulation.reporters.append(
        CheckpointReporter(str(checkpoint), args.checkpoint_interval)
    )

    total_target = args.equilibration_steps + args.production_steps
    if simulation.currentStep > total_target:
        raise RuntimeError(
            f"Checkpoint step {simulation.currentStep} exceeds target {total_target}"
        )

    started = time.time()
    equilibration_remaining = max(
        0, args.equilibration_steps - simulation.currentStep
    )
    if equilibration_remaining:
        simulation.step(equilibration_remaining)

    simulation.reporters.append(
        DCDReporter(
            str(args.output / "trajectory.dcd"),
            args.report_interval,
            append=resumed and (args.output / "trajectory.dcd").exists(),
        )
    )
    production_remaining = total_target - simulation.currentStep
    if production_remaining:
        simulation.step(production_remaining)

    final_state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getEnergy=True,
        enforcePeriodicBox=True,
    )
    (args.output / "final_state.xml").write_text(XmlSerializer.serialize(final_state))
    simulation.saveCheckpoint(str(checkpoint))
    runtime = time.time() - started
    completion = {
        "seed": args.seed,
        "equilibration_steps": args.equilibration_steps,
        "production_steps": args.production_steps,
        "report_interval": args.report_interval,
        "time_step_fs": args.timestep_fs,
        "production_time_ns": args.production_steps * args.timestep_fs / 1_000_000.0,
        "runtime_seconds_this_invocation": runtime,
        "resumed_from_checkpoint": resumed,
        "final_step": simulation.currentStep,
        "final_potential_energy_kj_mol": final_state.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        ),
        "platform": platform.getName(),
        "platform_properties": platform_properties,
    }
    done_path.write_text(json.dumps(completion, indent=2) + "\n")
    print(json.dumps(completion, indent=2))


if __name__ == "__main__":
    main()

