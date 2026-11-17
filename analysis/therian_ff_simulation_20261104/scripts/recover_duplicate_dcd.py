#!/usr/bin/env python3
"""Rebuild one canonical replica directory after concurrent duplicate writers.

The complete source directory is preserved.  For each reported step, this
utility retains the first state row and its corresponding production DCD frame.
It refuses non-regular step sequences and never overwrites an output directory.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import mdtraj as md


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--checkpoint-step", required=True, type=int)
    parser.add_argument("--report-interval", required=True, type=int)
    parser.add_argument("--equilibration-steps", required=True, type=int)
    args = parser.parse_args()

    if args.output.exists():
        raise RuntimeError(f"Refusing to overwrite output: {args.output}")

    state_path = args.source / "state.tsv"
    lines = state_path.read_text().splitlines()
    headers = [line for line in lines if line.startswith("#")]
    if len(headers) != 1:
        raise RuntimeError(f"Expected one state header; found {len(headers)}")

    records: list[tuple[int, str]] = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        records.append((int(line.split("\t", 1)[0]), line))

    first_record: dict[int, str] = {}
    selected_dcd_indices: list[int] = []
    production_record_index = 0
    for step, line in records:
        is_first = step not in first_record
        if is_first:
            first_record[step] = line
        if step > args.equilibration_steps:
            if is_first:
                selected_dcd_indices.append(production_record_index)
            production_record_index += 1

    expected_steps = list(
        range(args.report_interval, args.checkpoint_step + 1, args.report_interval)
    )
    observed_steps = sorted(first_record)
    if observed_steps != expected_steps:
        raise RuntimeError(
            "Unique state steps do not form the expected regular sequence: "
            f"expected {len(expected_steps)}, observed {len(observed_steps)}"
        )

    trajectory = md.load_dcd(
        str(args.source / "trajectory.dcd"), top=str(args.topology)
    )
    if trajectory.n_frames != production_record_index:
        raise RuntimeError(
            f"DCD/state mismatch: {trajectory.n_frames} frames versus "
            f"{production_record_index} production state rows"
        )

    args.output.mkdir(parents=True)
    clean_state = [headers[0], *(first_record[step] for step in observed_steps)]
    (args.output / "state.tsv").write_text("\n".join(clean_state) + "\n")
    trajectory[selected_dcd_indices].save_dcd(str(args.output / "trajectory.dcd"))

    for name in (
        "checkpoint.chk",
        "checkpoint.step37500.pre_resume.chk",
        "checkpoint.step285000.pre_sanitize.chk",
        "driver.log",
    ):
        source = args.source / name
        if source.exists():
            shutil.copy2(source, args.output / name)

    provenance = {
        "source_directory": str(args.source.resolve()),
        "selection_rule": "first state row and corresponding DCD frame per step",
        "checkpoint_step": args.checkpoint_step,
        "report_interval": args.report_interval,
        "equilibration_steps": args.equilibration_steps,
        "source_state_rows": len(records),
        "canonical_state_rows": len(observed_steps),
        "removed_duplicate_state_rows": len(records) - len(observed_steps),
        "source_dcd_frames": trajectory.n_frames,
        "canonical_dcd_frames": len(selected_dcd_indices),
        "removed_duplicate_dcd_frames": (
            trajectory.n_frames - len(selected_dcd_indices)
        ),
    }
    (args.output / "RECOVERY_PROVENANCE.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
