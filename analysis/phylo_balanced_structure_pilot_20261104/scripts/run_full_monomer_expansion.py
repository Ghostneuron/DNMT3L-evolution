#!/usr/bin/env python3
"""Resume-safe three-seed monomer expansion for all 24 balanced-panel inputs."""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import time
from pathlib import Path

from run_sentinel_predictions import archive_prediction, now, validate_prediction


SEEDS = [4101, 4102, 4103]


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_status(path, new_row, target_order):
    rows = read_tsv(path) if path.exists() else []
    key = (new_row["target"], new_row["seed"])
    rows = [row for row in rows if (row["target"], row["seed"]) != key]
    rows.append(new_row)
    rows.sort(key=lambda row: (target_order.index(row["target"]), int(row["seed"])))
    fields = ["target", "seed", "status", "source", "started_at", "finished_at", "elapsed_seconds", "log", "output"]
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(".tsv.partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-root", type=Path, required=True)
    parser.add_argument("--boltz", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()

    panel = read_tsv(args.panel_root / "metadata/panel_manifest.tsv")
    targets = [row["safe_id"] for row in panel]
    if len(targets) != 24 or len(set(targets)) != 24:
        raise ValueError("Expected exactly 24 unique balanced-panel targets")
    sentinel_status = read_tsv(args.panel_root / "metadata/sentinel_run_status.tsv")
    sentinel_monomers = {
        (row["target"], int(row["seed"])): args.panel_root / row["output"]
        for row in sentinel_status if row["mode"] == "monomer" and row["status"] == "VALIDATED"
    }
    status_path = args.panel_root / "metadata/full_monomer_run_status.tsv"
    env = os.environ.copy()
    env["BOLTZ_CACHE"] = str(args.cache)
    env["MPLCONFIGDIR"] = "${TMPDIR}/dnmt3l_mplconfig"

    for target in targets:
        for seed in SEEDS:
            reused = sentinel_monomers.get((target, seed))
            if reused is not None:
                validate_prediction(reused)
                write_status(status_path, {
                    "target": target, "seed": str(seed), "status": "VALIDATED_REUSED", "source": "sentinel",
                    "started_at": "", "finished_at": now(), "elapsed_seconds": "", "log": "",
                    "output": str(reused.relative_to(args.panel_root)),
                }, targets)
                continue

            destination = args.panel_root / f"outputs/full_monomer/{target}/seed{seed}"
            try:
                validate_prediction(destination)
                write_status(status_path, {
                    "target": target, "seed": str(seed), "status": "VALIDATED", "source": "full_monomer_expansion",
                    "started_at": "", "finished_at": now(), "elapsed_seconds": "", "log": "",
                    "output": str(destination.relative_to(args.panel_root)),
                }, targets)
                continue
            except (ValueError, FileNotFoundError):
                pass

            input_yaml = args.panel_root / f"inputs_local_msa/monomer_yaml/{target}.yaml"
            scratch_run = args.scratch / f"full_monomer_{target}_seed{seed}"
            log_path = args.panel_root / f"logs/full_monomer/{target}_seed{seed}.log"
            log_path.parent.mkdir(parents=True, exist_ok=True)
            command = [
                str(args.boltz), "predict", str(input_yaml), "--out_dir", str(scratch_run),
                "--accelerator", "mps", "--devices", "1", "--seed", str(seed), "--model", "boltz2",
                "--recycling_steps", "3", "--sampling_steps", "200", "--diffusion_samples", "1",
                "--max_parallel_samples", "1", "--output_format", "mmcif", "--write_full_pae",
                "--write_full_pde", "--num_workers", "0",
            ]
            started = now()
            start_clock = time.monotonic()
            write_status(status_path, {
                "target": target, "seed": str(seed), "status": "RUNNING", "source": "full_monomer_expansion",
                "started_at": started, "finished_at": "", "elapsed_seconds": "",
                "log": str(log_path.relative_to(args.panel_root)), "output": str(destination.relative_to(args.panel_root)),
            }, targets)
            with log_path.open("w") as log:
                log.write("COMMAND\t" + " ".join(command) + "\n")
                log.flush()
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
            elapsed = time.monotonic() - start_clock
            if result.returncode:
                write_status(status_path, {
                    "target": target, "seed": str(seed), "status": f"FAILED_{result.returncode}", "source": "full_monomer_expansion",
                    "started_at": started, "finished_at": now(), "elapsed_seconds": f"{elapsed:.1f}",
                    "log": str(log_path.relative_to(args.panel_root)), "output": str(destination.relative_to(args.panel_root)),
                }, targets)
                raise SystemExit(result.returncode)
            archive_prediction(scratch_run, destination, {
                "mode": "monomer", "target": target, "seed": seed,
                "software": "boltz-community==2.13.0; Boltz-2 checkpoint",
                "command": command, "started_at": started, "elapsed_seconds": elapsed,
                "origin": "balanced 24-input monomer expansion",
            })
            write_status(status_path, {
                "target": target, "seed": str(seed), "status": "VALIDATED", "source": "full_monomer_expansion",
                "started_at": started, "finished_at": now(), "elapsed_seconds": f"{elapsed:.1f}",
                "log": str(log_path.relative_to(args.panel_root)), "output": str(destination.relative_to(args.panel_root)),
            }, targets)
            print(f"VALIDATED monomer {target} seed {seed} in {elapsed:.1f}s", flush=True)


if __name__ == "__main__":
    main()
