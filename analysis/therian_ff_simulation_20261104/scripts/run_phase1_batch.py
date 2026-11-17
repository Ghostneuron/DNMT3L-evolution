#!/usr/bin/env python3
"""Run the prespecified Phase-1 explicit-solvent replica batch in two stages."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    temporary.replace(path)


def command_for(
    python: Path,
    runner: Path,
    root: Path,
    variant: str,
    seed: int,
    threads: int,
) -> list[str]:
    return [
        str(python),
        str(runner),
        "--prepared",
        str(root / "prepared_hmr" / variant),
        "--output",
        str(root / "production_hmr" / f"{variant}_seed{seed}"),
        "--seed",
        str(seed),
        "--equilibration-steps",
        "25000",
        "--production-steps",
        "1250000",
        "--report-interval",
        "2500",
        "--checkpoint-interval",
        "2500",
        "--platform",
        "CPU",
        "--threads",
        str(threads),
        "--timestep-fs",
        "4.0",
    ]


def run_stage(stage_name: str, jobs: list[dict], status: dict, status_path: Path) -> None:
    processes: list[tuple[dict, subprocess.Popen, object]] = []
    status["stage"] = stage_name
    status["updated_utc"] = now()
    atomic_json(status_path, status)
    for job in jobs:
        output = Path(job["output"])
        output.mkdir(parents=True, exist_ok=True)
        if (output / "COMPLETE.json").exists():
            job["status"] = "already_complete"
            continue
        log_handle = (output / "driver.log").open("a")
        log_handle.write(f"\n[{now()}] launch: {' '.join(job['command'])}\n")
        log_handle.flush()
        process = subprocess.Popen(
            job["command"],
            stdin=subprocess.DEVNULL,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            cwd=status["root"],
        )
        job.update({"status": "running", "pid": process.pid, "launched_utc": now()})
        processes.append((job, process, log_handle))
    status["updated_utc"] = now()
    atomic_json(status_path, status)

    while processes:
        remaining = []
        for job, process, log_handle in processes:
            code = process.poll()
            if code is None:
                remaining.append((job, process, log_handle))
                continue
            log_handle.write(f"[{now()}] exit code: {code}\n")
            log_handle.close()
            job["exit_code"] = code
            job["finished_utc"] = now()
            job["status"] = "complete" if code == 0 else "failed"
        processes = remaining
        status["updated_utc"] = now()
        atomic_json(status_path, status)
        if processes:
            time.sleep(60)
    failures = [job for job in jobs if job.get("status") == "failed"]
    if failures:
        raise RuntimeError(f"Stage {stage_name} failed: {[x['id'] for x in failures]}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    runner = root / "scripts" / "run_replica.py"
    status_path = root / "production_hmr" / "BATCH_STATUS.json"
    status_path.parent.mkdir(parents=True, exist_ok=True)

    specs = [
        ("stage1_four_way", "amniote_YEE", 1101, 2),
        ("stage1_four_way", "therian_FTD", 2101, 2),
        ("stage1_four_way", "amniote_YEE", 1102, 2),
        ("stage1_four_way", "therian_FTD", 2102, 2),
        ("stage2_paired", "amniote_YEE", 1103, 5),
        ("stage2_paired", "therian_FTD", 2103, 5),
    ]
    jobs = []
    for stage, variant, seed, threads in specs:
        command = command_for(args.python.resolve(), runner, root, variant, seed, threads)
        jobs.append(
            {
                "id": f"{variant}_seed{seed}",
                "stage": stage,
                "variant": variant,
                "seed": seed,
                "threads": threads,
                "output": str(root / "production_hmr" / f"{variant}_seed{seed}"),
                "command": command,
                "status": "pending",
            }
        )
    status = {
        "analysis_id": "DNMT3L_FF_THERIAN_MD_20261104_V1_PHASE1",
        "driver_pid": os.getpid(),
        "driver_started_utc": now(),
        "updated_utc": now(),
        "root": str(root),
        "python": str(args.python.resolve()),
        "equilibration_ps_per_replica": 100.0,
        "production_ns_per_replica": 5.0,
        "replicas_per_variant": 3,
        "jobs": jobs,
        "stage": "initializing",
        "status": "running",
    }
    atomic_json(status_path, status)
    try:
        run_stage(
            "stage1_four_way",
            [x for x in jobs if x["stage"] == "stage1_four_way"],
            status,
            status_path,
        )
        run_stage(
            "stage2_paired",
            [x for x in jobs if x["stage"] == "stage2_paired"],
            status,
            status_path,
        )
    except Exception as error:
        status.update({"status": "failed", "error": repr(error), "updated_utc": now()})
        atomic_json(status_path, status)
        raise
    status.update({"status": "complete", "stage": "complete", "updated_utc": now()})
    atomic_json(status_path, status)
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
