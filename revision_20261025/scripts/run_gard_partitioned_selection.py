#!/usr/bin/env python3
"""Run FEL, MEME, and FUBAR on the completed GARD-partitioned Sauropsida alignment."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
import subprocess
import time


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--hyphy", default="hyphy")
    parser.add_argument("--libpath", type=Path, help="HyPhy resource directory for a noninstalled executable")
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("--mpi-workers", type=int, default=1)
    args = parser.parse_args()
    alignment = args.alignment.resolve()
    outdir = args.outdir.resolve()
    if not alignment.exists() or alignment.stat().st_size == 0:
        raise ValueError(f"GARD alignment snapshot is absent or empty: {alignment}")
    outdir.mkdir(parents=True, exist_ok=True)

    if args.mpi_workers < 1:
        raise ValueError("--mpi-workers must be at least one")
    prefix = (
        [args.mpirun, "-np", str(args.mpi_workers), args.hyphy]
        if args.mpi_workers > 1 else [args.hyphy]
    )
    if args.libpath:
        prefix.append(f"LIBPATH={args.libpath.resolve()}")

    commands = {
        "FUBAR": prefix + [
            "fubar", "--alignment", str(alignment), "--code", "Universal",
            "--grid", "20", "--method", "Variational-Bayes", "--concentration_parameter", "0.5",
            "--cache", str(outdir / "Sauropsida_MACSE_GARD.FUBAR.cache"),
            "--output", str(outdir / "Sauropsida_MACSE_GARD_FUBAR.json"),
        ],
        "FEL": prefix + [
            "fel", "--alignment", str(alignment), "--code", "Universal",
            "--branches", "All", "--srv", "Yes", "--multiple-hits", "None", "--pvalue", "0.1",
            "--output", str(outdir / "Sauropsida_MACSE_GARD_FEL.json"),
        ],
        "MEME": prefix + [
            "meme", "--alignment", str(alignment), "--code", "Universal",
            "--branches", "All", "--multiple-hits", "None", "--pvalue", "0.1",
            "--output", str(outdir / "Sauropsida_MACSE_GARD_MEME.json"),
        ],
    }
    manifest_rows = []
    for method, command in commands.items():
        output = Path(command[-1])
        if output.exists() and output.stat().st_size:
            raise FileExistsError(f"Refusing to overwrite completed {method} output: {output}")
        log = outdir / f"Sauropsida_MACSE_GARD_{method}.log"
        started = time.time()
        with log.open("w") as handle:
            completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, text=True, check=False)
        elapsed = time.time() - started
        if completed.returncode:
            raise subprocess.CalledProcessError(completed.returncode, command)
        if not output.exists() or output.stat().st_size == 0:
            raise ValueError(f"{method} completed without a nonempty output")
        manifest_rows.append({
            "Method": method,
            "Status": "complete",
            "Elapsed_Seconds": round(elapsed, 3),
            "Input": alignment.name,
            "Input_SHA256": sha256(alignment),
            "Output": output.name,
            "Output_SHA256": sha256(output),
            "Log": log.name,
            "Command": " ".join(command),
        })
    with (outdir / "gard_partitioned_selection_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_rows[0].keys())
        writer.writeheader()
        writer.writerows(manifest_rows)


if __name__ == "__main__":
    main()
