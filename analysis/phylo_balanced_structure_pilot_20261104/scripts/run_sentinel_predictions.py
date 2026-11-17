#!/usr/bin/env python3
"""Run and archive the resumable three-seed structure sentinel panel."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


TARGETS = [
    "Homo_sapiens",
    "ANCESTRAL_AMNIOTA",
    "Geotrypetes_seraphini",
    "Alligator_sinensis",
    "Monodelphis_domestica",
    "Mus_musculus",
    "ANCESTRAL_THERIA",
]
MODES = ["monomer", "complex"]
SEEDS = [4101, 4102, 4103]
SMOKE_SOURCES = {
    ("monomer", "Homo_sapiens", 4101): "smoke_local_msa_monomer_human_seed4101",
    ("monomer", "ANCESTRAL_AMNIOTA", 4101): "smoke_local_msa_monomer_amniota_seed4101",
    ("complex", "Homo_sapiens", 4101): "smoke_local_msa_complex_human_seed4101",
    ("complex", "ANCESTRAL_AMNIOTA", 4101): "smoke_local_msa_complex_amniota_seed4101",
}


def now():
    return datetime.now(timezone.utc).isoformat()


def atomic_text(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    partial.write_text(text)
    partial.replace(path)


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def prediction_dir(run_dir: Path, target: str):
    matches = list(run_dir.glob(f"boltz_results_*/predictions/{target}"))
    if len(matches) != 1:
        raise ValueError(f"Expected one prediction directory for {target} in {run_dir}, found {len(matches)}")
    return matches[0]


def validate_prediction(path: Path):
    def files(pattern):
        return [item for item in path.glob(pattern) if not item.name.startswith("._")]

    expected = {
        "cif": files("*.cif"),
        "confidence": files("confidence_*.json"),
        "pae": files("pae_*.npz"),
        "pde": files("pde_*.npz"),
        "plddt": files("plddt_*.npz"),
    }
    if any(len(files) != 1 for files in expected.values()):
        raise ValueError(f"Incomplete prediction artifacts: {path}")
    confidence = json.loads(expected["confidence"][0].read_text())
    if not confidence.get("structure_checks_passed"):
        raise ValueError(f"Structure checks failed: {path}")
    return confidence


def archive_prediction(source_run: Path, destination: Path, metadata: dict):
    source = prediction_dir(source_run, metadata["target"])
    confidence = validate_prediction(source)
    destination.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for source_file in source.iterdir():
        if not source_file.is_file() or source_file.name.startswith("._"):
            continue
        destination_file = destination / source_file.name
        shutil.copyfile(source_file, destination_file)
        hashes[source_file.name] = sha256(destination_file)
    validate_prediction(destination)
    metadata = {**metadata, "confidence": confidence, "output_sha256": hashes, "archived_at": now()}
    atomic_text(destination / "run_metadata.json", json.dumps(metadata, indent=2, sort_keys=True) + "\n")


def append_status(path: Path, row: dict):
    rows = []
    if path.exists():
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    key = (row["mode"], row["target"], row["seed"])
    rows = [old for old in rows if (old["mode"], old["target"], old["seed"]) != key]
    rows.append(row)
    rows.sort(key=lambda item: (MODES.index(item["mode"]), TARGETS.index(item["target"]), int(item["seed"])))
    fields = ["mode", "target", "seed", "status", "started_at", "finished_at", "elapsed_seconds", "log", "output"]
    partial = path.with_suffix(".tsv.partial")
    path.parent.mkdir(parents=True, exist_ok=True)
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

    status_path = args.panel_root / "metadata/sentinel_run_status.tsv"
    version = subprocess.run([str(args.boltz), "--help"], text=True, capture_output=True, check=False)
    software_identity = "boltz-community==2.13.0; Boltz-2 checkpoint"
    env = os.environ.copy()
    env["BOLTZ_CACHE"] = str(args.cache)
    env["MPLCONFIGDIR"] = "${TMPDIR}/dnmt3l_mplconfig"

    # Preserve all completed smoke artifacts before launching new work.
    for (mode, target, seed), smoke_name in SMOKE_SOURCES.items():
        destination = args.panel_root / f"outputs/sentinel/{mode}/{target}/seed{seed}"
        try:
            validate_prediction(destination)
            continue
        except (ValueError, FileNotFoundError):
            source_run = args.scratch / smoke_name
            if source_run.exists():
                archive_prediction(source_run, destination, {
                    "mode": mode, "target": target, "seed": seed,
                    "software": software_identity, "origin": "validated smoke run",
                })
                append_status(status_path, {
                    "mode": mode, "target": target, "seed": str(seed), "status": "VALIDATED",
                    "started_at": "", "finished_at": now(), "elapsed_seconds": "",
                    "log": "", "output": str(destination.relative_to(args.panel_root)),
                })

    for mode in MODES:
        for target in TARGETS:
            for seed in SEEDS:
                destination = args.panel_root / f"outputs/sentinel/{mode}/{target}/seed{seed}"
                try:
                    validate_prediction(destination)
                    append_status(status_path, {
                        "mode": mode, "target": target, "seed": str(seed), "status": "VALIDATED",
                        "started_at": "", "finished_at": now(), "elapsed_seconds": "",
                        "log": "", "output": str(destination.relative_to(args.panel_root)),
                    })
                    continue
                except (ValueError, FileNotFoundError):
                    pass

                smoke_name = SMOKE_SOURCES.get((mode, target, seed))
                if smoke_name:
                    source_run = args.scratch / smoke_name
                    if source_run.exists():
                        archive_prediction(source_run, destination, {
                            "mode": mode, "target": target, "seed": seed,
                            "software": software_identity, "origin": "validated smoke run",
                        })
                        append_status(status_path, {
                            "mode": mode, "target": target, "seed": str(seed), "status": "VALIDATED",
                            "started_at": "", "finished_at": now(), "elapsed_seconds": "",
                            "log": "", "output": str(destination.relative_to(args.panel_root)),
                        })
                        continue

                input_yaml = args.panel_root / f"inputs_local_msa/{mode}_yaml/{target}.yaml"
                scratch_run = args.scratch / f"sentinel_{mode}_{target}_seed{seed}"
                log_path = args.panel_root / f"logs/sentinel/{mode}_{target}_seed{seed}.log"
                log_path.parent.mkdir(parents=True, exist_ok=True)
                command = [
                    str(args.boltz), "predict", str(input_yaml),
                    "--out_dir", str(scratch_run), "--accelerator", "mps", "--devices", "1",
                    "--seed", str(seed), "--model", "boltz2", "--recycling_steps", "3",
                    "--sampling_steps", "200", "--diffusion_samples", "1",
                    "--max_parallel_samples", "1", "--output_format", "mmcif",
                    "--write_full_pae", "--write_full_pde", "--num_workers", "0",
                ]
                started = now()
                start_clock = time.monotonic()
                append_status(status_path, {
                    "mode": mode, "target": target, "seed": str(seed), "status": "RUNNING",
                    "started_at": started, "finished_at": "", "elapsed_seconds": "",
                    "log": str(log_path.relative_to(args.panel_root)), "output": str(destination.relative_to(args.panel_root)),
                })
                with log_path.open("w") as log:
                    log.write("COMMAND\t" + " ".join(command) + "\n")
                    log.flush()
                    result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
                elapsed = time.monotonic() - start_clock
                if result.returncode != 0:
                    append_status(status_path, {
                        "mode": mode, "target": target, "seed": str(seed), "status": f"FAILED_{result.returncode}",
                        "started_at": started, "finished_at": now(), "elapsed_seconds": f"{elapsed:.1f}",
                        "log": str(log_path.relative_to(args.panel_root)), "output": str(destination.relative_to(args.panel_root)),
                    })
                    raise SystemExit(result.returncode)
                archive_prediction(scratch_run, destination, {
                    "mode": mode, "target": target, "seed": seed, "software": software_identity,
                    "command": command, "started_at": started, "elapsed_seconds": elapsed,
                    "origin": "sentinel runner",
                })
                append_status(status_path, {
                    "mode": mode, "target": target, "seed": str(seed), "status": "VALIDATED",
                    "started_at": started, "finished_at": now(), "elapsed_seconds": f"{elapsed:.1f}",
                    "log": str(log_path.relative_to(args.panel_root)), "output": str(destination.relative_to(args.panel_root)),
                })
                print(f"VALIDATED {mode} {target} seed {seed} in {elapsed:.1f}s", flush=True)


if __name__ == "__main__":
    main()
