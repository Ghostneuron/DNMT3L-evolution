#!/usr/bin/env python3
"""Resumably download and validate the admitted Bgee germline H5AD panel."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def stamp():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, value):
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(path: Path):
    subprocess.run(["h5dump", "-H", str(path)], stdout=subprocess.DEVNULL, check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--checksum-dir", type=Path, required=True)
    parser.add_argument("--status", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.checksum_dir.mkdir(parents=True, exist_ok=True)
    args.status.parent.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(args.output_dir).free < 2 * 1024**3:
        raise SystemExit("Fewer than 2 GiB free; refusing panel download")

    with args.manifest.open(newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    state = {
        "analysis": "Murat_Bgee_germline_H5AD_panel",
        "started_utc": stamp(),
        "updated_utc": stamp(),
        "manifest": str(args.manifest.resolve()),
        "jobs": [],
    }
    for record in records:
        state["jobs"].append({
            "species": record["species"],
            "filename": record["filename"],
            "url": record["url"],
            "status": "pending",
        })
    atomic_json(args.status, state)

    for job in state["jobs"]:
        output = args.output_dir / job["filename"]
        partial = output.with_suffix(output.suffix + ".partial")
        checksum = args.checksum_dir / f"{output.stem}.sha256"
        try:
            if output.exists():
                validate(output)
                job["status"] = "validated_existing"
            else:
                job["status"] = "downloading"
                state["updated_utc"] = stamp()
                atomic_json(args.status, state)
                subprocess.run([
                    "curl", "--fail", "--location", "--retry", "4", "--retry-delay", "5",
                    "--continue-at", "-", "--output", str(partial), job["url"],
                ], check=True)
                validate(partial)
                partial.replace(output)
                job["status"] = "downloaded_validated"
            job["size_bytes"] = output.stat().st_size
            job["sha256"] = sha256(output)
            temporary_checksum = checksum.with_suffix(checksum.suffix + ".partial")
            temporary_checksum.write_text(f"{job['sha256']}  {output.name}\n")
            temporary_checksum.replace(checksum)
            job["completed_utc"] = stamp()
        except Exception as exc:
            job["status"] = "failed"
            job["error"] = repr(exc)
            state["updated_utc"] = stamp()
            atomic_json(args.status, state)
            raise
        state["updated_utc"] = stamp()
        atomic_json(args.status, state)
        print(f"{job['species']}: {job['status']} {job['size_bytes']} bytes", flush=True)

    state["status"] = "complete"
    state["completed_utc"] = stamp()
    state["updated_utc"] = stamp()
    atomic_json(args.status, state)


if __name__ == "__main__":
    main()
