from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
import shlex
import shutil
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "analysis" / "selection_reanalysis_337"
INPUT_DIR = ROOT / "inputs"
TREE_DIR = ROOT / "trees"
RUN_DIR = ROOT / "runs"
LOG_DIR = PROJECT / "logs" / "selection_reanalysis_337"
MANIFEST = ROOT / "RUN_MANIFEST.tsv"
MANIFEST_LOCK = ROOT / ".run_manifest.lock"
FUNCTIONAL_POSITIONS = {88, 90, 105, 106, 107, 204, 349, 351, 352, 355, 358, 377}
CROSSWALK = PROJECT / "results" / "selection_337_human_coordinate_crosswalk.tsv"
FIELDS = [
    "run_id",
    "status",
    "analysis",
    "input",
    "input_sha256",
    "tree",
    "command",
    "started_at",
    "completed_at",
    "exit_code",
    "pid",
    "host",
    "log",
    "output",
    "output_sha256",
    "notes",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def relative(path: Path) -> str:
    return str(path.relative_to(PROJECT))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_manifest() -> list[dict[str, str]]:
    if not MANIFEST.exists():
        return []
    with MANIFEST.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_manifest(rows: list[dict[str, object]]) -> None:
    temporary = MANIFEST.with_name(f"{MANIFEST.name}.tmp.{os.getpid()}")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(MANIFEST)


def upsert(row: dict[str, object]) -> None:
    with MANIFEST_LOCK.open("a+") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        rows = read_manifest()
        rows = [existing for existing in rows if existing["run_id"] != row["run_id"]]
        rows.append(row)
        rows.sort(key=lambda item: item["run_id"])
        write_manifest(rows)
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def validate_json(path: Path) -> None:
    if not path.exists() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty JSON output: {path}")
    with path.open() as handle:
        payload = json.load(handle)
    if "input" not in payload or "data partitions" not in payload:
        raise ValueError(f"HyPhy JSON lacks expected completed-output keys: {path}")


def functional_site_pattern(analysis: str) -> str:
    with CROSSWALK.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sites = sorted(
        int(row["filtered_codon_site"])
        for row in rows
        if row["analysis_id"] == analysis and int(row["human_position"]) in FUNCTIONAL_POSITIONS
    )
    if len(sites) != len(FUNCTIONAL_POSITIONS):
        raise ValueError(
            f"{analysis}: expected {len(FUNCTIONAL_POSITIONS)} retained functional sites; found {len(sites)}"
        )
    return ",".join(str(site) for site in sites)


def command_for(model: str, analysis: str, alignment: Path, tree: Path, output: Path) -> list[str]:
    hyphy = shutil.which("hyphy")
    if not hyphy:
        raise RuntimeError("hyphy is not available on PATH")
    common = ["--alignment", str(alignment), "--tree", str(tree), "--code", "Universal"]
    if model == "FUBAR":
        return [
            hyphy,
            "fubar",
            *common,
            "--grid",
            "20",
            "--method",
            "Variational-Bayes",
            "--concentration_parameter",
            "0.5",
            "--output",
            str(output),
        ]
    if model == "MEME_STANDARD":
        return [
            hyphy,
            "meme",
            *common,
            "--branches",
            "All",
            "--pvalue",
            "0.1",
            "--multiple-hits",
            "None",
            "--output",
            str(output),
        ]
    if model == "MEME_STANDARD_FUNCTIONAL":
        return [
            hyphy,
            "meme",
            *common,
            "--branches",
            "All",
            "--pvalue",
            "0.1",
            "--multiple-hits",
            "None",
            "--limit-to-sites",
            functional_site_pattern(analysis),
            "--output",
            str(output),
        ]
    if model == "MEME_MULTIHIT_FUNCTIONAL":
        return [
            hyphy,
            "meme",
            *common,
            "--branches",
            "All",
            "--pvalue",
            "0.1",
            "--multiple-hits",
            "Double+Triple",
            "--site-multihit",
            "Estimate",
            "--limit-to-sites",
            functional_site_pattern(analysis),
            "--output",
            str(output),
        ]
    raise ValueError(model)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis", required=True, choices=["CONSENSUS", "MACSE", "MAFFT", "PRANK", "MUSCLE5", "COBALT"])
    parser.add_argument(
        "--model",
        required=True,
        choices=[
            "FUBAR",
            "MEME_STANDARD",
            "MEME_STANDARD_FUNCTIONAL",
            "MEME_MULTIHIT_FUNCTIONAL",
        ],
    )
    parser.add_argument("--tree-kind", choices=["COMMON", "NCBI_TAXONOMY"], default="COMMON")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    RUN_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    if args.tree_kind == "NCBI_TAXONOMY":
        if args.analysis != "CONSENSUS":
            raise ValueError("The NCBI taxonomy topology is defined for CONSENSUS only")
        alignment = INPUT_DIR / "Placentalia_CONSENSUS_human_coordinate_occupancy80_deduplicated_ncbi_taxonomy.fasta"
        tree = TREE_DIR / "Placentalia_CONSENSUS_selection_ncbi_taxonomy_topology.nwk"
    else:
        alignment = INPUT_DIR / f"Placentalia_{args.analysis}_human_coordinate_occupancy80_deduplicated.fasta"
        tree = TREE_DIR / f"Placentalia_{args.analysis}_selection_common_topology.nwk"
    suffix = "_NCBI_TAXONOMY" if args.tree_kind == "NCBI_TAXONOMY" else ""
    run_id = f"{args.analysis}_{args.model}{suffix}"
    run_path = RUN_DIR / run_id
    run_path.mkdir(parents=True, exist_ok=True)
    output = run_path / f"{run_id}.json"
    log = LOG_DIR / f"{run_id}.log"

    if not alignment.exists() or not tree.exists():
        raise FileNotFoundError(f"Required input missing: alignment={alignment.exists()}, tree={tree.exists()}")
    previous = {row["run_id"]: row for row in read_manifest()}.get(run_id)
    if previous and previous.get("status") == "COMPLETE" and not args.force:
        validate_json(output)
        print(f"{run_id} already complete and validated; use --force to rerun")
        return

    command = command_for(args.model, args.analysis, alignment, tree, output)
    command_text = shlex.join(command)
    row: dict[str, object] = {
        "run_id": run_id,
        "status": "STARTING",
        "analysis": args.analysis,
        "input": relative(alignment),
        "input_sha256": sha256(alignment),
        "tree": relative(tree),
        "command": command_text,
        "started_at": utc_now(),
        "completed_at": "",
        "exit_code": "",
        "pid": "",
        "host": socket.gethostname(),
        "log": relative(log),
        "output": relative(output),
        "output_sha256": "",
        "notes": (
            f"OMP_NUM_THREADS={args.threads}; tree_kind={args.tree_kind}; "
            "decision rules frozen before output inspection"
        ),
    }
    upsert(row)

    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(args.threads)
    with log.open("w") as handle:
        handle.write(f"STARTED_UTC\t{row['started_at']}\n")
        handle.write(f"HOST\t{row['host']}\n")
        handle.write(f"COMMAND\t{command_text}\n")
        handle.flush()
        process = subprocess.Popen(
            command,
            cwd=PROJECT,
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
        row["pid"] = process.pid
        row["status"] = "RUNNING"
        upsert(row)
        print(f"{run_id} running as PID {process.pid}; log={relative(log)}", flush=True)
        exit_code = process.wait()
        handle.write(f"\nEXIT_CODE\t{exit_code}\n")
        handle.write(f"COMPLETED_UTC\t{utc_now()}\n")

    row["exit_code"] = exit_code
    row["completed_at"] = utc_now()
    if exit_code == 0:
        try:
            validate_json(output)
            row["status"] = "COMPLETE"
            row["output_sha256"] = sha256(output)
        except Exception as error:
            row["status"] = "FAILED_VALIDATION"
            row["notes"] = f"{row['notes']}; validation error: {error}"
    else:
        row["status"] = "FAILED"
    upsert(row)
    print(f"{run_id}: {row['status']} (exit {exit_code})", flush=True)
    if row["status"] != "COMPLETE":
        sys.exit(1)


if __name__ == "__main__":
    main()
