from __future__ import annotations

import os

import csv
import hashlib
import json
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
INPUTS = PROJECT / "analysis_v2" / "selection_inputs"
TREES = PROJECT / "analysis_v2" / "species_trees"
OUT = PROJECT / "analysis_v2" / "species_tree_hyphy"
LOGS = PROJECT / "analysis_v2" / "logs" / "species_tree_hyphy"
DATASETS = tuple(f"{scope}_MACSE" for scope in ("Placentalia", "Euarchontoglires", "Laurasiatheria", "Sauropsida"))
METHODS = ("fubar", "fel", "meme", "busted")
MAX_WORKERS = 3


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def valid_json(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        json.loads(path.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError):
        return False
    return True


def command_for(method: str, alignment: Path, tree: Path, output: Path) -> list[str]:
    base = [
        os.environ.get("HYPHY_BIN", "hyphy"),
        method,
        "--alignment",
        str(alignment),
        "--tree",
        str(tree),
        "--code",
        "Universal",
    ]
    if method == "fubar":
        return base + ["--grid", "20", "--method", "Variational-Bayes", "--concentration_parameter", "0.5", "--output", str(output)]
    if method in {"fel", "meme"}:
        return base + ["--branches", "All", "--pvalue", "0.1", "--multiple-hits", "None", "--output", str(output)]
    return base + ["--branches", "All", "--multiple-hits", "None", "--output", str(output)]


def run_job(dataset: str, method: str) -> dict[str, object]:
    alignment = INPUTS / f"{dataset}_selection_input.fasta"
    tree = TREES / f"{dataset}_NCBI_taxonomy_tree.nwk"
    output_dir = OUT / dataset
    output_dir.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{method.upper()}_species_tree.json"
    command = command_for(method, alignment, tree, output)
    start = time.monotonic()
    status = "reused_existing"
    if not valid_json(output):
        status = "newly_computed"
        with (LOGS / f"{dataset}_{method}.log").open("w") as log_handle:
            log_handle.write("COMMAND: " + " ".join(command) + "\n\n")
            log_handle.flush()
            completed = subprocess.run(command, stdout=log_handle, stderr=subprocess.STDOUT, text=True, check=False)
        if completed.returncode:
            raise subprocess.CalledProcessError(completed.returncode, command)
        if not valid_json(output):
            raise ValueError(f"Invalid HyPhy JSON: {output}")
    return {
        "Dataset": dataset,
        "Method": method.upper(),
        "Status": status,
        "Elapsed_Seconds": round(time.monotonic() - start, 3),
        "Output_Path": str(output.relative_to(PROJECT)),
        "Output_SHA256": sha256(output),
        "Command": " ".join(command),
    }


def main() -> None:
    jobs = [(dataset, method) for dataset in DATASETS for method in METHODS]
    rows: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(run_job, dataset, method): (dataset, method) for dataset, method in jobs}
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(f"{row['Dataset']} {row['Method']}: {row['Status']} ({row['Elapsed_Seconds']} s)", flush=True)
    rows.sort(key=lambda row: (str(row["Dataset"]), str(row["Method"])))
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "species_tree_hyphy_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
