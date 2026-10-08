from __future__ import annotations

import os

import csv
import hashlib
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from Bio import Phylo


PROJECT = Path(__file__).resolve().parents[1]
INPUTS = PROJECT / "analysis_v2" / "selection_inputs"
OUT = PROJECT / "analysis_v2" / "phylogenies"
LOGS = PROJECT / "analysis_v2" / "logs" / "phylogenies"
MAX_WORKERS = 3


def read_names(path: Path) -> list[str]:
    return [line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_tree(alignment: Path) -> dict[str, object]:
    label = alignment.name.removesuffix("_selection_input.fasta")
    output_dir = OUT / label
    output_dir.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / label
    tree_path = Path(f"{prefix}.treefile")
    command = [
        os.environ.get("IQTREE2_BIN", "iqtree2"),
        "-s",
        str(alignment),
        "-st",
        "CODON",
        "-m",
        "MGK+F3X4+R4",
        "-fast",
        "--alrt",
        "1000",
        "-T",
        "2",
        "--seed",
        "20261006",
        "--prefix",
        str(prefix),
    ]
    start = time.monotonic()
    status = "reused_existing"
    if not tree_path.exists():
        status = "newly_inferred"
        with (LOGS / f"{label}.log").open("w") as log_handle:
            log_handle.write("COMMAND: " + " ".join(command) + "\n\n")
            completed = subprocess.run(
                command,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
        if completed.returncode:
            raise subprocess.CalledProcessError(completed.returncode, command)

    alignment_names = set(read_names(alignment))
    tree = Phylo.read(tree_path, "newick")
    tree_names = {tip.name for tip in tree.get_terminals()}
    if alignment_names != tree_names:
        raise ValueError(
            f"Tip mismatch for {label}: missing={sorted(alignment_names-tree_names)}, extra={sorted(tree_names-alignment_names)}"
        )
    return {
        "Dataset": label,
        "Sequence_Count": len(alignment_names),
        "Alignment_Path": str(alignment.relative_to(PROJECT)),
        "Tree_Path": str(tree_path.relative_to(PROJECT)),
        "Tree_SHA256": sha256(tree_path),
        "Status": status,
        "Elapsed_Seconds": round(time.monotonic() - start, 3),
        "Command": " ".join(command),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    alignments = sorted(path for path in INPUTS.glob("*_selection_input.fasta") if not path.name.startswith("._"))
    if len(alignments) != 12:
        raise ValueError(f"Expected 12 selection inputs, found {len(alignments)}")
    rows: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(build_tree, alignment): alignment for alignment in alignments}
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(f"{row['Dataset']}: {row['Status']} ({row['Elapsed_Seconds']} s)", flush=True)
    rows.sort(key=lambda row: str(row["Dataset"]))
    with (OUT / "tree_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
