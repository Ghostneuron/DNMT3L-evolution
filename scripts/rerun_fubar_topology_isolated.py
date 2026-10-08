from __future__ import annotations

import os

import csv
import hashlib
import json
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
V2 = PROJECT / "analysis_v2"
INPUTS = V2 / "selection_inputs"
GENE_TREES = V2 / "phylogenies"
SPECIES_TREES = V2 / "species_trees"
GENE_OUT = V2 / "hyphy"
SPECIES_OUT = V2 / "species_tree_hyphy"
ISOLATED = V2 / "fubar_isolated_inputs"
DEPRECATED = V2 / "deprecated" / "fubar_shared_cache_20261006"
LOGS = V2 / "logs" / "fubar_topology_isolated"
MAX_WORKERS = 5


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


def archive_initial_output(output: Path, topology: str, dataset: str) -> Path | None:
    if not output.exists():
        return None
    DEPRECATED.mkdir(parents=True, exist_ok=True)
    archived = DEPRECATED / f"{dataset}_{topology}_{output.name}"
    if not archived.exists():
        shutil.copy2(output, archived)
    return archived


def archive_shared_input_caches() -> Path:
    archive_dir = DEPRECATED / "selection_input_caches"
    archive_dir.mkdir(parents=True, exist_ok=True)
    manifest = archive_dir / "shared_cache_archive_manifest.csv"
    rows: list[dict[str, object]] = []
    for cache in sorted(INPUTS.glob("*.FUBAR.cache")):
        if cache.name.startswith("._"):
            continue
        source_hash = sha256(cache)
        archived = archive_dir / cache.name
        if archived.exists() and sha256(archived) != source_hash:
            archived = archive_dir / f"{cache.name}.{source_hash[:12]}"
        if not archived.exists():
            shutil.copy2(cache, archived)
        rows.append(
            {
                "Original_Cache": str(cache.relative_to(PROJECT)),
                "Original_Cache_SHA256": source_hash,
                "Archived_Cache": str(archived.relative_to(PROJECT)),
                "Archived_Cache_SHA256": sha256(archived),
            }
        )
        cache.unlink()
    if rows:
        with manifest.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    elif not manifest.exists():
        raise FileNotFoundError("No shared FUBAR caches were found and no prior archive manifest exists")
    return manifest


def refresh_parent_manifest(path: Path, topology_rows: list[dict[str, object]], provenance: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if "Superseded_Command" not in fieldnames:
        fieldnames.append("Superseded_Command")
    if "Final_Provenance" not in fieldnames:
        fieldnames.append("Final_Provenance")
    replacements = {str(row["Dataset"]): row for row in topology_rows}
    replaced: set[str] = set()
    for row in rows:
        dataset = str(row.get("Dataset", ""))
        if str(row.get("Method", "")).upper() != "FUBAR" or dataset not in replacements:
            row.setdefault("Superseded_Command", "")
            row.setdefault("Final_Provenance", "")
            continue
        final = replacements[dataset]
        row["Superseded_Command"] = row.get("Command", "")
        row["Status"] = "recomputed_topology_isolated"
        row["Elapsed_Seconds"] = str(final["Elapsed_Seconds"])
        row["Output_Path"] = str(final["Output"])
        row["Output_SHA256"] = str(final["Output_SHA256"])
        row["Command"] = str(final["Command"])
        row["Final_Provenance"] = str(provenance.relative_to(PROJECT))
        replaced.add(dataset)
    if replaced != set(replacements):
        missing = sorted(set(replacements) - replaced)
        raise ValueError(f"Could not refresh FUBAR rows in {path}: {missing}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_job(topology: str, dataset: str) -> dict[str, object]:
    original_alignment = INPUTS / f"{dataset}_selection_input.fasta"
    if topology == "gene_tree":
        tree = GENE_TREES / dataset / f"{dataset}.treefile"
        output = GENE_OUT / dataset / "FUBAR_all_sites_standard.json"
    else:
        tree = SPECIES_TREES / f"{dataset}_NCBI_taxonomy_tree.nwk"
        output = SPECIES_OUT / dataset / "FUBAR_species_tree.json"
    if not tree.exists():
        raise FileNotFoundError(tree)
    tree_hash = sha256(tree)
    isolated_dir = ISOLATED / topology / dataset / tree_hash[:12]
    isolated_dir.mkdir(parents=True, exist_ok=True)
    alignment = isolated_dir / f"{dataset}_{topology}_alignment.fasta"
    shutil.copy2(original_alignment, alignment)
    output.parent.mkdir(parents=True, exist_ok=True)
    archived = archive_initial_output(output, topology, dataset)

    command = [
        os.environ.get("HYPHY_BIN", "hyphy"),
        "fubar",
        "--alignment",
        str(alignment),
        "--tree",
        str(tree),
        "--code",
        "Universal",
        "--grid",
        "20",
        "--method",
        "Variational-Bayes",
        "--concentration_parameter",
        "0.5",
        "--output",
        str(output),
    ]
    LOGS.mkdir(parents=True, exist_ok=True)
    log = LOGS / f"{dataset}_{topology}.log"
    start = time.monotonic()
    with log.open("w") as handle:
        handle.write("COMMAND: " + " ".join(command) + "\n\n")
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, text=True, check=False)
    if completed.returncode:
        raise subprocess.CalledProcessError(completed.returncode, command)
    if not valid_json(output):
        raise ValueError(f"FUBAR did not produce valid JSON: {output}")
    cache = Path(str(alignment) + ".FUBAR.cache")
    if not cache.exists():
        raise FileNotFoundError(cache)
    return {
        "Topology": topology,
        "Dataset": dataset,
        "Elapsed_Seconds": round(time.monotonic() - start, 3),
        "Original_Alignment": str(original_alignment.relative_to(PROJECT)),
        "Original_Alignment_SHA256": sha256(original_alignment),
        "Isolated_Alignment": str(alignment.relative_to(PROJECT)),
        "Tree": str(tree.relative_to(PROJECT)),
        "Tree_SHA256": tree_hash,
        "Cache": str(cache.relative_to(PROJECT)),
        "Cache_SHA256": sha256(cache),
        "Output": str(output.relative_to(PROJECT)),
        "Output_SHA256": sha256(output),
        "Archived_Initial_Output": str(archived.relative_to(PROJECT)) if archived else "",
        "Archived_Initial_Output_SHA256": sha256(archived) if archived else "",
        "Command": " ".join(command),
    }


def main() -> None:
    datasets = sorted(
        path.name.removesuffix("_selection_input.fasta")
        for path in INPUTS.glob("*_selection_input.fasta")
        if not path.name.startswith("._")
    )
    if len(datasets) != 12:
        raise ValueError(f"Expected 12 datasets, found {len(datasets)}")
    species_datasets = [f"{scope}_MACSE" for scope in ("Placentalia", "Euarchontoglires", "Laurasiatheria", "Sauropsida")]
    jobs = [("gene_tree", dataset) for dataset in datasets] + [("species_tree", dataset) for dataset in species_datasets]

    shared_cache_manifest = archive_shared_input_caches()
    print(f"Archived shared-path caches under {shared_cache_manifest.parent.relative_to(PROJECT)}", flush=True)

    rows: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(run_job, topology, dataset): (topology, dataset) for topology, dataset in jobs}
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(f"{row['Topology']} {row['Dataset']}: {row['Elapsed_Seconds']} s", flush=True)
    rows.sort(key=lambda row: (str(row["Topology"]), str(row["Dataset"])))
    ISOLATED.mkdir(parents=True, exist_ok=True)
    manifest = ISOLATED / "fubar_topology_isolation_manifest.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    refresh_parent_manifest(
        GENE_OUT / "hyphy_manifest.csv",
        [row for row in rows if row["Topology"] == "gene_tree"],
        manifest,
    )
    refresh_parent_manifest(
        SPECIES_OUT / "species_tree_hyphy_manifest.csv",
        [row for row in rows if row["Topology"] == "species_tree"],
        manifest,
    )
    print(f"Wrote {manifest.relative_to(PROJECT)}")


if __name__ == "__main__":
    main()
