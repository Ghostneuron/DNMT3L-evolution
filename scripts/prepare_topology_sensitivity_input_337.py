from __future__ import annotations

import csv
import hashlib
import shutil
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
INPUTS = PROJECT / "analysis" / "selection_reanalysis_337" / "inputs"
SOURCE = INPUTS / "Placentalia_CONSENSUS_human_coordinate_occupancy80_deduplicated.fasta"
OUTPUT = INPUTS / "Placentalia_CONSENSUS_human_coordinate_occupancy80_deduplicated_ncbi_taxonomy.fasta"
MANIFEST = PROJECT / "results" / "selection_337_topology_sensitivity_input.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    shutil.copyfile(SOURCE, OUTPUT)
    source_hash = sha256(SOURCE)
    output_hash = sha256(OUTPUT)
    if source_hash != output_hash:
        raise ValueError("Topology-sensitivity input copy is not byte-identical")
    row = {
        "source": str(SOURCE.relative_to(PROJECT)),
        "source_sha256": source_hash,
        "isolated_input": str(OUTPUT.relative_to(PROJECT)),
        "isolated_input_sha256": output_hash,
        "reason": "unique filename prevents FUBAR cache reuse across tree topologies",
    }
    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    print(row)


if __name__ == "__main__":
    main()
