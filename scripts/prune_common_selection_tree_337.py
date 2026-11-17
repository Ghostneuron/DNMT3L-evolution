from __future__ import annotations

import csv
import hashlib
from pathlib import Path

from Bio import Phylo


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "analysis" / "selection_reanalysis_337"
INPUT_DIR = ROOT / "inputs"
TREE_DIR = ROOT / "trees"
SOURCE_TREE = TREE_DIR / "Placentalia_CONSENSUS_common_MGK_F3X4_R4.treefile"
SUMMARY = PROJECT / "results" / "selection_337_tree_pruning_summary.tsv"


def fasta_names(path: Path) -> list[str]:
    names = [line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")]
    if len(names) != len(set(names)):
        raise ValueError(f"Duplicate FASTA identifiers in {path}")
    return names


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    if not SOURCE_TREE.exists():
        raise FileNotFoundError(SOURCE_TREE)

    source = Phylo.read(SOURCE_TREE, "newick")
    source_tips = [tip.name for tip in source.get_terminals()]
    if len(source_tips) != len(set(source_tips)):
        raise ValueError("The common source tree contains duplicate tip labels")

    rows: list[dict[str, object]] = []
    for alignment in sorted(INPUT_DIR.glob("Placentalia_*_human_coordinate_occupancy80_deduplicated.fasta")):
        label = alignment.name.split("_")[1]
        wanted = set(fasta_names(alignment))
        source_set = set(source_tips)
        missing = sorted(wanted - source_set)
        if missing:
            raise ValueError(f"{label}: alignment tips absent from source tree: {missing}")

        tree = Phylo.read(SOURCE_TREE, "newick")
        for tip in sorted(source_set - wanted):
            tree.prune(tip)

        observed = {tip.name for tip in tree.get_terminals()}
        if observed != wanted:
            raise ValueError(
                f"{label}: pruned tree mismatch; missing={sorted(wanted-observed)}, "
                f"extra={sorted(observed-wanted)}"
            )

        output = TREE_DIR / f"Placentalia_{label}_selection_common_topology.nwk"
        Phylo.write(tree, output, "newick")
        rows.append(
            {
                "analysis_id": label,
                "source_tree": str(SOURCE_TREE.relative_to(PROJECT)),
                "source_tip_count": len(source_tips),
                "selection_tip_count": len(wanted),
                "tips_pruned": len(source_tips) - len(wanted),
                "tip_set_verified": True,
                "selection_tree": str(output.relative_to(PROJECT)),
                "selection_tree_sha256": sha256(output),
            }
        )

    with SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"{row['analysis_id']}: {row['selection_tip_count']} tips, "
            f"verified={row['tip_set_verified']}"
        )


if __name__ == "__main__":
    main()
