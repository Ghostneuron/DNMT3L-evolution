from __future__ import annotations

import csv
import hashlib
import shutil
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "source_package" / "analysis" / "orthology_rebuild_20260826"
ORTHOLOGY = SOURCE / "orthology_audit"
OUT = PROJECT / "analysis_v2" / "inputs"

SCOPES = ("Placentalia", "Euarchontoglires", "Laurasiatheria", "Sauropsida")
ANCHOR_SCOPES = {"Laurasiatheria", "Sauropsida"}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks).upper()))
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    return records


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 90):
                handle.write(sequence[start : start + 90] + "\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    placental_cds = dict(
        read_fasta(ORTHOLOGY / "DNMT3L_selection_ready_Placentalia_cds.fasta")
    )
    placental_protein = dict(
        read_fasta(ORTHOLOGY / "DNMT3L_selection_ready_Placentalia_protein.fasta")
    )
    human_cds = placental_cds["Homo_sapiens"]
    human_protein = placental_protein["Homo_sapiens"]

    rows: list[dict[str, object]] = []
    for scope in SCOPES:
        for molecule in ("cds", "protein"):
            source = ORTHOLOGY / f"DNMT3L_selection_ready_{scope}_{molecule}.fasta"
            records = read_fasta(source)
            anchor_added = scope in ANCHOR_SCOPES
            if anchor_added:
                anchor_sequence = human_cds if molecule == "cds" else human_protein
                records.append(("Homo_sapiens_coordinate_anchor", anchor_sequence))
            destination = OUT / f"{scope}_{molecule}_with_coordinate_anchor.fasta"
            write_fasta(records, destination)
            rows.append(
                {
                    "Scope": scope,
                    "Molecule": molecule,
                    "Biological_Sequence_Count": len(records) - int(anchor_added),
                    "Coordinate_Anchor_Added": anchor_added,
                    "Alignment_Sequence_Count": len(records),
                    "Source_Path": str(source.relative_to(PROJECT)),
                    "Output_Path": str(destination.relative_to(PROJECT)),
                    "Output_SHA256": sha256(destination),
                }
            )

    manifest = OUT / "input_manifest.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    shutil.copy2(
        SOURCE / "clade_alignments" / "Homo_sapiens_anchor_DNMT3L_cds.fasta",
        OUT / "archived_human_coordinate_anchor_cds.fasta",
    )
    print(f"Prepared {len(rows)} inputs in {OUT}")


if __name__ == "__main__":
    main()
