from __future__ import annotations

import os

import csv
import subprocess
from pathlib import Path

from Bio import Phylo


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "source_package" / "analysis" / "orthology_rebuild_20260826" / "gene_family_tree"
OUT = PROJECT / "analysis_v2" / "orthology_sensitivity"
REFERENCE = "DNMT3L__Homo_sapiens"


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
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def exact_split(tree, target: set[str]) -> tuple[bool, float | None]:
    all_tips = {tip.name for tip in tree.get_terminals()}
    for clade in tree.find_clades(order="level"):
        descendants = {tip.name for tip in clade.get_terminals()}
        if descendants == target or all_tips - descendants == target:
            return True, clade.confidence
    return False, None


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    family = read_fasta(SOURCE / "DNMT3_family_candidates_plus_paralog_references.fasta")
    dnmt3c_records = read_fasta(OUT / "Mouse_DNMT3C_P0DOY1.fasta")
    if len(dnmt3c_records) != 1:
        raise ValueError("Expected one downloaded DNMT3C record")
    family.append(("DNMT3C__Mus_musculus__P0DOY1", dnmt3c_records[0][1]))
    combined = OUT / "DNMT3_family_plus_DNMT3C.fasta"
    write_fasta(family, combined)

    aligned = OUT / "DNMT3_family_plus_DNMT3C_MAFFT.fasta"
    with aligned.open("w") as output, (OUT / "mafft.log").open("w") as log:
        command = [os.environ.get("MAFFT_BIN", "mafft"), "--auto", "--thread", "8", str(combined)]
        completed = subprocess.run(command, stdout=output, stderr=log, text=True, check=False)
    if completed.returncode:
        raise subprocess.CalledProcessError(completed.returncode, command)

    records = read_fasta(aligned)
    by_name = dict(records)
    reference = by_name[REFERENCE]
    retained = [index for index, residue in enumerate(reference) if residue != "-"]
    if len(retained) != 387:
        raise ValueError(f"Expected 387 reference positions, found {len(retained)}")
    trimmed_records = [(name, "".join(sequence[index] for index in retained)) for name, sequence in records]
    trimmed = OUT / "DNMT3_family_plus_DNMT3C_human_columns.fasta"
    write_fasta(trimmed_records, trimmed)

    prefix = OUT / "DNMT3_family_plus_DNMT3C_LG_F_R4"
    command = [
        os.environ.get("IQTREE2_BIN", "iqtree2"),
        "-s",
        str(trimmed),
        "-m",
        "LG+F+R4",
        "-fast",
        "--alrt",
        "1000",
        "-T",
        "4",
        "--seed",
        "20261006",
        "--prefix",
        str(prefix),
    ]
    with (OUT / "iqtree.log").open("w") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True, check=False)
    if completed.returncode:
        raise subprocess.CalledProcessError(completed.returncode, command)

    tree_path = Path(f"{prefix}.treefile")
    tree = Phylo.read(tree_path, "newick")
    tips = [tip.name for tip in tree.get_terminals()]
    groups = {
        group: {tip for tip in tips if tip.startswith(f"{group}__")}
        for group in ("DNMT3L", "DNMT3A", "DNMT3B", "DNMT3C")
    }
    non_l = groups["DNMT3A"] | groups["DNMT3B"] | groups["DNMT3C"]
    summary_rows: list[dict[str, object]] = []
    for label, members in (
        ("DNMT3L", groups["DNMT3L"]),
        ("DNMT3A", groups["DNMT3A"]),
        ("DNMT3B_plus_DNMT3C", groups["DNMT3B"] | groups["DNMT3C"]),
        ("DNMT3A_plus_DNMT3B_plus_DNMT3C", non_l),
    ):
        separated, support = exact_split(tree, members)
        summary_rows.append(
            {
                "Group": label,
                "N_Tips": len(members),
                "Exact_Root_Independent_Split": separated,
                "SH_aLRT_Support": support,
            }
        )
    with (OUT / "group_separation_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    distance_rows: list[dict[str, object]] = []
    for candidate in sorted(groups["DNMT3L"]):
        nearest_l = min(
            ((tree.distance(candidate, other), other) for other in groups["DNMT3L"] if other != candidate),
            default=(float("inf"), None),
        )
        nearest_non_l = min(
            ((tree.distance(candidate, other), other) for other in non_l),
            default=(float("inf"), None),
        )
        distance_rows.append(
            {
                "Candidate_Tip": candidate,
                "Nearest_DNMT3L_Tip": nearest_l[1],
                "Nearest_DNMT3L_Distance": nearest_l[0],
                "Nearest_DNMT3A_B_or_C_Tip": nearest_non_l[1],
                "Nearest_DNMT3A_B_or_C_Distance": nearest_non_l[0],
                "Nearest_Group_Is_DNMT3L": nearest_l[0] < nearest_non_l[0],
                "Paralog_to_DNMT3L_Distance_Ratio": nearest_non_l[0] / nearest_l[0] if nearest_l[0] > 0 else float("inf"),
            }
        )
    distance_rows.sort(key=lambda row: float(row["Paralog_to_DNMT3L_Distance_Ratio"]))
    with (OUT / "candidate_distance_audit.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(distance_rows[0]))
        writer.writeheader()
        writer.writerows(distance_rows)

    print("Group separation")
    for row in summary_rows:
        print(row)
    failures = [row for row in distance_rows if not bool(row["Nearest_Group_Is_DNMT3L"])]
    print(f"Candidates nearest to a non-DNMT3L reference: {len(failures)}")
    print(f"DNMT3C nearest tip: {min((tree.distance(next(iter(groups['DNMT3C'])), tip), tip) for tip in groups['DNMT3B'])}")


if __name__ == "__main__":
    main()
