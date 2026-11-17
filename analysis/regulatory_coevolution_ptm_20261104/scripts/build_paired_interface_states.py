#!/usr/bin/env python3
"""Build a species-matched DNMT3A-DNMT3L interface-state inventory."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


DNMT3A_SITES = (725, 726, 729, 732, 772)
DNMT3L_SITES = (261, 291, 293)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def translate_codon(codon: str) -> str:
    codon = codon.upper().replace("!", "N")
    if len(codon) != 3 or "-" in codon or any(base not in "ACGT" for base in codon):
        return "NA"
    return str(Seq(codon).translate())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dnmt3a-alignment", type=Path, required=True)
    parser.add_argument("--dnmt3a-crosswalk", type=Path, required=True)
    parser.add_argument(
        "--dnmt3a-column-field",
        default="filtered_codon_site_1based",
        choices=("filtered_codon_site_1based", "original_macse_codon_1based"),
    )
    parser.add_argument("--dnmt3l-alignment", type=Path, required=True)
    parser.add_argument("--dnmt3l-matrix", type=Path, required=True)
    parser.add_argument("--states-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    args = parser.parse_args()

    dnmt3a = {record.id: str(record.seq).upper() for record in SeqIO.parse(args.dnmt3a_alignment, "fasta")}
    dnmt3l = {record.id: str(record.seq).upper() for record in SeqIO.parse(args.dnmt3l_alignment, "fasta")}
    crosswalk = {
        int(row["human_DNMT3A1_aa"]): int(row[args.dnmt3a_column_field])
        for row in read_tsv(args.dnmt3a_crosswalk)
        if row["human_DNMT3A1_aa"] and int(row["human_DNMT3A1_aa"]) in DNMT3A_SITES
    }
    if set(crosswalk) != set(DNMT3A_SITES):
        raise SystemExit(f"Incomplete DNMT3A crosswalk: {crosswalk}")

    metadata = {row["safe_id"]: row for row in read_tsv(args.dnmt3l_matrix)}
    matched = sorted(set(dnmt3a) & set(dnmt3l) & set(metadata))
    rows = []
    for safe_id in matched:
        meta = metadata[safe_id]
        if meta.get("Selection_Ready") != "True":
            continue
        row: dict[str, object] = {
            "safe_id": safe_id,
            "species": meta["species"],
            "clade": meta["clade"],
            "dnmt3a_alignment": args.dnmt3a_alignment.name,
            "dnmt3l_alignment": args.dnmt3l_alignment.name,
        }
        for position in DNMT3A_SITES:
            column = crosswalk[position]
            codon = dnmt3a[safe_id][(column - 1) * 3 : column * 3]
            row[f"DNMT3A1_{position}_codon"] = codon
            row[f"DNMT3A1_{position}_state"] = translate_codon(codon)
        for position in DNMT3L_SITES:
            state = dnmt3l[safe_id][position - 1]
            row[f"DNMT3L_NP_{position}_state"] = state if state not in "-?X" else "NA"
        rows.append(row)

    write_tsv(args.states_output, rows)

    summary = []
    for gene, positions in (("DNMT3A1", DNMT3A_SITES), ("DNMT3L_NP", DNMT3L_SITES)):
        for position in positions:
            field = f"{gene}_{position}_state"
            counts = Counter(str(row[field]) for row in rows)
            resolved = {key: value for key, value in counts.items() if key != "NA"}
            summary.append({
                "gene_coordinate": field,
                "matched_species": len(rows),
                "resolved_species": sum(resolved.values()),
                "distinct_resolved_states": len(resolved),
                "state_counts": ";".join(f"{state}:{count}" for state, count in sorted(resolved.items())),
                "missing_or_unresolved": counts.get("NA", 0),
                "primary_coevolution_informative": "True" if len(resolved) > 1 else "False",
                "interpretation": (
                    "variable in the matched panel"
                    if len(resolved) > 1
                    else "constant in the matched panel; cannot contribute a paired state transition"
                ),
            })
    write_tsv(args.summary_output, summary)


if __name__ == "__main__":
    main()
