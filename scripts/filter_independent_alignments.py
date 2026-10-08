from __future__ import annotations

import csv
import hashlib
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

from Bio.Seq import Seq


PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / "analysis_v2" / "raw_alignments"
FILTERED = PROJECT / "analysis_v2" / "filtered_alignments"
SELECTION = PROJECT / "analysis_v2" / "selection_inputs"
SUMMARIES = PROJECT / "analysis_v2" / "summaries"

SCOPES = ("Placentalia", "Euarchontoglires", "Laurasiatheria", "Sauropsida")
ALIGNERS = ("MACSE", "MAFFT", "PRANK")
MIN_OCCUPANCY = 0.80


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


def codon_state(codon: str) -> tuple[str, str]:
    codon = codon.upper().replace("!", "-")
    if codon == "---":
        return "gap", "-"
    if len(codon) == 3 and set(codon) <= set("ACGT"):
        amino_acid = str(Seq(codon).translate())
        return ("stop", "*") if amino_acid == "*" else ("resolved", amino_acid)
    return "ambiguous", "X"


def normalized_codon(codon: str) -> str:
    codon = codon.upper().replace("!", "-")
    state, _ = codon_state(codon)
    if state == "resolved":
        return codon
    if state == "gap":
        return "---"
    return "NNN"


def map_reference_positions(sequence: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    reference_position = 0
    for site, start in enumerate(range(0, len(sequence), 3), start=1):
        state, _ = codon_state(sequence[start : start + 3])
        if state == "resolved":
            reference_position += 1
            mapping[reference_position] = site
    return mapping


def domain_for(position: int) -> str:
    if position <= 40:
        return "N_terminal"
    if position <= 173:
        return "ADD"
    if position <= 177:
        return "linker"
    if position <= 380:
        return "C_terminal_MTase_like"
    return "C_terminal_tail"


def structural_region_for(position: int) -> str:
    if 341 <= position <= 355:
        return "switching_helix"
    if 356 <= position <= 360:
        return "switching_helix_adjacent_loop"
    return "other"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    FILTERED.mkdir(parents=True, exist_ok=True)
    SELECTION.mkdir(parents=True, exist_ok=True)
    SUMMARIES.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, object]] = []
    site_map_rows: list[dict[str, object]] = []
    duplicate_rows: list[dict[str, object]] = []
    residue_maps: dict[tuple[str, str], dict[tuple[str, int], tuple[str, str]]] = {}

    for scope in SCOPES:
        anchor = "Homo_sapiens_coordinate_anchor" if scope in {"Laurasiatheria", "Sauropsida"} else "Homo_sapiens"
        for aligner in ALIGNERS:
            label = f"{scope}_{aligner}"
            raw_path = RAW / scope / aligner / f"{label}_NT.fasta"
            records = read_fasta(raw_path)
            names = [name for name, _ in records]
            if len(names) != len(set(names)):
                raise ValueError(f"Duplicate FASTA IDs in {raw_path}")
            lengths = {len(sequence) for _, sequence in records}
            if len(lengths) != 1:
                raise ValueError(f"Unequal sequence lengths in {raw_path}")
            alignment_nt = lengths.pop()
            if alignment_nt % 3:
                raise ValueError(f"Non-codon alignment length in {raw_path}")
            by_name = dict(records)
            if anchor not in by_name:
                raise ValueError(f"Missing anchor {anchor} in {raw_path}")
            reference_map = map_reference_positions(by_name[anchor])
            biological_records = [(name, sequence) for name, sequence in records if name != anchor or scope not in {"Laurasiatheria", "Sauropsida"}]

            retained_positions: list[int] = []
            codons_by_species: dict[str, list[str]] = {name: [] for name, _ in biological_records}
            residue_map: dict[tuple[str, int], tuple[str, str]] = {}
            raw_site_stats: dict[int, dict[str, object]] = {}

            for human_position, raw_site in sorted(reference_map.items()):
                resolved = 0
                residue_counts: Counter[str] = Counter()
                for name, sequence in biological_records:
                    codon = sequence[3 * (raw_site - 1) : 3 * raw_site]
                    state, residue = codon_state(codon)
                    residue_map[(name, human_position)] = (state, residue)
                    if state == "resolved":
                        resolved += 1
                        residue_counts[residue] += 1
                occupancy = resolved / len(biological_records)
                entropy = 0.0
                if resolved:
                    for count in residue_counts.values():
                        frequency = count / resolved
                        entropy -= frequency * math.log2(frequency)
                raw_site_stats[human_position] = {
                    "resolved": resolved,
                    "occupancy": occupancy,
                    "residue_counts": residue_counts,
                    "entropy": entropy,
                }
                if occupancy >= MIN_OCCUPANCY:
                    retained_positions.append(human_position)
                    for name, sequence in biological_records:
                        codon = sequence[3 * (raw_site - 1) : 3 * raw_site]
                        codons_by_species[name].append(normalized_codon(codon))

            output_records = [(name, "".join(codons_by_species[name])) for name, _ in biological_records]
            filtered_path = FILTERED / f"{label}_occupancy80.fasta"
            write_fasta(output_records, filtered_path)

            for filtered_site, human_position in enumerate(retained_positions, start=1):
                stats = raw_site_stats[human_position]
                resolved = int(stats["resolved"])
                occupancy = float(stats["occupancy"])
                residue_counts = stats["residue_counts"]
                site_map_rows.append(
                    {
                        "Scope": scope,
                        "Aligner": aligner,
                        "Filtered_Site": filtered_site,
                        "Human_Position": human_position,
                        "Domain": domain_for(human_position),
                        "Structural_Region": structural_region_for(human_position),
                        "Resolved_N": resolved,
                        "Occupancy": occupancy,
                        "Amino_Acid_State_Count": len(residue_counts),
                        "Major_Amino_Acid_Frequency": max(residue_counts.values()) / resolved if resolved else "",
                        "Amino_Acid_Shannon_Entropy": stats["entropy"],
                        "Amino_Acid_Counts": ";".join(f"{residue}:{count}" for residue, count in sorted(residue_counts.items())),
                    }
                )

            names_by_sequence: dict[str, list[str]] = defaultdict(list)
            for name, sequence in output_records:
                names_by_sequence[sequence].append(name)
            retained_names: set[str] = set()
            for sequence, sequence_names in names_by_sequence.items():
                representative = "Homo_sapiens" if "Homo_sapiens" in sequence_names else sorted(sequence_names)[0]
                retained_names.add(representative)
                sequence_hash = hashlib.sha256(sequence.encode()).hexdigest()
                for name in sorted(sequence_names):
                    duplicate_rows.append(
                        {
                            "Scope": scope,
                            "Aligner": aligner,
                            "Sequence_SHA256": sequence_hash,
                            "Representative": representative,
                            "Sequence_Name": name,
                            "Status": "retained" if name == representative else "removed_exact_duplicate",
                            "Duplicate_Group_Size": len(sequence_names),
                        }
                    )
            selection_records = [(name, sequence) for name, sequence in output_records if name in retained_names]
            selection_path = SELECTION / f"{label}_selection_input.fasta"
            write_fasta(selection_records, selection_path)

            missing_fractions = []
            for _, sequence in output_records:
                states = [codon_state(sequence[start : start + 3])[0] for start in range(0, len(sequence), 3)]
                missing_fractions.append(sum(state != "resolved" for state in states) / len(states))
            focal_site = retained_positions.index(358) + 1 if 358 in retained_positions else ""
            summary_rows.append(
                {
                    "Scope": scope,
                    "Aligner": aligner,
                    "Raw_Sequence_Count_Including_Anchor": len(records),
                    "Biological_Sequence_Count": len(output_records),
                    "Raw_Alignment_Codons": alignment_nt // 3,
                    "Human_Positions_Mapped": len(reference_map),
                    "Retained_Codons": len(retained_positions),
                    "First_Human_Position": min(retained_positions),
                    "Last_Human_Position": max(retained_positions),
                    "h358_Filtered_Site": focal_site,
                    "h358_Raw_Occupancy": float(raw_site_stats.get(358, {"occupancy": 0.0})["occupancy"]),
                    "Median_Sequence_Missing_Pct": 100 * median(missing_fractions),
                    "Exact_Duplicates_Removed": len(output_records) - len(selection_records),
                    "Selection_Sequence_Count": len(selection_records),
                    "Filtered_SHA256": sha256(filtered_path),
                    "Selection_Input_SHA256": sha256(selection_path),
                }
            )
            residue_maps[(scope, aligner)] = residue_map
            print(f"{label}: {len(retained_positions)} codons, {len(selection_records)} unique sequences")

    concordance_rows: list[dict[str, object]] = []
    for scope in SCOPES:
        maps = {aligner: residue_maps[(scope, aligner)] for aligner in ALIGNERS}
        species = sorted(set.intersection(*(set(name for name, _ in data) for data in maps.values())))
        positions = sorted(set.intersection(*(set(position for _, position in data) for data in maps.values())))
        for position in positions:
            all_resolved = 0
            all_exact = 0
            any_resolved = 0
            pairwise_exact = Counter({"MACSE_MAFFT": 0, "MACSE_PRANK": 0, "MAFFT_PRANK": 0})
            pairwise_resolved = Counter({"MACSE_MAFFT": 0, "MACSE_PRANK": 0, "MAFFT_PRANK": 0})
            for name in species:
                calls = {aligner: maps[aligner][(name, position)] for aligner in ALIGNERS}
                resolved_calls = {aligner: residue for aligner, (state, residue) in calls.items() if state == "resolved"}
                if resolved_calls:
                    any_resolved += 1
                if len(resolved_calls) == 3:
                    all_resolved += 1
                    if len(set(resolved_calls.values())) == 1:
                        all_exact += 1
                for first, second in (("MACSE", "MAFFT"), ("MACSE", "PRANK"), ("MAFFT", "PRANK")):
                    key = f"{first}_{second}"
                    if first in resolved_calls and second in resolved_calls:
                        pairwise_resolved[key] += 1
                        if resolved_calls[first] == resolved_calls[second]:
                            pairwise_exact[key] += 1
            concordance_rows.append(
                {
                    "Scope": scope,
                    "Human_Position": position,
                    "Domain": domain_for(position),
                    "Structural_Region": structural_region_for(position),
                    "Species_N": len(species),
                    "Any_Aligner_Resolved_N": any_resolved,
                    "All_Three_Resolved_N": all_resolved,
                    "All_Three_Resolved_Fraction": all_resolved / len(species),
                    "All_Three_Exact_AA_N": all_exact,
                    "All_Three_Exact_Among_Resolved": all_exact / all_resolved if all_resolved else "",
                    "All_Three_Exact_Among_All_Taxa": all_exact / len(species),
                    "Alignment_Stable_Coverage80_Concordance90": (
                        all_resolved / len(species) >= 0.80
                        and all_resolved > 0
                        and all_exact / all_resolved >= 0.90
                    ),
                    "MACSE_MAFFT_Exact_Among_Resolved": pairwise_exact["MACSE_MAFFT"] / pairwise_resolved["MACSE_MAFFT"] if pairwise_resolved["MACSE_MAFFT"] else "",
                    "MACSE_PRANK_Exact_Among_Resolved": pairwise_exact["MACSE_PRANK"] / pairwise_resolved["MACSE_PRANK"] if pairwise_resolved["MACSE_PRANK"] else "",
                    "MAFFT_PRANK_Exact_Among_Resolved": pairwise_exact["MAFFT_PRANK"] / pairwise_resolved["MAFFT_PRANK"] if pairwise_resolved["MAFFT_PRANK"] else "",
                }
            )

    write_csv(SUMMARIES / "filtered_alignment_summary.csv", summary_rows)
    write_csv(SUMMARIES / "retained_site_map.csv", site_map_rows)
    write_csv(SUMMARIES / "exact_duplicate_manifest.csv", duplicate_rows)
    write_csv(SUMMARIES / "cross_aligner_concordance_by_site.csv", concordance_rows)


if __name__ == "__main__":
    main()
