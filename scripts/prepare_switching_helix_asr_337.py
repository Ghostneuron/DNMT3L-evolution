#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from Bio import Phylo, SeqIO


PROJECT = Path(__file__).resolve().parents[1]
MATRIX = PROJECT / "results/DNMT3L_337_orthology_evidence_matrix.tsv"
DATA_REPORT = (
    PROJECT
    / "data/raw/ncbi_dataset_DNMT3L_orthologs_20261102/ncbi_dataset/data/data_report.jsonl"
)
FAMILY_ALIGNMENT = (
    PROJECT
    / "analysis/gene_family_tree_337/DNMT3_family_337_MAFFT_human_DNMT3L_columns.fasta"
)
FAMILY_TREE = (
    PROJECT / "analysis/gene_family_tree_337/DNMT3_family_337_LG_F_R4_fast.treefile"
)
ATLAS = PROJECT / "results/functional_site_states_by_record.tsv"
MACSE_ALIGNMENTS = {
    "Placentalia": PROJECT / "analysis/codon_alignments_337/Placentalia/MACSE/Placentalia_MACSE_AA.fasta",
    "Sauropsida": PROJECT / "analysis/codon_alignments_337/Sauropsida/MACSE/Sauropsida_MACSE_AA.fasta",
    "Amphibia": PROJECT / "analysis/codon_alignments_337/Amphibia/MACSE/Amphibia_MACSE_AA.fasta",
}
ROOT = PROJECT / "analysis/switching_helix_asr_337"
INPUTS = ROOT / "inputs"
TREES = ROOT / "trees"
TAXONOMY = ROOT / "taxonomy"
RESULTS = PROJECT / "results"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_fasta(records: dict[str, str], path: Path) -> None:
    with path.open("w") as handle:
        for name in sorted(records):
            handle.write(f">{name}\n")
            sequence = records[name]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def reference_coordinate_states(path: Path) -> dict[tuple[str, int], str]:
    records = {record.id: str(record.seq).upper().replace("!", "X") for record in SeqIO.parse(path, "fasta")}
    reference_name = "Homo_sapiens" if "Homo_sapiens" in records else "Homo_sapiens_anchor"
    reference = records[reference_name]
    columns: dict[int, int] = {}
    position = 0
    for column, state in enumerate(reference):
        if state in {"-", "*"}:
            continue
        position += 1
        columns[position] = column
    if position != 387:
        raise ValueError(f"{path} maps {position}, not 387, human positions")
    output: dict[tuple[str, int], str] = {}
    for name, sequence in records.items():
        if name == "Homo_sapiens_anchor":
            continue
        for position in range(341, 356):
            output[(name, position)] = sequence[columns[position]]
    return output


def main() -> None:
    for directory in (INPUTS, TREES, TAXONOMY):
        directory.mkdir(parents=True, exist_ok=True)

    with MATRIX.open(newline="") as handle:
        selected_rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["Selection_Ready"] == "True"
        ]
    if len(selected_rows) != 289:
        raise ValueError(f"Expected 289 selection-ready records, observed {len(selected_rows)}")
    selected = {row["safe_id"]: row for row in selected_rows}

    tax_by_gene: dict[str, tuple[str, str]] = {}
    with DATA_REPORT.open() as handle:
        for line in handle:
            record = json.loads(line)
            tax_by_gene[str(record["geneId"])] = (str(record["taxId"]), record["taxname"])

    family_records = {
        record.id.removeprefix("DNMT3L__"): str(record.seq).upper()
        for record in SeqIO.parse(FAMILY_ALIGNMENT, "fasta")
        if record.id.startswith("DNMT3L__")
    }
    if set(family_records) != set(selected):
        raise ValueError(
            "Family alignment mismatch: "
            f"missing={sorted(set(selected)-set(family_records))}, "
            f"extra={sorted(set(family_records)-set(selected))}"
        )
    lengths = {len(sequence) for sequence in family_records.values()}
    if lengths != {387}:
        raise ValueError(f"Expected 387 human-coordinate columns, observed {sorted(lengths)}")
    human = family_records["Homo_sapiens"]
    if "-" in human or human != selected["Homo_sapiens"]["protein_seq"]:
        raise ValueError("Human-coordinate projection does not reproduce NP_037501.2")

    alignment_path = INPUTS / "DNMT3L_289_human_coordinate_MAFFT.fasta"
    write_fasta(family_records, alignment_path)

    macse_states: dict[tuple[str, int], str] = {}
    for path in MACSE_ALIGNMENTS.values():
        for key, state in reference_coordinate_states(path).items():
            if key in macse_states and macse_states[key] != state:
                raise ValueError(f"Nested MACSE alignments disagree for {key}: {macse_states[key]} versus {state}")
            macse_states[key] = state

    map_rows: list[dict[str, str]] = []
    for safe_id in sorted(selected):
        row = selected[safe_id]
        taxid, taxname = tax_by_gene[row["gene_id"]]
        map_rows.append(
            {
                "safe_id": safe_id,
                "species": row["species"],
                "gene_id": row["gene_id"],
                "taxid": taxid,
                "taxname": taxname,
                "clade": row["clade"],
            }
        )
    if len({row["taxid"] for row in map_rows}) != len(map_rows):
        repeats: dict[str, list[str]] = {}
        for row in map_rows:
            repeats.setdefault(row["taxid"], []).append(row["safe_id"])
        repeats = {key: value for key, value in repeats.items() if len(value) > 1}
        raise ValueError(f"Repeated taxonomy IDs: {repeats}")

    map_path = TAXONOMY / "DNMT3L_289_taxid_map.tsv"
    with map_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(map_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(map_rows)
    ids_path = TAXONOMY / "DNMT3L_289_taxids.txt"
    ids_path.write_text(",".join(row["taxid"] for row in map_rows) + "\n")

    tree = Phylo.read(FAMILY_TREE, "newick")
    keep_names = {f"DNMT3L__{name}" for name in selected}
    for tip in list(tree.get_terminals()):
        if tip.name not in keep_names:
            tree.prune(tip)
    for tip in tree.get_terminals():
        tip.name = tip.name.removeprefix("DNMT3L__")
    observed_tips = {tip.name for tip in tree.get_terminals()}
    if observed_tips != set(selected):
        raise ValueError("Pruned DNMT3L gene tree has an incorrect tip set")
    gene_tree_path = TREES / "DNMT3L_289_gene_tree_topology.nwk"
    Phylo.write(tree, gene_tree_path, "newick")

    atlas_lookup: dict[tuple[str, int], tuple[str, str]] = {}
    with ATLAS.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["feature_id"] in {"SWITCH_Q349", "SWITCH_K351", "SWITCH_Q352", "SWITCH_K355"}:
                atlas_lookup[(row["safe_id"], int(row["human_np_position"]))] = (
                    row["query_residue"],
                    row["mapping_status"],
                )

    site_rows: list[dict[str, object]] = []
    mismatch_count = 0
    strict_records = dict(family_records)
    functional_positions = {349, 351, 352, 355}
    for safe_id in sorted(family_records):
        sequence = family_records[safe_id]
        for position in range(341, 356):
            state = sequence[position - 1]
            atlas_state, atlas_status = atlas_lookup.get((safe_id, position), ("", "NOT_ATLAS_SITE"))
            macse_state = macse_states.get((safe_id, position), "")
            agrees = ""
            if atlas_state:
                agrees = state == atlas_state
                mismatch_count += int(not agrees)
            resolved_states = [value for value in (state, atlas_state, macse_state) if value and value not in {"-", "X", "?"}]
            strict_concordant = bool(resolved_states) and len(set(resolved_states)) == 1
            if position in functional_positions and not strict_concordant:
                edited = list(strict_records[safe_id])
                edited[position - 1] = "X"
                strict_records[safe_id] = "".join(edited)
            site_rows.append(
                {
                    "safe_id": safe_id,
                    "species": selected[safe_id]["species"],
                    "clade": selected[safe_id]["clade"],
                    "human_np_position": position,
                    "mafft_human_coordinate_state": state,
                    "resolved": state not in {"-", "X", "?"},
                    "macse_clade_state": macse_state,
                    "atlas_query_state": atlas_state,
                    "atlas_mapping_status": atlas_status,
                    "mafft_atlas_agreement": agrees,
                    "strict_cross_method_concordance": strict_concordant,
                }
            )
    strict_alignment_path = INPUTS / "DNMT3L_289_human_coordinate_strict_focal_masked.fasta"
    write_fasta(strict_records, strict_alignment_path)
    site_path = RESULTS / "switching_helix_337_observed_states.tsv"
    with site_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(site_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(site_rows)

    manifest_rows = [
        {
            "artifact": str(path.relative_to(PROJECT)),
            "sha256": sha256(path),
            "role": role,
        }
        for path, role in (
            (MATRIX, "orthology/QC gate"),
            (DATA_REPORT, "NCBI taxonomy identifiers"),
            (FAMILY_ALIGNMENT, "source MAFFT alignment projected to human coordinates"),
            (FAMILY_TREE, "source DNMT3-family gene tree"),
            (alignment_path, "289-tip ASR protein alignment"),
            (strict_alignment_path, "289-tip ASR protein alignment with discordant focal calls masked"),
            (map_path, "289-tip taxonomy map"),
            (ids_path, "NCBI taxonomy query"),
            (gene_tree_path, "DNMT3L gene-tree topology sensitivity"),
            (site_path, "observed switching-helix state audit"),
        )
    ]
    manifest_path = ROOT / "INPUT_MANIFEST.tsv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"Prepared {len(selected)} ASR tips and {len(site_rows)} switching-helix calls")
    print(f"Atlas/MAFFT mismatches at the four atlas sites: {mismatch_count}")
    for position in sorted(functional_positions):
        masked = sum(strict_records[name][position - 1] == "X" for name in strict_records)
        print(f"Strict focal masks at NP{position}: {masked}")
    print(f"Alignment: {alignment_path.relative_to(PROJECT)}")
    print(f"Taxonomy query: {ids_path.relative_to(PROJECT)}")
    print(f"Gene tree: {gene_tree_path.relative_to(PROJECT)}")


if __name__ == "__main__":
    main()
