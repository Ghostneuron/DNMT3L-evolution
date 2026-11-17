#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, NeighborSearch, ShrakeRupley


PROJECT = Path(__file__).resolve().parents[1]
LEGACY = PROJECT.parent / "Post_Genetica_Revision_20261025/GitHub_Release_Stage_v1.1.0/analysis_v2"
PDB = LEGACY / "structure/9MPP.cif"
CROSSWALK = LEGACY / "summaries/NP_037501_2_to_9MPP_Q9UJW3_crosswalk.csv"
INVENTORY = PROJECT / "results/therian_stem_full_substitution_inventory.tsv"
OUT = PROJECT / "results/therian_stem_interface_atomic_context.tsv"
SUMMARY = PROJECT / "results/therian_stem_interface_residue_context.tsv"
MANIFEST = PROJECT / "analysis/therian_interface_shift_337/STRUCTURE_CONTEXT_MANIFEST.json"
FOCAL = {261, 291, 293, 313, 349}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def partner_category(chain: str) -> str:
    if chain in {"K", "M"}:
        return "DNMT3A"
    if chain in set("ABCDEFGH"):
        return "histone"
    if chain in {"I", "J"}:
        return "DNA"
    return "other"


def main() -> None:
    inventory = pd.read_csv(INVENTORY, sep="\t")
    focal = inventory.loc[inventory["human_np_position"].isin(FOCAL)].copy()
    if set(focal["human_np_position"]) != FOCAL:
        raise ValueError("Focal interface positions are missing from the inventory")

    crosswalk = pd.read_csv(CROSSWALK)
    np_to_entity = {
        int(row.NP_037501_2_Position): int(row.Q9UJW3_9MPP_Entity_Position)
        for row in crosswalk.itertuples(index=False)
        if pd.notna(row.Q9UJW3_9MPP_Entity_Position)
    }
    model = MMCIFParser(QUIET=True).get_structure("9MPP", PDB)[0]
    ShrakeRupley().compute(model, level="R")
    partner_atoms = [
        atom
        for chain_id in tuple("ABCDEFGHIJKM")
        for atom in model[chain_id].get_atoms()
        if atom.element != "H"
    ]
    search = NeighborSearch(partner_atoms)

    atomic_rows: list[dict[str, object]] = []
    residue_rows: list[dict[str, object]] = []
    for item in focal.itertuples(index=False):
        position = int(item.human_np_position)
        entity_position = np_to_entity[position]
        for chain_id in ("L", "N"):
            residue = model[chain_id][(" ", entity_position, " ")]
            aggregates: dict[tuple[str, int, str, str], dict[str, object]] = defaultdict(
                lambda: {"minimum_heavy_atom_distance_A": np.inf, "atom_pair_count_within_5A": 0}
            )
            for focal_atom in residue.get_atoms():
                if focal_atom.element == "H":
                    continue
                for partner_atom in search.search(focal_atom.coord, 5.0, level="A"):
                    partner_residue = partner_atom.get_parent()
                    partner_chain = partner_residue.get_parent().id
                    key = (
                        partner_chain,
                        int(partner_residue.id[1]),
                        partner_residue.resname,
                        partner_category(partner_chain),
                    )
                    distance = float(np.linalg.norm(focal_atom.coord - partner_atom.coord))
                    aggregates[key]["minimum_heavy_atom_distance_A"] = min(
                        aggregates[key]["minimum_heavy_atom_distance_A"], distance
                    )
                    aggregates[key]["atom_pair_count_within_5A"] += 1
            categories = sorted({key[3] for key in aggregates})
            residue_rows.append({
                "human_np_position": position,
                "state_pair": item.state_pair,
                "evidence_class": (
                    "robust_therian_stem_change"
                    if bool(item.robust_therian_stem_change)
                    else "alignment_sensitive_change"
                ),
                "dnmt3l_chain": chain_id,
                "pdb_entity_position": entity_position,
                "pdb_residue_name": residue.resname,
                "residue_SASA_A2": float(residue.sasa),
                "contact_categories_within_5A": ";".join(categories),
                "contact_partner_residue_count": len(aggregates),
                "minimum_interface_distance_A": min(
                    (value["minimum_heavy_atom_distance_A"] for value in aggregates.values()),
                    default=np.nan,
                ),
            })
            for (partner_chain, partner_position, partner_name, category), values in sorted(aggregates.items()):
                atomic_rows.append({
                    "human_np_position": position,
                    "state_pair": item.state_pair,
                    "evidence_class": (
                        "robust_therian_stem_change"
                        if bool(item.robust_therian_stem_change)
                        else "alignment_sensitive_change"
                    ),
                    "dnmt3l_chain": chain_id,
                    "dnmt3l_entity_position": entity_position,
                    "dnmt3l_residue_name": residue.resname,
                    "partner_category": category,
                    "partner_chain": partner_chain,
                    "partner_residue_position": partner_position,
                    "partner_residue_name": partner_name,
                    "minimum_heavy_atom_distance_A": values["minimum_heavy_atom_distance_A"],
                    "atom_pair_count_within_5A": values["atom_pair_count_within_5A"],
                })

    atomic = pd.DataFrame(atomic_rows).sort_values(
        ["human_np_position", "dnmt3l_chain", "minimum_heavy_atom_distance_A"]
    )
    residue_summary = pd.DataFrame(residue_rows).sort_values(["human_np_position", "dnmt3l_chain"])
    atomic.to_csv(OUT, sep="\t", index=False)
    residue_summary.to_csv(SUMMARY, sep="\t", index=False)

    manifest = {
        "run_id": "THERIAN_INTERFACE_ATOMIC_CONTEXT_337_V1",
        "status": "VALIDATED",
        "command": "python scripts/map_therian_interface_changes_337.py",
        "interpretation_limit": "Observed human 9MPP contact localization only; no mutation energy or functional effect is estimated.",
        "input_hashes": {str(path.relative_to(PROJECT)) if path.is_relative_to(PROJECT) else str(path): sha256(path)
                         for path in (PDB, CROSSWALK, INVENTORY)},
        "output_hashes": {str(path.relative_to(PROJECT)): sha256(path) for path in (OUT, SUMMARY)},
        "validator": {
            "focal_positions": sorted(FOCAL),
            "resolved_dnmt3l_instances": len(residue_summary),
            "expected_dnmt3l_instances": 2 * len(FOCAL),
            "atomic_partner_rows": len(atomic),
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(residue_summary.to_string(index=False))


if __name__ == "__main__":
    main()
