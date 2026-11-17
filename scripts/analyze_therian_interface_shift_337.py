#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, NeighborSearch
from scipy.stats import fisher_exact


PROJECT = Path(__file__).resolve().parents[1]
ASR = PROJECT / "analysis/switching_helix_asr_337"
ANALYSIS = PROJECT / "analysis/therian_interface_shift_337"
RESULTS = PROJECT / "results"
LEGACY = PROJECT.parent / "Post_Genetica_Revision_20261025/GitHub_Release_Stage_v1.1.0/analysis_v2"

TAXONOMY_STATE = ASR / "runs/TAXONOMY_STRICT_LG_F_R4/DNMT3L_289_taxonomy_strict.state"
GENE_STATE = ASR / "runs/GENE_TREE_STRICT_LG_F_R4/DNMT3L_289_gene_tree_strict.state"
GLOBAL_ALIGNMENT = ASR / "inputs/DNMT3L_289_human_coordinate_MAFFT.fasta"
MATRIX = RESULTS / "DNMT3L_337_orthology_evidence_matrix.tsv"
MANIFEST = PROJECT / "metadata/functional_site_manifest.tsv"
PDB = LEGACY / "structure/9MPP.cif"
CROSSWALK = LEGACY / "summaries/NP_037501_2_to_9MPP_Q9UJW3_crosswalk.csv"
MACSE = {
    "Placentalia": PROJECT / "analysis/codon_alignments_337/Placentalia/MACSE/Placentalia_MACSE_AA.fasta",
    "Sauropsida": PROJECT / "analysis/codon_alignments_337/Sauropsida/MACSE/Sauropsida_MACSE_AA.fasta",
    "Amphibia": PROJECT / "analysis/codon_alignments_337/Amphibia/MACSE/Amphibia_MACSE_AA.fasta",
}

AA_CLASS = {
    "A": "nonpolar", "V": "nonpolar", "L": "nonpolar", "I": "nonpolar", "M": "nonpolar",
    "F": "aromatic", "W": "aromatic", "Y": "aromatic",
    "S": "polar", "T": "polar", "N": "polar", "Q": "polar", "C": "polar",
    "K": "positive", "R": "positive", "H": "positive",
    "D": "negative", "E": "negative", "G": "special", "P": "special",
}
AA_CHARGE = {aa: 0 for aa in AA_CLASS}
AA_CHARGE.update({"K": 1, "R": 1, "H": 0.1, "D": -1, "E": -1})
UNRESOLVED = {"-", "X", "?", "*", "!"}
PP_THRESHOLD = 0.80
SEED = 20261104
N_PERMUTATIONS = 100_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_iqtree_states(path: Path, nodes: dict[str, str]) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", comment="#")
    frame = frame.loc[frame["Node"].isin(nodes)].copy()
    frame["lineage"] = frame["Node"].map(nodes)
    frame["posterior"] = frame.apply(lambda row: float(row[f"p_{row['State']}"]), axis=1)
    return frame[["Node", "lineage", "Site", "State", "posterior"]]


def alignment_states(path: Path) -> tuple[dict[str, str], dict[int, int]]:
    records = {
        record.id: str(record.seq).upper().replace("!", "X")
        for record in SeqIO.parse(path, "fasta")
    }
    reference_name = "Homo_sapiens" if "Homo_sapiens" in records else "Homo_sapiens_anchor"
    if reference_name not in records:
        raise ValueError(f"No human coordinate anchor in {path}")
    reference = records[reference_name]
    columns: dict[int, int] = {}
    position = 0
    for column, state in enumerate(reference):
        if state in {"-", "*"}:
            continue
        position += 1
        columns[position] = column
    if position != 387:
        raise ValueError(f"{path} maps {position}, not 387 human positions")
    records.pop("Homo_sapiens_anchor", None)
    return records, columns


def state_at(sequence: str, columns: dict[int, int], position: int) -> str:
    return sequence[columns[position]]


def resolved(state: str) -> bool:
    return state not in UNRESOLVED


def alignment_audit() -> pd.DataFrame:
    global_records, global_columns = alignment_states(GLOBAL_ALIGNMENT)
    matrix = pd.read_csv(MATRIX, sep="\t")
    placental_subclades = {"Euarchontoglires", "Laurasiatheria", "Afrotheria", "Xenarthra"}
    clade_by_id = {
        row.safe_id: ("Placentalia" if row.clade in placental_subclades else row.clade)
        for row in matrix.itertuples(index=False)
    }
    selected = matrix.loc[matrix["Selection_Ready"].astype(str).eq("True")]
    expected = set(selected["safe_id"])
    if set(global_records) != expected:
        raise ValueError("Global ASR alignment and selected orthology matrix disagree")

    macse_data: dict[str, tuple[dict[str, str], dict[int, int]]] = {
        clade: alignment_states(path) for clade, path in MACSE.items()
    }
    rows: list[dict[str, object]] = []
    for position in range(1, 388):
        row: dict[str, object] = {"human_np_position": position}
        overall_paired = 0
        overall_agree = 0
        clade_gate = True
        failure_reasons: list[str] = []
        for clade in ["Amphibia", "Sauropsida", "Marsupialia", "Placentalia"]:
            names = [name for name in global_records if clade_by_id[name] == clade]
            global_calls = [state_at(global_records[name], global_columns, position) for name in names]
            global_resolved = sum(resolved(state) for state in global_calls)
            global_occupancy = global_resolved / len(names)
            row[f"{clade}_global_resolved"] = global_resolved
            row[f"{clade}_global_occupancy"] = global_occupancy
            clade_gate &= global_occupancy >= 0.80
            if global_occupancy < 0.80:
                failure_reasons.append(f"{clade}_global_occupancy_lt_0.80")

            if clade in macse_data:
                macse_records, macse_columns = macse_data[clade]
                biological = [name for name in names if name in macse_records]
                macse_calls = [state_at(macse_records[name], macse_columns, position) for name in biological]
                macse_resolved = sum(resolved(state) for state in macse_calls)
                macse_occupancy = macse_resolved / len(biological)
                paired = 0
                agree = 0
                for name in biological:
                    left = state_at(global_records[name], global_columns, position)
                    right = state_at(macse_records[name], macse_columns, position)
                    if resolved(left) and resolved(right):
                        paired += 1
                        agree += int(left == right)
                agreement = agree / paired if paired else np.nan
                row[f"{clade}_macse_resolved"] = macse_resolved
                row[f"{clade}_macse_occupancy"] = macse_occupancy
                row[f"{clade}_paired_calls"] = paired
                row[f"{clade}_agreement"] = agreement
                clade_gate &= macse_occupancy >= 0.80
                if macse_occupancy < 0.80:
                    failure_reasons.append(f"{clade}_macse_occupancy_lt_0.80")
                if paired >= 4:
                    clade_gate &= agreement >= 0.90
                    if agreement < 0.90:
                        failure_reasons.append(f"{clade}_agreement_lt_0.90")
                overall_paired += paired
                overall_agree += agree
            else:
                row[f"{clade}_macse_resolved"] = ""
                row[f"{clade}_macse_occupancy"] = ""
                row[f"{clade}_paired_calls"] = ""
                row[f"{clade}_agreement"] = ""
        overall_agreement = overall_agree / overall_paired if overall_paired else np.nan
        row["cross_method_paired_calls"] = overall_paired
        row["cross_method_agreement"] = overall_agreement
        if overall_agreement < 0.95:
            failure_reasons.append("overall_agreement_lt_0.95")
        row["alignment_gate_pass"] = bool(clade_gate and overall_agreement >= 0.95)
        row["alignment_gate_failure_reason"] = ";".join(failure_reasons)
        rows.append(row)
    return pd.DataFrame(rows)


def structure_features() -> pd.DataFrame:
    crosswalk = pd.read_csv(CROSSWALK)
    entity_to_np = {
        int(row.Q9UJW3_9MPP_Entity_Position): int(row.NP_037501_2_Position)
        for row in crosswalk.itertuples(index=False)
        if pd.notna(row.Q9UJW3_9MPP_Entity_Position)
    }
    model = MMCIFParser(QUIET=True).get_structure("9MPP", PDB)[0]
    dnmt3l_chains = [model[name] for name in ("L", "N")]
    dnmt3a_atoms = [atom for name in ("K", "M") for atom in model[name].get_atoms() if atom.element != "H"]
    nucleosome_atoms = [atom for name in tuple("ABCDEFGHIJ") for atom in model[name].get_atoms() if atom.element != "H"]
    dnmt3a_search = NeighborSearch(dnmt3a_atoms)
    nucleosome_search = NeighborSearch(nucleosome_atoms)

    features: dict[int, dict[str, object]] = {
        position: {
            "human_np_position": position,
            "structure_resolved": False,
            "dnmt3a_contact_5A": False,
            "nucleosome_contact_5A": False,
            "resolved_dnmt3l_copies": 0,
        }
        for position in range(1, 388)
    }
    for chain in dnmt3l_chains:
        seen: set[int] = set()
        for residue in chain:
            if residue.id[0] != " ":
                continue
            entity_position = int(residue.id[1])
            position = entity_to_np.get(entity_position)
            if position is None:
                continue
            features[position]["structure_resolved"] = True
            if position not in seen:
                features[position]["resolved_dnmt3l_copies"] += 1
                seen.add(position)
            atoms = [atom for atom in residue.get_atoms() if atom.element != "H"]
            if any(dnmt3a_search.search(atom.coord, 5.0, level="A") for atom in atoms):
                features[position]["dnmt3a_contact_5A"] = True
            if any(nucleosome_search.search(atom.coord, 5.0, level="A") for atom in atoms):
                features[position]["nucleosome_contact_5A"] = True
    frame = pd.DataFrame(features.values())
    frame["interface_union_5A"] = frame["dnmt3a_contact_5A"] | frame["nucleosome_contact_5A"]
    return frame


def bh_adjust(p_values: list[float]) -> list[float]:
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    adjusted = np.empty_like(values)
    running = 1.0
    n = len(values)
    for rank_index in range(n - 1, -1, -1):
        index = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[index] * n / rank)
        adjusted[index] = running
    return adjusted.tolist()


def fisher_row(name: str, universe: pd.DataFrame, mask: pd.Series, primary: bool) -> dict[str, object]:
    changed = universe["robust_therian_stem_change"].astype(bool)
    feature = mask.loc[universe.index].astype(bool)
    a = int((changed & feature).sum())
    b = int((~changed & feature).sum())
    c = int((changed & ~feature).sum())
    d = int((~changed & ~feature).sum())
    odds, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
    return {
        "test": name,
        "primary": primary,
        "universe_positions": len(universe),
        "feature_positions": int(feature.sum()),
        "changes_in_feature": a,
        "changes_outside_feature": c,
        "odds_ratio": odds,
        "fisher_one_sided_p": p_value,
    }


def main() -> None:
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    taxonomy = read_iqtree_states(TAXONOMY_STATE, {"NCBI_32524": "Amniota", "NCBI_32525": "Theria"})
    gene = read_iqtree_states(GENE_STATE, {"Node93": "Amniota", "Node92": "Theria"})
    alignment = alignment_audit()
    structure = structure_features()

    taxonomy_wide = taxonomy.pivot(index="Site", columns="lineage", values=["State", "posterior"])
    gene_wide = gene.pivot(index="Site", columns="lineage", values=["State", "posterior"])
    atlas = pd.read_csv(MANIFEST, sep="\t")
    atlas_positions = set(pd.to_numeric(atlas["human_np_position"], errors="coerce").dropna().astype(int))
    blosum = substitution_matrices.load("BLOSUM62")

    rows: list[dict[str, object]] = []
    for position in range(1, 388):
        t_a = str(taxonomy_wide.loc[position, ("State", "Amniota")])
        t_t = str(taxonomy_wide.loc[position, ("State", "Theria")])
        g_a = str(gene_wide.loc[position, ("State", "Amniota")])
        g_t = str(gene_wide.loc[position, ("State", "Theria")])
        t_a_pp = float(taxonomy_wide.loc[position, ("posterior", "Amniota")])
        t_t_pp = float(taxonomy_wide.loc[position, ("posterior", "Theria")])
        g_a_pp = float(gene_wide.loc[position, ("posterior", "Amniota")])
        g_t_pp = float(gene_wide.loc[position, ("posterior", "Theria")])
        node_pp_pass = min(t_a_pp, t_t_pp, g_a_pp, g_t_pp) >= PP_THRESHOLD
        topology_state_agreement = t_a == g_a and t_t == g_t
        reconstructed_change = t_a != t_t
        align_pass = bool(alignment.loc[position - 1, "alignment_gate_pass"])
        robust = node_pp_pass and topology_state_agreement and reconstructed_change and align_pass
        alignment_sensitive = node_pp_pass and topology_state_agreement and reconstructed_change and not align_pass
        eligible = node_pp_pass and topology_state_agreement and align_pass
        if 41 <= position <= 173:
            broad_region = "ADD"
        elif 178 <= position <= 380:
            broad_region = "C_terminal"
        else:
            broad_region = "other"
        old_state = t_a if topology_state_agreement else ""
        new_state = t_t if topology_state_agreement else ""
        rows.append({
            "human_np_position": position,
            "taxonomy_amniota_state": t_a,
            "taxonomy_amniota_pp": t_a_pp,
            "taxonomy_theria_state": t_t,
            "taxonomy_theria_pp": t_t_pp,
            "gene_tree_amniota_state": g_a,
            "gene_tree_amniota_pp": g_a_pp,
            "gene_tree_theria_state": g_t,
            "gene_tree_theria_pp": g_t_pp,
            "node_pp_pass": node_pp_pass,
            "topology_state_agreement": topology_state_agreement,
            "eligible_position": eligible,
            "reconstructed_change": reconstructed_change,
            "robust_therian_stem_change": robust,
            "alignment_sensitive_change": alignment_sensitive,
            "amniota_state": old_state,
            "theria_state": new_state,
            "state_pair": f"{old_state}>{new_state}" if old_state and new_state and old_state != new_state else "",
            "amniota_class": AA_CLASS.get(old_state, ""),
            "theria_class": AA_CLASS.get(new_state, ""),
            "class_change": bool(old_state and new_state and AA_CLASS.get(old_state) != AA_CLASS.get(new_state)),
            "charge_change": AA_CHARGE.get(new_state, np.nan) - AA_CHARGE.get(old_state, np.nan) if old_state and new_state else np.nan,
            "blosum62": float(blosum[old_state, new_state]) if old_state and new_state else np.nan,
            "broad_region": broad_region,
            "add_domain": 41 <= position <= 173,
            "c_terminal": 178 <= position <= 380,
            "switching_helix": 341 <= position <= 355,
            "adjacent_loop": 356 <= position <= 360,
            "functional_atlas_site": position in atlas_positions,
        })
    inventory = pd.DataFrame(rows).merge(alignment, on="human_np_position").merge(structure, on="human_np_position")

    eligible = inventory.loc[inventory["eligible_position"]].copy()
    structure_eligible = eligible.loc[eligible["structure_resolved"]].copy()
    cterm_eligible = eligible.loc[eligible["c_terminal"]].copy()
    tests = [
        fisher_row("interface_union_5A", structure_eligible, inventory["interface_union_5A"], True),
        fisher_row("switching_helix_vs_other_C_terminal", cterm_eligible, inventory["switching_helix"], True),
        fisher_row("ADD_domain", eligible, inventory["add_domain"], False),
        fisher_row("C_terminal_region", eligible, inventory["c_terminal"], False),
        fisher_row("adjacent_loop", eligible, inventory["adjacent_loop"], False),
        fisher_row("DNMT3A_contact_5A", structure_eligible, inventory["dnmt3a_contact_5A"], False),
        fisher_row("nucleosome_contact_5A", structure_eligible, inventory["nucleosome_contact_5A"], False),
        fisher_row("functional_atlas_sites", eligible, inventory["functional_atlas_site"], False),
    ]
    enrichment = pd.DataFrame(tests)
    enrichment["BH_q"] = bh_adjust(enrichment["fisher_one_sided_p"].tolist())

    robust = inventory.loc[inventory["robust_therian_stem_change"]].copy()
    rng = np.random.default_rng(SEED)
    observed_interface = int(robust["interface_union_5A"].sum())
    strata = defaultdict(lambda: {"positions": [], "n_changes": 0})
    for row in eligible.itertuples(index=False):
        key = (row.broad_region, bool(row.structure_resolved))
        strata[key]["positions"].append(int(row.human_np_position))
        strata[key]["n_changes"] += int(row.robust_therian_stem_change)
    interface_lookup = dict(zip(inventory["human_np_position"], inventory["interface_union_5A"]))
    simulated = np.zeros(N_PERMUTATIONS, dtype=np.int16)
    for draw in range(N_PERMUTATIONS):
        total = 0
        for values in strata.values():
            count = values["n_changes"]
            if count:
                chosen = rng.choice(values["positions"], size=count, replace=False)
                total += sum(bool(interface_lookup[int(position)]) for position in chosen)
        simulated[draw] = total
    empirical_p = (1 + int((simulated >= observed_interface).sum())) / (N_PERMUTATIONS + 1)

    primary_interface = enrichment.loc[enrichment["test"].eq("interface_union_5A")].iloc[0]
    if observed_interface >= 2 and primary_interface.BH_q < 0.05 and empirical_p < 0.05:
        decision = "GO_COORDINATED_INTERFACE_REMODELING"
    elif observed_interface >= 1:
        decision = "NARROW_INTERFACE_HYPOTHESIS"
    elif len(robust):
        decision = "DESCRIPTIVE_BRANCH_CHANGE"
    else:
        decision = "STOP_FUNCTIONAL_SHIFT_CLAIM"

    inventory_path = RESULTS / "therian_stem_full_substitution_inventory.tsv"
    robust_path = RESULTS / "therian_stem_robust_substitutions.tsv"
    alignment_path = RESULTS / "therian_stem_alignment_robustness.tsv"
    structure_path = RESULTS / "therian_stem_structure_contacts.tsv"
    enrichment_path = RESULTS / "therian_stem_interface_enrichment.tsv"
    summary_path = RESULTS / "therian_stem_interface_shift_summary.tsv"
    inventory.to_csv(inventory_path, sep="\t", index=False)
    robust.to_csv(robust_path, sep="\t", index=False)
    alignment.to_csv(alignment_path, sep="\t", index=False)
    structure.to_csv(structure_path, sep="\t", index=False)
    enrichment.to_csv(enrichment_path, sep="\t", index=False)
    summary = pd.DataFrame([
        {
            "decision": decision,
            "eligible_positions": len(eligible),
            "robust_therian_stem_changes": len(robust),
            "alignment_sensitive_changes": int(inventory["alignment_sensitive_change"].sum()),
            "structure_resolved_eligible_positions": len(structure_eligible),
            "robust_interface_union_changes": observed_interface,
            "primary_interface_fisher_p": primary_interface.fisher_one_sided_p,
            "primary_interface_BH_q": primary_interface.BH_q,
            "stratified_permutation_p": empirical_p,
            "permutation_draws": N_PERMUTATIONS,
            "seed": SEED,
        }
    ])
    summary.to_csv(summary_path, sep="\t", index=False)

    inputs = [TAXONOMY_STATE, GENE_STATE, GLOBAL_ALIGNMENT, MATRIX, MANIFEST, PDB, CROSSWALK, *MACSE.values()]
    run_record = {
        "run_id": "THERIAN_INTERFACE_SHIFT_337_V1",
        "status": "VALIDATED",
        "command": "python scripts/analyze_therian_interface_shift_337.py",
        "decision_rules": "analysis/therian_interface_shift_337/DECISION_RULES.md",
        "seed": SEED,
        "permutation_draws": N_PERMUTATIONS,
        "input_hashes": {str(path): sha256(path) for path in inputs},
        "output_hashes": {str(path.relative_to(PROJECT)): sha256(path) for path in [inventory_path, robust_path, alignment_path, structure_path, enrichment_path, summary_path]},
        "validator": {
            "inventory_rows": len(inventory),
            "expected_inventory_rows": 387,
            "eligible_positions": len(eligible),
            "robust_changes": len(robust),
            "decision": decision,
        },
    }
    (ANALYSIS / "RUN_MANIFEST.json").write_text(json.dumps(run_record, indent=2) + "\n")

    print(json.dumps(summary.iloc[0].to_dict(), indent=2))
    print("Robust changes:")
    for row in robust.itertuples(index=False):
        print(
            f"NP{row.human_np_position} {row.state_pair} region={row.broad_region} "
            f"helix={row.switching_helix} DNMT3A_contact={row.dnmt3a_contact_5A} "
            f"nucleosome_contact={row.nucleosome_contact_5A}"
        )


if __name__ == "__main__":
    main()
