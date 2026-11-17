#!/usr/bin/env python3
"""Confidence-gated structural analysis of the completed DNMT3L sentinel panel."""

from __future__ import annotations

import csv
import itertools
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.SVDSuperimposer import SVDSuperimposer


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/sentinel_analysis"
DOMAINS = {"ADD": (49, 181), "CTD": (178, 380)}
FOCAL = [261, 291, 293, 313, 349]
ANCESTORS = {"ANCESTRAL_AMNIOTA", "ANCESTRAL_THERIA"}


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def one(path, pattern):
    found = [item for item in path.glob(pattern) if not item.name.startswith("._")]
    if len(found) != 1:
        raise ValueError((path, pattern, len(found)))
    return found[0]


def superpose(reference, moving):
    sup = SVDSuperimposer()
    sup.set(np.asarray(reference, float), np.asarray(moving, float))
    sup.run()
    rotation, translation = sup.get_rotran()
    transformed = np.asarray(moving, float) @ rotation + translation
    return float(sup.get_rms()), transformed, rotation, translation


def jaccard(left, right):
    union = left | right
    return 1.0 if not union else len(left & right) / len(union)


def residue_atoms(residue):
    return [atom for atom in residue if atom.element != "H"]


def contact_pairs(chain_a, chain_b, cutoff=5.0):
    atoms_b = [atom for residue in chain_b for atom in residue_atoms(residue)]
    search = NeighborSearch(atoms_b)
    b_index = {id(residue): index for index, residue in enumerate(chain_b)}
    pairs = set()
    for a_index, residue in enumerate(chain_a):
        for atom in residue_atoms(residue):
            for hit in search.search(atom.coord, cutoff, "A"):
                pairs.add((a_index, b_index[id(hit.get_parent())]))
    return pairs


def min_chain_distance(residue, other_chain):
    other = np.asarray([atom.coord for other_residue in other_chain for atom in residue_atoms(other_residue)])
    current = np.asarray([atom.coord for atom in residue_atoms(residue)])
    if not len(current) or not len(other):
        return np.nan
    return float(np.sqrt(((current[:, None, :] - other[None, :, :]) ** 2).sum(axis=2)).min())


def main():
    status = read_tsv(ROOT / "metadata/sentinel_run_status.tsv")
    if len(status) != 42 or any(row["status"] != "VALIDATED" for row in status):
        raise ValueError("Sentinel completion gate has not passed")
    manifest = {row["safe_id"]: row for row in read_tsv(ROOT / "metadata/panel_manifest.tsv")}
    crosswalk = defaultdict(dict)
    reverse_b = defaultdict(dict)
    for row in read_tsv(ROOT / "metadata/human_coordinate_crosswalk.tsv"):
        target = row["safe_id"]
        human = int(row["human_position"])
        crosswalk[target][human] = {
            "monomer": int(row["monomer_input_position"]),
            "b": int(row["complex_chain_B_position"]) if row["complex_chain_B_position"] else None,
            "pp": float(row["ancestral_posterior_probability"]) if row["ancestral_posterior_probability"] else None,
            "residue": row["monomer_residue"],
        }
        if row["complex_chain_B_position"]:
            reverse_b[target][int(row["complex_chain_B_position"]) - 1] = human

    predictions = {}
    metric_rows = []
    focal_rows = []
    for row in status:
        mode, target, seed = row["mode"], row["target"], int(row["seed"])
        root = ROOT / row["output"]
        confidence = json.loads(one(root, "confidence_*.json").read_text())
        plddt = np.load(one(root, "plddt_*.npz"))["plddt"]
        pae = np.load(one(root, "pae_*.npz"))["pae"]
        structure = MMCIFParser(QUIET=True).get_structure("prediction", one(root, "*.cif"))[0]
        chains = {chain.id: list(chain.get_residues()) for chain in structure}
        predictions[(mode, target, seed)] = {"confidence": confidence, "plddt": plddt, "pae": pae, "chains": chains}

        base = {
            "mode": mode, "target": target, "clade": manifest[target]["lineage_group"], "seed": seed,
            "confidence_score": f"{confidence['confidence_score']:.6f}",
            "mean_plddt": f"{float(plddt.mean()):.6f}",
            "structure_checks_passed": str(confidence["structure_checks_passed"]),
        }
        if mode == "monomer":
            for domain, (start, end) in DOMAINS.items():
                mapped = [crosswalk[target][h]["monomer"] - 1 for h in range(start, end + 1) if h in crosswalk[target]]
                values = plddt[mapped]
                base[f"{domain}_mapped_n"] = len(mapped)
                base[f"{domain}_plddt_ge_0.70_fraction"] = f"{float((values >= .70).mean()):.6f}"
                base[f"{domain}_confidence_pass"] = str(float((values >= .70).mean()) >= .80)
            base.update({"iptm": "", "interchain_pae": "", "interface_residue_n": "", "interface_plddt_median": "", "interface_pae_median": "", "interface_confidence_pass": ""})
        else:
            chain_a, chain_b = chains["A"], chains["B"]
            pairs = contact_pairs(chain_a, chain_b)
            a_interface = sorted({a for a, _ in pairs})
            b_interface = sorted({b for _, b in pairs})
            interface_values = np.r_[plddt[a_interface], plddt[290 + np.asarray(b_interface, int)]]
            contact_pae = np.asarray([value for a, b in pairs for value in (pae[a, 290 + b], pae[290 + b, a])])
            interface_pass = confidence["iptm"] >= .70 and np.median(interface_values) >= .70 and np.median(contact_pae) <= 10
            base.update({
                "ADD_mapped_n": "", "ADD_plddt_ge_0.70_fraction": "", "ADD_confidence_pass": "",
                "CTD_mapped_n": "", "CTD_plddt_ge_0.70_fraction": "", "CTD_confidence_pass": "",
                "iptm": f"{confidence['iptm']:.6f}",
                "interchain_pae": f"{confidence['complex_ipae']:.6f}",
                "interface_residue_n": len(a_interface) + len(b_interface),
                "interface_plddt_median": f"{float(np.median(interface_values)):.6f}",
                "interface_pae_median": f"{float(np.median(contact_pae)):.6f}",
                "interface_confidence_pass": str(interface_pass),
            })
            for human in FOCAL:
                mapping = crosswalk[target].get(human)
                if not mapping or mapping["b"] is None:
                    continue
                b_index = mapping["b"] - 1
                local_pairs = {(a, b) for a, b in pairs if b == b_index}
                focal_rows.append({
                    "target": target, "clade": manifest[target]["lineage_group"], "seed": seed,
                    "human_position": human, "residue": mapping["residue"],
                    "ancestral_posterior_probability": "" if mapping["pp"] is None else f"{mapping['pp']:.5f}",
                    "plddt": f"{plddt[290 + b_index]:.6f}",
                    "minimum_DNMT3A_heavy_atom_distance_A": f"{min_chain_distance(chain_b[b_index], chain_a):.4f}",
                    "DNMT3A_contact_residue_n_at_5A": len({a for a, _ in local_pairs}),
                    "contact_DNMT3A_positions": ";".join(str(623 + a) for a, _ in sorted(local_pairs)),
                })
        metric_rows.append(base)

    write_tsv(OUT / "per_prediction_confidence.tsv", metric_rows)
    write_tsv(OUT / "focal_position_complex_geometry.tsv", focal_rows)

    targets = sorted({row["target"] for row in status})
    seeds = sorted({int(row["seed"]) for row in status})
    within_rows = []
    monomer_within_max = defaultdict(dict)
    complex_within_max = {}
    for target in targets:
        for domain, (start, end) in DOMAINS.items():
            comparisons = []
            for seed_a, seed_b in itertools.combinations(seeds, 2):
                left = predictions[("monomer", target, seed_a)]
                right = predictions[("monomer", target, seed_b)]
                human_positions = []
                for human in range(start, end + 1):
                    mapping = crosswalk[target].get(human)
                    if not mapping:
                        continue
                    index = mapping["monomer"] - 1
                    if left["plddt"][index] < .70 or right["plddt"][index] < .70:
                        continue
                    if target in ANCESTORS and mapping["pp"] < .80:
                        continue
                    human_positions.append(human)
                left_coords = [left["chains"]["A"][crosswalk[target][h]["monomer"] - 1]["CA"].coord for h in human_positions]
                right_coords = [right["chains"]["A"][crosswalk[target][h]["monomer"] - 1]["CA"].coord for h in human_positions]
                rmsd, *_ = superpose(left_coords, right_coords)
                comparisons.append(rmsd)
                within_rows.append({
                    "mode": "monomer", "target": target, "region": domain, "seed_a": seed_a, "seed_b": seed_b,
                    "mapped_high_confidence_n": len(human_positions), "rmsd_A": f"{rmsd:.6f}",
                    "contact_pair_jaccard": "", "contact_residue_jaccard": "",
                })
            monomer_within_max[target][domain] = max(comparisons)

        comparisons = []
        for seed_a, seed_b in itertools.combinations(seeds, 2):
            left = predictions[("complex", target, seed_a)]
            right = predictions[("complex", target, seed_b)]
            a_indices = [i for i in range(290) if left["plddt"][i] >= .70 and right["plddt"][i] >= .70]
            b_human = [h for h, mapping in crosswalk[target].items() if mapping["b"] is not None and left["plddt"][290 + mapping["b"] - 1] >= .70 and right["plddt"][290 + mapping["b"] - 1] >= .70 and (target not in ANCESTORS or mapping["pp"] >= .80)]
            ref_a = [left["chains"]["A"][i]["CA"].coord for i in a_indices]
            mov_a = [right["chains"]["A"][i]["CA"].coord for i in a_indices]
            _, _, rotation, translation = superpose(ref_a, mov_a)
            ref_b = np.asarray([left["chains"]["B"][crosswalk[target][h]["b"] - 1]["CA"].coord for h in b_human])
            mov_b = np.asarray([right["chains"]["B"][crosswalk[target][h]["b"] - 1]["CA"].coord for h in b_human])
            transformed_b = mov_b @ rotation + translation
            orientation_rmsd = float(np.sqrt(np.mean(np.sum((ref_b - transformed_b) ** 2, axis=1))))
            internal_rmsd, *_ = superpose(ref_b, mov_b)
            pairs_left = contact_pairs(left["chains"]["A"], left["chains"]["B"])
            pairs_right = contact_pairs(right["chains"]["A"], right["chains"]["B"])
            canonical_left = {(623 + a, reverse_b[target][b]) for a, b in pairs_left if b in reverse_b[target]}
            canonical_right = {(623 + a, reverse_b[target][b]) for a, b in pairs_right if b in reverse_b[target]}
            residues_left = {("A", a) for a, _ in canonical_left} | {("B", b) for _, b in canonical_left}
            residues_right = {("A", a) for a, _ in canonical_right} | {("B", b) for _, b in canonical_right}
            comparisons.append(orientation_rmsd)
            within_rows.append({
                "mode": "complex", "target": target, "region": "DNMT3A_aligned_DNMT3L_CTD", "seed_a": seed_a, "seed_b": seed_b,
                "mapped_high_confidence_n": len(b_human), "rmsd_A": f"{orientation_rmsd:.6f}",
                "internal_DNMT3L_CTD_rmsd_A": f"{internal_rmsd:.6f}",
                "contact_pair_jaccard": f"{jaccard(canonical_left, canonical_right):.6f}",
                "contact_residue_jaccard": f"{jaccard(residues_left, residues_right):.6f}",
            })
        complex_within_max[target] = max(comparisons)
    write_tsv(OUT / "within_target_seed_dispersion.tsv", within_rows)

    contrast_rows = []
    for target_a, target_b in itertools.combinations(targets, 2):
        for domain, (start, end) in DOMAINS.items():
            values = []
            ns = []
            for seed in seeds:
                left = predictions[("monomer", target_a, seed)]
                right = predictions[("monomer", target_b, seed)]
                human_positions = []
                for human in range(start, end + 1):
                    ma, mb = crosswalk[target_a].get(human), crosswalk[target_b].get(human)
                    if not ma or not mb:
                        continue
                    ia, ib = ma["monomer"] - 1, mb["monomer"] - 1
                    if left["plddt"][ia] < .70 or right["plddt"][ib] < .70:
                        continue
                    if target_a in ANCESTORS and ma["pp"] < .80:
                        continue
                    if target_b in ANCESTORS and mb["pp"] < .80:
                        continue
                    human_positions.append(human)
                left_coords = [left["chains"]["A"][crosswalk[target_a][h]["monomer"] - 1]["CA"].coord for h in human_positions]
                right_coords = [right["chains"]["A"][crosswalk[target_b][h]["monomer"] - 1]["CA"].coord for h in human_positions]
                rmsd, *_ = superpose(left_coords, right_coords)
                values.append(rmsd)
                ns.append(len(human_positions))
            threshold = 2 * max(monomer_within_max[target_a][domain], monomer_within_max[target_b][domain])
            contrast_rows.append({
                "mode": "monomer", "target_a": target_a, "target_b": target_b, "region": domain,
                "matched_seed_rmsd_mean_A": f"{np.mean(values):.6f}", "matched_seed_rmsd_min_A": f"{min(values):.6f}",
                "matched_seed_rmsd_max_A": f"{max(values):.6f}", "mapped_high_confidence_n_min": min(ns),
                "two_x_larger_within_seed_max_A": f"{threshold:.6f}",
                "all_three_between_exceed_threshold": str(all(value > threshold for value in values)),
            })

        values = []
        ns = []
        for seed in seeds:
            left = predictions[("complex", target_a, seed)]
            right = predictions[("complex", target_b, seed)]
            a_indices = [i for i in range(290) if left["plddt"][i] >= .70 and right["plddt"][i] >= .70]
            human_positions = []
            for human in sorted(set(crosswalk[target_a]) & set(crosswalk[target_b])):
                ma, mb = crosswalk[target_a][human], crosswalk[target_b][human]
                if ma["b"] is None or mb["b"] is None:
                    continue
                if left["plddt"][290 + ma["b"] - 1] < .70 or right["plddt"][290 + mb["b"] - 1] < .70:
                    continue
                if target_a in ANCESTORS and ma["pp"] < .80:
                    continue
                if target_b in ANCESTORS and mb["pp"] < .80:
                    continue
                human_positions.append(human)
            ref_a = [left["chains"]["A"][i]["CA"].coord for i in a_indices]
            mov_a = [right["chains"]["A"][i]["CA"].coord for i in a_indices]
            _, _, rotation, translation = superpose(ref_a, mov_a)
            ref_b = np.asarray([left["chains"]["B"][crosswalk[target_a][h]["b"] - 1]["CA"].coord for h in human_positions])
            mov_b = np.asarray([right["chains"]["B"][crosswalk[target_b][h]["b"] - 1]["CA"].coord for h in human_positions])
            transformed_b = mov_b @ rotation + translation
            values.append(float(np.sqrt(np.mean(np.sum((ref_b - transformed_b) ** 2, axis=1)))))
            ns.append(len(human_positions))
        threshold = 2 * max(complex_within_max[target_a], complex_within_max[target_b])
        contrast_rows.append({
            "mode": "complex", "target_a": target_a, "target_b": target_b, "region": "DNMT3A_aligned_DNMT3L_CTD",
            "matched_seed_rmsd_mean_A": f"{np.mean(values):.6f}", "matched_seed_rmsd_min_A": f"{min(values):.6f}",
            "matched_seed_rmsd_max_A": f"{max(values):.6f}", "mapped_high_confidence_n_min": min(ns),
            "two_x_larger_within_seed_max_A": f"{threshold:.6f}",
            "all_three_between_exceed_threshold": str(all(value > threshold for value in values)),
        })
    write_tsv(OUT / "between_target_contrasts.tsv", contrast_rows)
    print(f"Wrote sentinel analysis to {OUT}")


if __name__ == "__main__":
    main()
