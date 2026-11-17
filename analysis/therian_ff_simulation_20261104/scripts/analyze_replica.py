#!/usr/bin/env python3
"""Calculate prespecified FF-interface metrics for one completed replica."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd


def atoms_for(topology, chain_id: str, residue_id: int, names: set[str]) -> list[int]:
    found = []
    for atom in topology.atoms:
        residue = atom.residue
        if (
            residue.chain.chain_id == chain_id
            and residue.resSeq == residue_id
            and atom.name in names
        ):
            found.append(atom.index)
    if not found:
        raise RuntimeError(f"No atoms found for {chain_id}:{residue_id} {sorted(names)}")
    return found


def minimum_pair_distance(
    traj, first: list[int], second: list[int], periodic: bool
) -> np.ndarray:
    pairs = np.array([(i, j) for i in first for j in second], dtype=int)
    distances = md.compute_distances(traj, pairs, periodic=periodic)
    return distances.min(axis=1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replica", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    trajectory_path = args.replica / "trajectory.dcd"
    completion_path = args.replica / "COMPLETE.json"
    if not completion_path.exists():
        raise RuntimeError(f"Replica has not passed completion gate: {args.replica}")
    traj = md.load(str(trajectory_path), top=str(args.topology))
    top = traj.topology
    periodic = traj.unitcell_lengths is not None

    protein_atoms = top.select("protein")
    backbone_atoms = top.select("protein and backbone")
    if periodic:
        traj.image_molecules(inplace=True)
    traj.superpose(traj, 0, atom_indices=backbone_atoms)

    aromatic_names = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}
    dnmt3l_261 = atoms_for(top, "B", 261, aromatic_names)
    dnmt3a_728 = atoms_for(top, "A", 728, aromatic_names)
    dnmt3a_768 = atoms_for(top, "A", 768, aromatic_names)
    arg_names = {"NE", "CZ", "NH1", "NH2"}
    # DNMT3A1 manuscript coordinates are four residues higher than the 2QRV
    # construct in this region: R729/E725/G726 correspond to PDB R725/E721/G722.
    dnmt3a_r725 = atoms_for(top, "A", 725, arg_names)
    dnmt3a_e721 = atoms_for(top, "A", 721, {"OE1", "OE2"})
    dnmt3a_g722 = atoms_for(top, "A", 722, {"N", "CA", "C", "O"})
    focal_291 = atoms_for(top, "B", 291, {"OE1", "OE2", "OG1"})
    focal_293 = atoms_for(top, "B", 293, {"OE1", "OE2", "OD1", "OD2"})

    metrics = pd.DataFrame(
        {
            "frame": np.arange(traj.n_frames),
            "time_ps": traj.time,
            "min_B261_A728_nm": minimum_pair_distance(
                traj, dnmt3l_261, dnmt3a_728, periodic
            ),
            "min_B261_A768_nm": minimum_pair_distance(
                traj, dnmt3l_261, dnmt3a_768, periodic
            ),
            "min_B291_A725_R729ref_nm": minimum_pair_distance(
                traj, focal_291, dnmt3a_r725, periodic
            ),
            "min_B291_A721_E725ref_nm": minimum_pair_distance(
                traj, focal_291, dnmt3a_e721, periodic
            ),
            "min_B293_A725_R729ref_nm": minimum_pair_distance(
                traj, focal_293, dnmt3a_r725, periodic
            ),
            "min_B293_A721_E725ref_nm": minimum_pair_distance(
                traj, focal_293, dnmt3a_e721, periodic
            ),
            "min_B293_A722_G726ref_nm": minimum_pair_distance(
                traj, focal_293, dnmt3a_g722, periodic
            ),
            "protein_backbone_rmsd_nm": md.rmsd(
                traj, traj, 0, atom_indices=backbone_atoms
            ),
        }
    )
    for column in [c for c in metrics if c.startswith("min_")]:
        metrics[column.replace("min_", "contact05_")] = (
            metrics[column] < 0.5
        ).astype(int)

    active_loop_ca = top.select(
        "chainid 0 and resid 81 to 102 and name CA"
    )
    # The explicit chain/resSeq audit below protects against mdtraj resid indexing drift.
    active_loop_ca = np.array(
        [
            atom.index
            for atom in top.atoms
            if atom.residue.chain.chain_id == "A"
            and 700 <= atom.residue.resSeq <= 721
            and atom.name == "CA"
        ],
        dtype=int,
    )
    if len(active_loop_ca) != 22:
        raise RuntimeError(f"Expected 22 active-loop CA atoms; found {len(active_loop_ca)}")
    loop_rmsf = md.rmsf(traj, traj, 0, atom_indices=active_loop_ca)

    args.output.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.output / "frame_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "chain": "A",
            "pdb_residue": list(range(700, 722)),
            "dnmt3a1_reference_residue": list(range(704, 726)),
            "ca_rmsf_nm": loop_rmsf,
        }
    ).to_csv(args.output / "active_loop_rmsf.tsv", sep="\t", index=False)

    summary = {
        "frames": traj.n_frames,
        "atoms": traj.n_atoms,
        "protein_atoms": int(len(protein_atoms)),
        "periodic": periodic,
        "metrics": {
            column: {
                "mean": float(metrics[column].mean()),
                "median": float(metrics[column].median()),
                "sd": float(metrics[column].std(ddof=1)),
            }
            for column in metrics
            if column not in {"frame", "time_ps"}
        },
        "active_loop_mean_ca_rmsf_nm": float(loop_rmsf.mean()),
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

