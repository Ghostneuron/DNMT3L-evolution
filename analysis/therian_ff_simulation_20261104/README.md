# DNMT3L therian-stem FF-interface simulation

## Scientific contrast

This analysis compares two reconstructed DNMT3L interface states against an
otherwise identical DNMT3A catalytic-domain partner:

- Amniote-like DNMT3L: Y261/E291/E293
- Therian-like DNMT3L: F261/T291/D293

The structural template is the first DNMT3A-DNMT3L heterodimer (chains A and B)
from PDB 2QRV. PDB 2QRV contains F261/N291/E293 in DNMT3L, so both primary
states are generated explicitly from the same coordinates.

## Prespecified interpretation

The primary question is whether replicated trajectories show a reproducible
difference in:

1. aromatic packing between DNMT3L 261 and DNMT3A F728/F768 (PDB numbering;
   corresponding to human DNMT3A1 F732/F772 in the project crosswalk);
2. contacts between DNMT3L 291/293 and DNMT3A1 E725/G726/R729
   (2QRV E721/G722/R725); and
3. dynamics of the adjacent DNMT3A active-site loop.

These simulations can support a trajectory-level structural hypothesis. They
cannot establish altered binding affinity, catalytic activity, imprinting, or
adaptation. Alchemical binding calculations are gated on stable, reproducible
endpoint simulations.

## Staged run contract

1. Validate the local OpenMM environment and source structure.
2. Prepare and minimize both variants with the same force field and solvent
   settings.
3. Run a short smoke trajectory for each variant.
4. Benchmark throughput and memory.
5. Launch matched independent replicas only if both smoke runs pass.
6. Analyze the frozen contact and dynamics readouts.
7. Decide whether an alchemical free-energy calculation is technically
   justified, with special handling of the charge-changing E291T mutation.

## Phase-1 production design

- Explicit TIP3P-FB solvent, 0.15 M NaCl, PME, 1.0-nm cutoff, 300 K and 1 bar.
- Protein-hydrogen mass repartitioning to 3 amu; water masses are unchanged.
- Four-femtosecond Langevin integration after matched 4-fs smoke validation.
- Three independent 5-ns production replicas per ancestral endpoint, each
  preceded by 100 ps of equilibration.
- The first four jobs use two CPU threads each; the final matched pair uses
  five threads each. Seeds are paired by replicate number but trajectories are
  independently randomized.
- Phase 1 is a convergence and effect-direction gate. No trajectory-level
  difference will be interpreted without cross-replica agreement.

## Retired implicit-solvent branch

A matched GBn2/NoCutoff screen was attempted and rejected. The therian endpoint
developed a NaN coordinate during bounded minimization, and the still-running
amniote minimization was stopped because an asymmetric protocol cannot support
the comparison. These failed preparations are retained under
`prepared_implicit/` as provenance and must not be analyzed as results.

Authoritative run state is recorded in `RUN_STATUS.md`; a process is not treated
as complete merely because its PID disappears.

