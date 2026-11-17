# Sentinel structural-screen decision

Frozen interpretation: 2026-11-05

## Validated evidence

- All 42 sentinel predictions passed output and dimensional validation: seven targets, three seeds, monomer and fixed-partner complex modes.
- All seven targets passed the ADD-domain confidence gate in all three seeds.
- The monomer C-terminal confidence gate passed in at least two seeds for human, mouse, reconstructed Amniota, and reconstructed Theria. It did not pass for Alligator, Monodelphis, or Geotrypetes.
- All seven fixed-partner complexes passed the interface-confidence gate in all three seeds. Minimum target-level ipTM values ranged from 0.820 to 0.890, minimum median interface pLDDT from 0.920 to 0.942, and maximum median contact-pair PAE from 2.08 to 3.06 A.
- Interface-contact residue overlap across seeds was high (minimum target-level Jaccard 0.833-0.907), although DNMT3L orientation dispersion after DNMT3A alignment reached 2.32 A for Geotrypetes and 1.98 A for mouse.

## Primary ancestral contrast

The Amniota-to-Theria comparison did not exceed the prespecified two-times-within-seed threshold:

- ADD monomer RMSD: mean 0.460 A; threshold 0.680 A.
- C-terminal monomer RMSD: mean 0.518 A; threshold 1.401 A.
- DNMT3A-aligned DNMT3L C-terminal RMSD: mean 1.962 A; threshold 3.106 A.

Therefore, these predictions do not support a reproducible global domain or complex-orientation shift on the therian stem.

## Focal ancestral states

- NP261 Y-to-F retained the same four predicted DNMT3A contacts (human DNMT3A positions 724, 729, 732, and 772) in all six ancestral predictions.
- NP291 E-to-T showed seed-variable marginal contact behavior and does not define a reproducible contact switch.
- NP293 E-to-D changed some predicted local contacts, but the global complex difference did not exceed seed dispersion; this remains a local structural hypothesis, not evidence of altered binding.
- NP313 H-to-C remained approximately 30 A from DNMT3A and is not tested mechanistically by this binary complex.
- NP349 K-to-Q remained approximately 12-17 A from DNMT3A with no direct DNMT3A contact. This binary complex does not test the proposed nucleosome-facing role.

## Expansion decision

- **GO: balanced monomer expansion to all 24 inputs.** Purpose: determine whether the reproducible Geotrypetes ADD divergence is shared by the two salamanders or is caecilian-specific, and map other lineage-local domain deviations. This is an exploratory localization screen, not a test of adaptation.
- **NO-GO: 24-input complex expansion at this stage.** The primary ancestral complex contrast did not exceed seed variability, so more species would add computation without resolving the therian-mechanism hypothesis.
- **NO-GO: new MD from the Boltz screen.** No Boltz contrast meets the MD-admission rule. The already-running 2QRV-based ancestral-interface MD remains an independent, prespecified analysis.

## Claim boundary

The sentinel results support confidence-qualified structural similarity or divergence under one prediction policy. They do not establish altered binding affinity, nucleosome engagement, methylation activity, genomic imprinting, or adaptive evolution.
