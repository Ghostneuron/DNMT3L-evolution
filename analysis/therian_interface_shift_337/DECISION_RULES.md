# Therian-stem DNMT3L interface-shift decision rules

Frozen on 2026-11-04 before generating the complete Amniota-to-Theria substitution inventory.

## Primary estimand

- Contrast: the amino-acid state reconstructed at the Amniota ancestor versus the state reconstructed at the Theria ancestor.
- Coordinate system: human DNMT3L NP_037501.2 positions 1-387.
- Primary topology: the branch-length-optimized NCBI-taxonomy topology used in the validated switching-helix ASR.
- Topology sensitivity: the validated DNMT3L gene-tree topology.
- Biological question: whether robust therian-stem substitutions are concentrated in prespecified DNMT3L functional interfaces more than expected among eligible human-coordinate positions.

## Substitution gate

A position is `ROBUST_THERIAN_STEM_CHANGE` only when all of the following hold:

1. Amniota and Theria each have maximum marginal posterior probability at least 0.80 on both topologies.
2. The reconstructed amino-acid states differ between Amniota and Theria.
3. The exact Amniota-to-Theria state pair agrees between the taxonomy and gene-tree topologies.
4. The position maps unambiguously to NP_037501.2.
5. The position passes the alignment-robustness gate below.

Changes meeting items 1-4 but not item 5 are retained as `ALIGNMENT_SENSITIVE_CHANGE` and cannot support enrichment or functional-shift claims.

## Alignment-robustness gate

- Global MAFFT and clade-specific MACSE states are compared by human reference coordinate.
- A position must have at least 80% resolved occupancy in Amphibia, Sauropsida, Marsupialia, and Placentalia in the global human-coordinate alignment.
- It must also have at least 80% resolved occupancy in the Amphibia, Sauropsida, and Placentalia MACSE alignments.
- Among taxa resolved by both methods, exact amino-acid agreement must be at least 95% overall and at least 90% within each clade having at least four paired calls.
- Marsupialia has no independent clade-specific MACSE alignment and is therefore occupancy-gated but not cross-aligner-confirmed. This limitation is reported.

## Prespecified regions and interfaces

- ADD domain: NP41-173.
- C-terminal methyltransferase-like region: NP178-380.
- Switching helix: NP341-355.
- Adjacent loop: NP356-360.
- Experimentally motivated atlas: the 12 positions frozen in `metadata/functional_site_manifest.tsv`.
- DNMT3A-contacting interface: any DNMT3L residue with a heavy atom within 5.0 Angstrom of a DNMT3A heavy atom in either resolved DNMT3L copy of PDB 9MPP.
- Nucleosome-contacting interface: any DNMT3L residue with a heavy atom within 5.0 Angstrom of histone or DNA heavy atoms in either resolved DNMT3L copy of PDB 9MPP.
- Structural coordinates are mapped through the verified NP_037501.2-to-9MPP crosswalk. Contact sets are defined before examining which positions changed.

## Enrichment tests

- Primary structural test: Fisher exact test for robust changes in the union of DNMT3A-contacting and nucleosome-contacting positions versus other eligible, structure-resolved positions.
- Primary regional test: Fisher exact test for robust changes in the switching helix versus other eligible C-terminal positions.
- Secondary tests: ADD, C-terminal, adjacent loop, DNMT3A contact, nucleosome contact, and the 12-site atlas. Benjamini-Hochberg correction is applied across all reported enrichment tests.
- A 100,000-draw permutation sensitivity samples the observed number of robust substitutions from eligible positions while preserving broad region (ADD, C-terminal, other) and structure-resolution counts. The statistic is the number falling in the interface union.
- Overlapping regions are not independent biological replications.

## Decision labels

- `GO_COORDINATED_INTERFACE_REMODELING`: at least two robust therian-stem substitutions fall in the prespecified interface union and the primary interface enrichment has BH q < 0.05 with an empirical permutation p < 0.05.
- `NARROW_INTERFACE_HYPOTHESIS`: one or more robust interface substitutions are present, but coordinated enrichment does not pass both thresholds.
- `DESCRIPTIVE_BRANCH_CHANGE`: robust substitutions exist but no interface concentration is supported.
- `STOP_FUNCTIONAL_SHIFT_CLAIM`: no robust therian-stem substitution survives the full gate.

Even a `GO` result supports coordinated sequence remodeling consistent with specialization. It does not establish altered binding, methylation, genomic imprinting, or adaptation.

## Structural-modeling gate

Ancestral complex modeling proceeds only if NP349 or at least one other switching-helix/interface position is a robust change. Static models will be treated as structural hypotheses, not functional measurements. Reciprocal substitutions must be evaluated in both reconstructed ancestral backgrounds, and whole-helix comparisons must be kept separate from single-site comparisons.

## Stop rules

- Do not add post hoc branch-site tests to the completed selection family.
- Do not treat the taxonomy and gene-tree results as independent replication.
- Do not infer imprinting coevolution from the single therian origin event.
- Do not call recurrent residue gains adaptive convergence without independent phenotypic and functional evidence.
