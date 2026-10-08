# Exploratory structure-informed residue survey

Date added: 2026-10-06, after inspection of the amino-acid constraint summaries and before inspection of the broadened HyPhy selection outputs.

## Rationale

The 2026 DNMT3A2-DNMT3L cryo-EM study identifies DNMT3L UniProt/PDB Q348 and Q351 as the residues homologous to the acidic-patch-interacting arginine fingers of DNMT3B3. Because NP_037501.2 contains an additional serine at position 333, these residues map to NP positions 349 and 352. Initial domain summaries indicated lineage-dependent variability across the switching helix, motivating a descriptive survey of these structurally annotated positions across all selection-ready records.

## Analysis boundary

- Map NP_037501.2 positions 341-360 to each of the 269 selection-ready proteins by the same global BLOSUM62 pairwise-alignment procedure used in the archived orthology audit.
- Report residue counts and major-state frequencies by the seven archived clade labels.
- Highlight NP349/UniProt-Q348 and NP352/UniProt-Q351 as structure-informed positions.
- Cross-check the focal state distributions against the independently generated MACSE, MAFFT, and PRANK summaries where those scopes are available.
- Treat the survey as descriptive and exploratory. It may motivate biochemical comparison of lineage-specific variants, but it cannot establish ancestral direction, adaptive selection, altered acidic-patch binding, or a methylation phenotype.

## Same-length structural-window sensitivity

Added 2026-10-06 before the broadened selection outputs were summarized:

- Compare the median entropy of the complete NP341-355 switching-helix window and NP356-360 adjacent-loop window with every complete, alignment-stable window of the same length within NP178-379.
- Report a conservative empirical upper-tail probability and apply Benjamini-Hochberg correction across every evaluable scope-region comparison.
- Treat this as a post hoc robustness analysis of the regional-variability observation. It does not convert variability into evidence of positive selection or mechanism.

## Tree-aware constraint sensitivity

Added 2026-10-06 before the broadened selection outputs were summarized:

- On each independently filtered MACSE data set, use its NCBI-taxonomy species topology and unordered-amino-acid Sankoff parsimony to calculate the minimum number of state changes per human-mapped position.
- Normalize each score by the number of resolved tips and restrict comparisons to positions meeting the cross-aligner stability rule.
- Re-evaluate ADD versus C-terminal constraint and the two structural-region contrasts, including an equal-length C-terminal window comparison for each complete focal region. This conservative sensitivity addresses shared ancestry but does not estimate substitution rates, ancestral direction, or selection.

## Coordinate correction

Added 2026-10-06 after verification against UniProt Q9UJW3 and the structural construct coordinates:

- UniProt Q9UJW3 residues 178-379 correspond to NP_037501.2 residues 178-380 because NP_037501.2 contains an additional serine at NP333.
- The C-terminal comparison interval and its equal-length background windows were therefore corrected to NP178-380 before final summaries were generated. No selection-model input or focal structural-region coordinate changed.

## Per-alignment domain-direction audit

Added 2026-10-06 after inspection of the aggregate domain-entropy summaries and while the broadened HyPhy runs were still in progress:

- Report ADD and C-terminal median entropy separately for every scope-by-aligner data set.
- Use this as a descriptive direction-of-effect audit only. Aligner agreement is technical robustness and is not counted as independent biological replication or as an additional hypothesis test.

## Alignment-overlap audit

Added 2026-10-06 after examining the number of taxa jointly resolved across the three aligners and before final selection outputs were summarized:

- Some positions had high conditional amino-acid agreement but only a small jointly resolved taxon set.
- A final alignment-stability label therefore requires both at least 80% joint taxon resolution and at least 90% exact amino-acid agreement within the joint set.
- This is a conservative eligibility correction. It cannot create a candidate and does not alter any codon-model fit.
