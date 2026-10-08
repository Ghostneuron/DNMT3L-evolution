# DNMT3L broadened evolutionary reanalysis plan

Date recorded: 2026-10-06, while the broadened HyPhy analyses were running and before the final sitewise outputs were available or summarized. This internal plan was not externally preregistered.

## Primary question

Which human-referenced DNMT3L codons, if any, retain evidence of diversifying selection after independent alignment, taxonomic partitioning, multiplicity correction, and model comparison?

The previously emphasized h358 site is retained as a historical case study. It is not privileged when identifying new candidates.

## Taxonomic scopes

1. Placentalia: pooled placental sensitivity analysis.
2. Euarchontoglires and Laurasiatheria: nonoverlapping placental partitions used to assess taxonomic reproducibility.
3. Sauropsida: a deeper lineage used to assess whether mammalian patterns generalize outside mammals.

Placentalia overlaps Euarchontoglires and Laurasiatheria and therefore is not counted as an independent biological replicate of either partition.

## Alignment and sequence rules

- Build or retain independent MACSE, protein-MAFFT plus codon back-translation, and codon-aware PRANK alignments for every scope.
- Map columns to NP_037501.2 with the human sequence present as either a biological record or a removable coordinate anchor.
- Retain human-mapped codons with at least 80% resolved occupancy in each alignment separately.
- Do not impose a cross-aligner consensus mask.
- Replace unresolved non-gap codons with `NNN`; preserve complete codon gaps as `---`.
- Collapse exact nucleotide duplicates after filtering, retaining human when it belongs to a duplicate group and otherwise retaining the alphabetically first identifier.
- Infer a separate codon maximum-likelihood tree for every filtered, deduplicated alignment.

## Alignment-confidence rule

A human-referenced position is alignment-stable within a scope when at least 90% of taxa resolved in all three alignments have the same amino acid in MACSE, MAFFT, and PRANK. Occupancy and the number of jointly resolved taxa will be reported so this conditional proportion is not mistaken for complete coverage.

## Selection models and multiplicity

- FUBAR: posterior probability of pervasive diversifying or purifying selection.
- FEL: maximum-likelihood test of pervasive selection; Benjamini-Hochberg correction across every tested codon within each alignment.
- MEME: test of episodic diversifying selection; Benjamini-Hochberg correction across every tested codon within each alignment.
- BUSTED: gene-wide test for episodic diversifying selection; Benjamini-Hochberg correction across the 12 scope-by-aligner analyses.

Exploratory sitewise thresholds use FUBAR posterior probability at least 0.90 and FEL or MEME q at most 0.10. Strict results use q at most 0.05.

## Recorded evidence classes

- Alignment-robust pervasive diversification: alignment-stable, with FUBAR positive posterior probability at least 0.90 and FEL beta greater than alpha at q at most 0.10 in the same alignment, in at least two of three aligners within a scope.
- Strict pervasive diversification: the same rule with FEL q at most 0.05.
- Alignment-robust episodic diversification: alignment-stable, with MEME q at most 0.10 in at least two of three aligners within a scope.
- Strict episodic diversification: the same rule with MEME q at most 0.05.
- Alignment-robust purifying selection: alignment-stable, with FUBAR negative posterior probability at least 0.90 and FEL beta less than alpha at q at most 0.05 in the same alignment, in at least two of three aligners within a scope.
- Cross-placental recurrence: the same human-referenced site meets a robust diversifying criterion separately in Euarchontoglires and Laurasiatheria.
- Deep-lineage recurrence: a site meeting a robust mammalian criterion also meets that criterion in Sauropsida. This is interpreted as recurrent statistical evidence, not proof of the same ancestral event or mechanism.

## Secondary analyses

- Compare the new independent, unmasked results with the archived consensus-masked placental results.
- Quantify amino-acid conservation independently of selection tests by human-mapped Shannon entropy, major-state frequency, and invariant-site fraction across the three aligners; compare the defined structural domains using rank-based tests.
- Summarize robust positive and purifying sites by the human ADD region (41-173), C-terminal methyltransferase-like region (178-379), and remaining segments.
- Map results to the experimentally defined DNMT3L switching helix (PDB/UniProt E340-K354, corresponding to NP_037501.2 positions 341-355 after the NP-specific serine insertion) and its immediately adjacent NP positions 356-360. Test enrichment only against alignment-stable eligible sites and treat this narrow-region analysis as secondary.
- Test domain enrichment only when a category contains enough sites for an interpretable contingency table; correct across tested domains.
- Run BUSTED with its alignment-error sink as a sensitivity analysis for any scope with robust standard BUSTED support, using MACSE as the representative alignment recorded before the final outputs were summarized.
- Build NCBI-taxonomy species topologies for the four MACSE datasets and repeat FUBAR, FEL, MEME, and BUSTED on those topologies. Species-tree results may qualify support obtained under the alignment-derived gene tree but may not create a primary site candidate.
- Run focused multiple-hit site sensitivities only for sites that satisfy a primary robust diversifying criterion. These tests may qualify a candidate but cannot create one.

## Interpretation guardrails

- Aligner agreement is a robustness check, not independent biological replication.
- A site must satisfy the full criterion; isolated nominal p-values or one-model calls are descriptive only.
- Structural location provides experimental context but does not establish a molecular mechanism.
- If no site survives the primary rules, the principal result will be the instability of site localization despite any robust gene-wide evidence, rather than a rescued single-residue claim.

## Dated diagnostic addition

Added 2026-10-06 while the broadened HyPhy runs were still in progress and before their final outputs were summarized:

- Recalculate pairwise codon-saturation diagnostics separately for all 12 independently filtered, duplicate-collapsed alignments.
- Evaluate every sequence pair or a deterministic sample of 3,000 pairs per alignment using third-codon-position transition/transversion behavior, Kimura two-parameter distance, and pairwise NG86 dN and dS.
- Use the diagnostic only to describe multiple-hit risk and qualify interpretation. It cannot create a site or gene-wide selection claim.

## Dated coordinate correction

Added 2026-10-06 after checking the NP_037501.2-to-UniProt Q9UJW3 sequence crosswalk and before final domain summaries were generated:

- The structural construct annotated as UniProt Q9UJW3 residues 178-379 maps to NP_037501.2 residues 178-380 because NP_037501.2 has an additional serine at NP333.
- Correct the C-terminal methyltransferase-like interval used for domain and equal-length-window summaries from NP178-379 to NP178-380. The archived sequence-quality gate remains reported as originally applied to NP178-379.
- This correction does not alter alignments, occupancy filtering, phylogenies, site coordinates, selection-model inputs, selection thresholds, or the NP341-360 structural-region definitions.

## Dated alignment-overlap correction

Added 2026-10-06 after auditing the joint-resolution counts and before final sitewise selection outputs were summarized:

- The original conditional concordance rule could label a position stable when only a small subset of shared taxa was resolved by all three aligners.
- Require at least 80% of the taxa shared among the three alignments to have resolved calls in all three, in addition to at least 90% exact amino-acid agreement within that joint set.
- This conservative coverage floor changes only downstream eligibility and robustness labels. It does not change independent occupancy masks, alignments, selection inputs, trees, HyPhy fits, p values, posterior probabilities, or multiplicity families.

## Dated conditional-site clarification

Added 2026-10-06 before any conditional MEME output was generated:

- Limit the double-plus-triple-nucleotide MEME sensitivity to sites meeting a strict alignment-robust diversifying criterion at `q <= 0.05`, rather than the broader exploratory `q <= 0.10` tier.
- Include historical h358 in each scope as a diagnostic audit case even when it does not meet the strict criterion.
- These conditional results may weaken or qualify a primary candidate but cannot create one; the h358 exception does not restore its earlier privileged status.
