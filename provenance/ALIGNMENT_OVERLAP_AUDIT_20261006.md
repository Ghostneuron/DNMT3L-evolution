# Cross-aligner taxon-overlap audit

Date: 2026-10-06

## Finding

The initial concordance rule asked whether at least 90% of taxa resolved by MACSE, MAFFT, and PRANK had the same amino acid. That conditional proportion can be misleading when the three aligners resolve largely different taxa at a position. For example, Laurasiatheria NP332 initially had exact agreement among the one taxon jointly resolved across all three aligners, despite 119 shared taxon identifiers. Several sauropsid positions had no jointly resolved taxon.

This was detected while the broadened HyPhy jobs were still running and before final sitewise selection outputs were summarized.

## Correction

A final alignment-stable position must satisfy both conditions:

1. At least 80% of taxa shared among the three alignments have a resolved amino-acid call in all three.
2. At least 90% of that jointly resolved set has the same amino acid in all three aligners.

The 80% floor matches the independent per-alignment occupancy threshold and prevents high conditional concordance from a small intersection. Aligner agreement remains a technical robustness check, not biological replication.

## Effect on eligibility

Under the initial agreement-only rule, stable-position counts were 381 in Placentalia, 376 in Euarchontoglires, 380 in Laurasiatheria, and 332 in Sauropsida. Adding the joint-resolution floor gives 366, 368, 377, and 309 positions, respectively, out of 387 human-reference positions evaluated per scope.

The correction changes only downstream eligibility, domain summaries, enrichment analyses, and robust-site labels. It does not alter an alignment, occupancy mask, duplicate collapse, selection input, phylogeny, HyPhy output, p value, posterior probability, or multiple-testing family.

## Focal structural positions

- NP349/structural Q348 passes both thresholds in all four scopes. Joint-resolution coverage ranges from 93.4% to 97.6%, and exact concordance within the joint set ranges from 97.3% to 100%.
- NP352/structural Q351 has 93.8% exact concordance in Sauropsida among the jointly resolved taxa, but joint coverage is 76.2% (32/42), below the final 80% floor. Its sauropsid state distribution is therefore reported as qualified descriptive evidence and cannot support a robust codon-selection call.
- Historical h358 passes both thresholds in all four scopes. Joint coverage ranges from 92.4% to 95.2%, and exact concordance ranges from 96.2% to 99.1%.
