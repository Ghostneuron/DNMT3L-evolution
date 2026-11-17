# Frozen decision rules for the 337-GeneID selection reanalysis

Frozen on 2026-11-02 before inspecting any HyPhy output from the updated sequence set.

## Scope and estimand

- Taxon scope: placental DNMT3L orthologues that passed protein-CDS identity, sequence QC, and DNMT3A/B/L family-tree validation.
- Primary alignment: the strict exact-codon consensus across MACSE, MAFFT, PRANK, MUSCLE5, and COBALT, restricted to human-coordinate positions with at least 80% post-consensus occupancy.
- Primary tree: one codon gene-tree topology inferred from the non-deduplicated strict-consensus alignment and pruned to each deduplicated selection input. Branch lengths are re-estimated by HyPhy.
- Primary site model: standard all-site MEME (`multiple-hits None`, all branches).
- Complementary model: FUBAR for pervasive selection. FUBAR and MEME answer different questions and are not treated as mutually validating replicates.
- MACSE, MAFFT, PRANK, MUSCLE5, and COBALT analyses are alignment sensitivities, not five independent biological tests.

## Coordinates and multiplicity

- Every result is mapped back to human DNMT3L NP_037501.2 coordinates with the frozen crosswalk.
- Standard MEME is run over every retained codon. Benjamini-Hochberg correction is applied across every MEME-tested codon within an alignment; `q < 0.05` is the primary discovery threshold.
- Nominal MEME `p < 0.05` without full-alignment `q < 0.05` is not evidence of episodic selection.
- FUBAR posterior probability `>= 0.90` is reported as model-specific support for pervasive diversifying or purifying selection; it is not described as a MEME-significant result.
- The prespecified functional-site family contains 12 unique human positions: 88, 90, 105, 106, 107, 204, 349, 351, 352, 355, 358, and 377. These sites receive focused reporting, but no reduced multiple-testing denominator will be used to rescue a result that fails the all-site MEME correction.
- h358 is exploratory because it originated in earlier analyses. It must meet the same all-site correction and robustness rules as any other site and will not be presented as an independent confirmation.

## Robustness labels

- **Robust episodic-selection support:** primary-consensus MEME `q < 0.05`, nominal `p < 0.05` in at least four of five raw-aligner sensitivities, and no reversal under the prespecified multiple-hit or topology sensitivity.
- **Primary-only/model-sensitive support:** primary-consensus MEME `q < 0.05` but the raw-aligner, multiple-hit, or topology criterion above is not met. This is reported as alignment/model dependent, not a robust site.
- **Sensitivity-only support:** one or more raw alignments meet a threshold but the strict-consensus primary does not. This cannot support a positive-selection claim.
- **No corrected support:** primary-consensus MEME `q >= 0.05`, regardless of nominal p-values.
- **Pervasive-selection support:** FUBAR posterior `>= 0.90` in the primary consensus and the same direction in at least four of five raw-aligner sensitivities. This label remains specific to FUBAR.

## Prespecified sensitivities

1. Repeat standard all-site MEME and FUBAR on each raw aligner after the same human-coordinate 80% occupancy rule and exact-sequence deduplication.
2. Run MEME with `Double+Triple` multiple hits for the 12 prespecified functional sites and for any full-alignment standard-MEME discovery. The standard all-site MEME/BH result remains primary.
3. Repeat focal interpretation on an NCBI-taxonomy/species-tree topology pruned to the same tips. The codon gene tree remains primary; topology disagreement downgrades a result.
4. Quantify codon saturation and inspect recombination/breakpoint sensitivity before interpreting isolated high-dN/dS sites.

## Claim boundaries

- Selection models test substitution patterns; they do not establish biochemical effects of Q/K or any other residue state.
- A conserved or phylogenetically clustered state can be biologically informative without positive-selection significance.
- No association with lifespan, methylation phenotype, or organismal trait will be inferred from this analysis.
- Failure to detect selection is not evidence that a site is functionally irrelevant.

