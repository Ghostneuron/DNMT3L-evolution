# DNMT3L regulatory deployment, matched-partner evolution, and PTM conservation

Started: 2026-11-04

This add-on asks why a deeply retained DNMT3L-class protein does not yield a single obvious sequence marker for therian genomic imprinting. It tests three non-exclusive alternatives before expanding structural modeling.

## Prespecified questions

### 1. Germline deployment

**Primary estimand:** within homologous male germ-cell states in the Murat et al. multi-species single-nucleus atlas, whether DNMT3L is detected and how its cell-state distribution compares descriptively with DNMT3A.

- Biological unit: donor/animal within species; nuclei are nested observations, not independent subjects.
- Primary endpoint: donor-aware detection and pseudobulk abundance by author-defined homologous germ-cell state.
- Confirmatory scope: processed, annotation-resolved data from the same study.
- Exploratory scope: female-germline promoter/isoform evidence from independent studies.
- Claim boundary: expression timing can support lineage-specific deployment; it cannot establish that DNMT3L caused imprinting.

### 2. Matched DNMT3A-DNMT3L evolution

**Primary estimand:** among species with independently validated DNMT3A and DNMT3L orthologues, whether prespecified interface residues show paired lineage transitions beyond the background expected from shared phylogeny and regional evolutionary rate.

- Primary DNMT3A coordinates: human DNMT3A1 E725, G726, R729, F732, and F772, sequence-verified against the structural construct before analysis.
- Primary DNMT3L coordinates: human NP_037501.2 positions 261, 291, and 293; other interface-proximal therian-stem changes are secondary.
- Primary comparison: branch-mapped state changes at prespecified contacts versus rate-, occupancy-, and region-matched non-interface positions.
- Raw across-species correlations are descriptive because both proteins share the species tree.
- Claim boundary: a significant paired pattern would prioritize coevolution; it would not demonstrate altered affinity or imprinting function.

### 3. PTM-capable residue conservation

**Primary estimand:** the lineage distribution and ancestral state history of experimentally observed mouse DNMT3L PTM-capable residues after source-to-human coordinate verification.

- Sites with direct functional testing, including mouse K238 and K412, are distinguished from mass-spectrometry-only candidates.
- Conservation of a modifiable residue does not demonstrate modification in another species.
- Motif prediction, if used, is exploratory and cannot replace empirical PTM evidence.

### 4. Comparative structures

Whole-protein or complex prediction will be admitted only after modules 1-3 nominate a species, lineage, or residue contrast with adequate sequence identity and coordinate support. Predicted local differences are hypotheses, not measurements of affinity, catalysis, or imprinting.

## Frozen decision order

1. Validate artifacts, metadata, identities, orthology, and coordinates.
2. Freeze inclusion, exclusion, and matched-background rules.
3. Run the smallest valid analysis.
4. Apply sensitivity analyses without promoting them to independent replication.
5. Update the manuscript only after an output passes its module-specific gates.

## Current scope decision

- `GO` for public-data inventory and small processed artifacts.
- `GO` for a matched-species DNMT3A-DNMT3L sequence inventory.
- `GO` for coordinate-qualified PTM state mapping.
- `STOP` using the downloaded adult single-cell panel as a primary imprint-establishment test: all nine files are processed/non-integer and developmentally stage-mismatched.
- `GO` for a descriptive adult annotation audit and for contextualizing it with the stage-resolved 2026 marmoset germ-cell study.
- `GO` for an experimentally anchored conservation test of the mouse DNMT3L G9a-contact motif L139/F140/I141.
- `WAIT` for comparative structure prediction until a sequence-based contrast passes its gates.
- `STOP` causal claims about the origin of imprinting from these observational analyses alone.

