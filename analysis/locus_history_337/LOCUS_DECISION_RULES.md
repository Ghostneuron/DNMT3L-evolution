# Vertebrate DNMT3L locus-evidence rules

These rules were formalized after the preliminary annotation and neighborhood screens but before construction of the canonical locus matrix. They grade evidence; they do not convert a missing annotation into a gene-loss claim.

## Presence grades

- `A_CONFIRMED_PRESENT`: an exact protein-CDS pair passes sequence/domain QC and the broad DNMT3-family tree places the record in the validated DNMT3L analysis set.
- `B1_SUPPORTED_PARTIAL_SYNTENIC`: an exact translated protein-CDS pair belongs to the NCBI DNMT3L ortholog group, is closer to validated DNMT3L sequences than to DNMT3A/B/C in the targeted human-coordinate tree, and a DNMT3L-core translated match occurs in the conserved neighborhood.
- `B2_SUPPORTED_PARTIAL_NONSYNTENIC`: the same annotation and targeted-tree evidence passes, but the conserved-neighborhood interval does not recover the annotated locus. This supports presence while showing that the tested marker block is not informative for that species.

## Non-detection grades

- `C_NOT_DETECTED_CONSERVED_NEIGHBORHOOD`: no DNMT3L annotation is recovered; at least two marker genes define a compact assembled interval; and neither the representative nor full-family translated search recovers a DNMT3L-core match at E <= 1e-10. This is a candidate absence, not a demonstrated loss.
- `D_UNRESOLVED_MARKER_CENTERED`: fewer than two markers define a conserved block. A negative marker-centered search is weak evidence and remains unresolved.

## Core-overlap safeguard

A translated hit to a DNMT3L-labelled subject is not accepted automatically. The subject interval must overlap the part of that record that maps to human NP_037501.2. A hit is considered core-compatible only when at least 40 subject residues and at least 30% of the aligned segment map to human-reference-occupied columns. This rule excludes matches confined to lineage-specific extensions or fused neighboring proteins.

## Loss inference

- No single-species `C` or `D` result is called a gene loss.
- Concordant `C` results in two bird species may support a candidate avian absence pattern, but not the timing or mechanism of loss.
- A complete vertebrate origin/loss reconstruction requires broader assembly sampling and explicit reconciliation; it is not inferred from this representative panel alone.
