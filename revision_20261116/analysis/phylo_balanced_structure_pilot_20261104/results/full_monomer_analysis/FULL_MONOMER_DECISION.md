# Full balanced-panel monomer decision

Frozen interpretation: 2026-11-05

## Scope

The analysis includes 24 phylogenetically balanced DNMT3L inputs and three Boltz-2 seeds per input. Seeds are technical reproducibility checks, not biological replicates. Comparisons use human-coordinate crosswalks, pLDDT at least 0.70, and ancestral posterior probability at least 0.80 where applicable.

## Confidence gates

- ADD: 23/24 targets passed the target gate. Dog passed only one of three seed-level confidence gates and is excluded from ADD comparisons.
- C-terminal monomer region: 10/24 targets passed. All six sauropsids failed, as did Geotrypetes and several mammalian tips. Broad vertebrate CTD inference from isolated monomers is therefore not supported.
- The strong CTD confidence previously observed in fixed-partner complexes indicates prediction-context sensitivity; low monomer confidence is not evidence that the biological domain is unstructured.

## Amphibian result

- The three confidence-qualified amphibian ADD comparisons ranged from 0.478 to 0.613 A mean matched-seed RMSD.
- None exceeded two times the larger within-target seed dispersion.
- Geotrypetes was not distinguishable from Lissotriton (0.613 A; threshold 1.171 A) or Triturus (0.542 A; threshold 1.012 A).
- The broader panel therefore does **not** support a caecilian-specific ADD fold. The earlier Geotrypetes contrasts reflected comparisons with selected non-amphibian tips plus its relatively low within-seed dispersion.
- The two confidence-qualified salamander CTDs were also similar (0.373 A) and did not exceed seed dispersion. Geotrypetes CTD did not pass its target gate.

## Lineage structure

- ADD distances were modest overall. Within-lineage medians were 0.542 A in Amphibia, 0.455 A in Sauropsida, 0.436 A in Marsupialia, and 0.560 A in Placentalia.
- Amphibian-marsupial ADD distances were larger descriptively (median 0.914 A), but only 4/12 individual pairs passed the seed-dispersion rule. This is not a clean clade-wide discontinuity.
- Reconstructed Amniota and Theria ADD structures differed by 0.460 A and did not pass the contrast gate. Their CTDs differed by 0.518 A and also did not pass.
- Only two ancestral-versus-extant ADD pairs passed: Geotrypetes versus Theria and Sarcophilus versus Amniota. Isolated tip contrasts do not define an ancestral functional transition.
- Confidence-qualified CTD variation was substantial within Placentalia; human-mouse remained a reproducible contrast. This within-clade heterogeneity argues against interpreting CTD geometry as a simple imprinting milestone.

## Decision

- **Stop additional monomer prediction:** the balanced panel is complete and answers the amphibian-localization question.
- **Keep the 24-input complex expansion stopped:** the primary Amniota-Theria complex contrast did not exceed seed variability.
- **Do not claim a therian-stem structural switch:** neither domain nor fixed-partner complex predictions support it.
- **Do not infer altered affinity, nucleosome binding, methylation, imprinting, or adaptation from these models.**
- **Use the structural screen as a rigorous negative result:** DNMT3L domain folds are broadly conserved under this prediction policy, while local sequence states remain candidates for biochemical—not global-fold—effects.
- **Next evidence gate:** complete and analyze the independent replicated 2QRV-based MD trajectories, then decide whether any local interface hypothesis survives trajectory-level testing.

## Manuscript implication

These results improve rigor but do not create a new causal-imprinting result. They support narrowing the manuscript from a proposed global structural transition to conserved domain architecture with lineage-structured local substitutions whose functional consequences remain unresolved.
