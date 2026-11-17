# Project status

Updated: 2026-11-04

## Current decision

`GO` for qualified sequence, PTM, and experimentally anchored interaction-motif inventories. `STOP` using the adult cross-species testis atlas as the primary regulatory test. `GO` for using stage-matched primate germ-cell and oocyte-promoter studies as published context. `WAIT` for new broad structure prediction.

## Module gates

| Module | Artifact | Metadata/identity | Scale/compute | Current decision |
|---|---|---|---|---|
| Multi-species male germline expression | Nine H5AD files PASS | DNMT3A 9/9; DNMT3L 6/9 testable | All matrices processed/non-integer; all samples adult or post-juvenile | STOP as primary test; retain descriptive audit only |
| Stage-matched primate germline deployment | Kojima et al. 2026 PASS | Marmoset developmental stages and public accessions PASS | Published scRNA/scBS and immunofluorescence | Use as prior evidence, not a new discovery |
| Female-germline promoter/isoform evolution | Published comparative study PASS | Study-level evidence PASS | Gene-level H5AD cannot test isoform usage | Treat the Eumuroida oocyte promoter as prior evidence; seek transcript-resolved data before extension |
| DNMT3A-DNMT3L coevolution | 265 matched species PASS | Five DNMT3A contacts and three DNMT3L sites PASS | PASS | Primary paired-change hypothesis not testable/supported because all five DNMT3A sites are invariant |
| PTM conservation | 41 mouse sites PASS | 31 concordant mappings; 10 unresolved | PASS | Report PTM-capable residue conservation, not modification conservation |
| G9a-contact motif | Mouse L139/F140/I141 functional evidence PASS | Maps concordantly to human NP105/106/107 | 289-orthologue atlas PASS | Report lineage-structured chemistry, not conserved binding |
| Comparative structure | WAIT | WAIT | UNKNOWN | Admit only after a sequence contrast is nominated |

## Active independent computation

The existing replicated 5-ns FF-interface MD batch in `../therian_ff_simulation_20261104/production_hmr/` is independent of this add-on and must not be overwritten, restarted, or used before its completion validator passes.

## Completed outputs

1. Adult germline atlas audit: all nine H5AD files (716,925,258 bytes total) were downloaded, validated, and checksummed. All are processed, non-integer matrices and all sampled animals are adult or post-juvenile, so counts-based inference and an imprint-establishment timing test are not admitted. DNMT3L is testable in six species; marmoset, platypus, and chicken are unresolved rather than expression-negative.
2. Matched partner inventory: 265 species in `results/dnmt3a_dnmt3l_matched_interface_states_vertebrate.tsv`. DNMT3A1 E725/G726/R729/F732/F772 are invariant in all 265 resolved pairs; DNMT3L NP261/NP291/NP293 vary. The prespecified direct-contact compensation model is therefore not supported at these five partner residues.
3. PTM source audit: 17 acetylation and 24 phosphorylation sites recovered from Nam et al. 2026; processed PRIDE project PXD032953 files retained with checksums.
4. PTM coordinate atlas: all 17 acetylation sites and 14 phosphorylation sites map concordantly through MAFFT and MACSE; seven phosphorylation sites map to human gaps and three are aligner-discordant. Outputs are `results/dnmt3l_mouse_human_ptm_crosswalk.tsv`, `results/dnmt3l_ptm_conservation_by_site.tsv`, and `results/dnmt3l_ptm_conservation_by_clade.tsv`.
5. Functional PTM contrast: mouse K238 maps to human NP204 and is retained as lysine in 64/289 resolved orthologues; mouse K412 maps to human NP377 and is lysine in 231/288. These values describe capacity for lysine modification only.
6. G9a-contact motif: mouse L139/F140/I141 maps to human NP105/106/107 (LLI). The exact mouse LFI triplet occurs in 28/289 orthologues, including 27 sampled muroid rodents and one amphibian. This is a comparative interaction-surface result, not evidence that G9a binding is lost elsewhere.
7. Regulatory evidence audit: the stage-resolved Kojima et al. 2026 marmoset study independently shows later DNMT3L deployment during prolonged postnatal de novo methylation. The 2023 Eumuroida oocyte-promoter study independently shows lineage-restricted cis-regulatory deployment. Both are prior published evidence.

## Next valid actions

1. Do not use the adult atlas to claim regulatory absence or expression differences; it failed the developmental-stage and raw-count gates.
2. Add the G9a-contact motif result only with its functional-evidence and phylogenetic non-independence limits.
3. Treat promoter/isoform evolution and stage-resolved primate deployment as published mechanistic context unless a new transcript-resolved comparative dataset is analyzed.
4. Keep new comparative structure prediction on hold. The existing replicated FF-interface MD batch is independent and remains in progress; it must pass its completion validator before interpretation.

