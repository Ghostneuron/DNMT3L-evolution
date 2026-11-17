# Germline regulatory evidence audit

Updated: 2026-11-04

## Decision

The available evidence supports **lineage- and stage-specific deployment of DNMT3L**, but it does not support a new causal claim that a particular DNMT3L protein substitution initiated genomic imprinting.

| Evidence source | Design and stage | What it contributes | Admissible use | Main limitation |
|---|---|---|---|---|
| Murat et al. cross-species spermatogenesis atlas (adult testis) | Single-nucleus RNA-seq from nine vertebrates; all sampled animals are adult or post-juvenile | Release-matched annotation audit and descriptive adult cell-state detection | Secondary evidence that adult DNMT3L measurement is sparse and annotation-sensitive | The downloaded matrices are processed, non-integer values rather than raw counts; adult sampling misses the principal imprint-establishment window; DNMT3L is unresolved in the marmoset H5AD and not testable in platypus or chicken |
| Kojima et al. 2026 primate germ-cell study | Marmoset scRNA-seq/scBS-seq from 103 cells across ten stages (1 day to 12 months), with 100 cells retained after QC; immunofluorescence across postnatal development | Stage-matched evidence that DNMT3L deployment accompanies later phases of prolonged postnatal de novo methylation | Primary published context for regulatory timing: DNMT3A was detected across methylation groups, whereas DNMT3L and DNMT3B were mainly detected in later groups; DNMT3L protein was absent in newborn testis and detected in a subset of later germ cells | This pattern is already published; it does not compare ancestral DNMT3L proteins or establish that DNMT3L caused the origin of imprinting |
| Behluli et al. 2023 oocyte-promoter study | Comparative oocyte RNA-seq and promoter sequence analysis | Establishes a lineage-restricted full-length oocyte DNMT3L promoter in the third intron of *Aire*, associated with the Eumuroida lineage and a TBPL2/TFIIA-compatible motif | Prior evidence that cis-regulatory evolution can change germline use without a pan-mammalian protein switch | This is not a new result of the present study; gene-level adult testis H5AD files cannot test promoter or isoform usage |
| Tatarakis et al. 2025 heterochromatin study | AlphaFold-Multimer nomination followed by recombinant binding, co-IP, localization, and mutagenesis | Maps a tested G9a-contact surface to mouse DNMT3L ADD residues L139/F140/I141 | Motivates the coordinate-qualified motif-conservation analysis below | Interaction was tested in mouse/ES-cell contexts; residue retention alone does not demonstrate binding in another species |

## Adult atlas gate results

- All nine H5AD files were downloaded, HDF5-validated, and checksummed (716,925,258 bytes total).
- All nine matrices are `float32` processed matrices. In a sample of 100,000 nonzero entries per file, 99.992%–99.999% were non-integers. Counts-based differential expression and pseudobulk inference are therefore not admitted.
- DNMT3A has a validated H5AD feature in all nine species. DNMT3L has a validated release-matched feature in human, chimpanzee, gorilla, macaque, mouse, and opossum.
- The marmoset DNMT3L identifier from the current annotation is absent from the older H5AD feature set. Platypus and chicken lack a validated DNMT3L feature. These three species are **not testable**, not expression-negative.
- Across the six testable species, adult DNMT3L detection is sparse (99.56%–99.98% of nuclei have zero stored value, depending on species). Because sampling is adult, this cannot be interpreted as absence during de novo methylation or imprint establishment.

## Experimentally anchored G9a-contact motif

- Mouse L139/F140/I141 maps concordantly in MAFFT and MACSE to human NP_037501.2 positions 105/106/107, where the human state is LLI.
- All 289 atlas orthologues resolve the triplet. The exact mouse LFI triplet occurs in 28/289: 27 sampled muroid rodents and one amphibian (*Geotrypetes seraphini*).
- The most frequent state is LLI (185/289), including 121/126 Laurasiatheria and the human sequence. Sauropsids are dominated by LIM (30/53), whereas marsupials are dominated by LLT (10/12).
- The lineage structure shows that the tested mouse contact chemistry is not universally conserved. It prioritizes direct comparison of mouse LFI and human/ancestral LLI, but it does **not** show that the interaction is absent, weaker, or imprint-specific in LLI species.
- The concentration of LFI in muroid rodents is notable beside the published lineage-restricted oocyte promoter, but the two observations are phylogenetically non-independent. Treat this as a candidate coupled regulatory/protein change, not as statistical evidence of coadaptation.

## Consequence for the manuscript

1. Regulatory timing and promoter choice provide a credible explanation for why no universal DNMT3L protein “imprinting switch” was recovered.
2. The adult cross-species atlas should not be presented as a new expression result in the main manuscript; it is best retained as an audit or supplement.
3. The G9a-contact motif atlas is an original, experimentally anchored comparative result and may enter a revised manuscript if framed as lineage-structured interaction-surface chemistry.
4. A mechanistic causal test still requires stage-matched promoter/isoform data in additional species or biochemical comparison of reconstructed proteins/complexes.

## Primary sources

- Kojima et al. 2026, *Nature Communications*, DOI: 10.1038/s41467-026-71405-z; public accessions DRA016143 and DRA016233; controlled human data JGAS000887 and JGAS000888.
- Behluli et al. 2023, *Epigenetics & Chromatin*, DOI: 10.1186/s13072-023-00518-2; PMID: 37924163.
- Tatarakis et al. 2025, *Molecular Cell*, DOI: 10.1016/j.molcel.2025.08.025; PMID: 40972525.
- Murat et al. 2023 cross-species spermatogenesis atlas, Bgee/CZI processed H5AD distribution.

