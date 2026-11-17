# Expression-support decision rules

Frozen before inspecting processed expression values on 2026-11-03.

## Estimand

For each non-mammalian locus already supported by the locus-history analysis, test whether a sequence- or identifier-compatible DNMT3L transcript has nonzero processed RNA-seq abundance in a named tissue and developmental stage.

## Biological unit

The biological unit is the independent animal or sequencing library represented by one processed-data column. Technical replicates are not counted as independent biological support.

## Admission gates

1. The species must have locus grade A, B1, or B2 in `results/vertebrate_dnmt3l_locus_matrix.tsv`.
2. The processed dataset must expose sample-level abundance or counts and recoverable tissue/stage metadata.
3. The expression row must map to the supported locus by an exact stable identifier or by a sequence comparison that distinguishes DNMT3L from DNMT3A, DNMT3B, and DNMT3C.
4. A dataset with only an ambiguous family-level match is reported as unresolved, not positive.
5. Expression is positive only when at least one admitted biological library has a processed value greater than zero. The number of positive and tested libraries must both be reported.

## Interpretation boundaries

- Positive expression supports transcription of the locus only in the sampled tissue/stage.
- Absence of processed reads is not evidence of genomic absence or universal transcriptional silence.
- Expression magnitudes are not compared across species, assemblies, quantification pipelines, tissues, or developmental stages.
- Expression does not establish biochemical activity, functional equivalence, residue effects, or reproductive function.
- Single-library evidence remains descriptive.
- Raw-data reprocessing is deferred unless processed artifacts fail for a supported locus and a bounded smoke test is justified.

## Frozen output states

- `SUPPORTED`: locus-compatible row and at least one nonzero biological library.
- `NOT_DETECTED`: locus-compatible row found, but all admitted libraries are zero.
- `IDENTIFIER_UNRESOLVED`: processed artifact exists but its rows cannot be linked specifically to the supported DNMT3L locus.
- `GROUP_SUMMARY_ONLY`: a locus-compatible processed row is available only as a pooled or group-level summary; it is retained as secondary descriptive evidence but does not pass the sample-level primary gate.
- `NO_PROCESSED_ARTIFACT`: no suitable processed matrix is available.
- `NOT_ELIGIBLE_LOCUS`: locus grade is C or D, so expression is not used to rescue locus identity.

`GROUP_SUMMARY_ONLY` was added after inspecting the artifact layout and before scripted summarization. It preserves, rather than relaxes, the prespecified sample-level admission gate.

