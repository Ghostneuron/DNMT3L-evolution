# DNMT3L evolution data and code

This repository contains the current, browsable data-and-code release supporting the manuscript **A muroid-specific DNMT3L G9a-association triplet among eutherians and conserved DNMT3A contacts across tetrapods**.

Published preprint: https://doi.org/10.5281/zenodo.24299393 (version 1.0.1)

Stable preprint concept DOI: https://doi.org/10.5281/zenodo.24299392

Published Zenodo dataset: https://doi.org/10.5281/zenodo.24298873 (version 1.0.1)

Stable dataset concept DOI: https://doi.org/10.5281/zenodo.20574053

## Main result and scope

All 27 sampled muroids with accepted DNMT3L NP105/NP106/NP107 mappings carry LFI, whereas none of the other 192 accepted eutherian mappings do. The amino-acid state is fixed in the sampled muroids, but the phenylalanine codon is not: 22 use TTC and five use TTT. The TTT carriers are not one sampled clade, and the sampled NCBI Taxonomy topology requires at least two synonymous codon changes within Muroidea.

The predicted G9a partner state at residues 699-703 is KLRFH in all 212 eutherians with a mapped EHMT2 orthologue. This partner-side result uses an AlphaFold-Multimer interface model (ModelArchive ma-dm-mmk9-101; rank-1 mean pLDDT 60.4, pTM 0.383, ipTM 0.533) and does not validate the interface. The DNMT3L triplet is supported by mutagenesis; its placement at the interface is model-derived.

## Contents

- `manuscript/`: preprint-ready manuscript and supplement in DOCX and PDF.
- `results/`: canonical machine-readable result tables, including the 219-eutherian triplet-codon analysis, the 212-species G9a partner analysis, matched DNMT3A-DNMT3L interface results, locus evidence, alignment tests, ancestral reconstructions, and selection summaries.
- `figures/`: the five main figures and two supplementary figures in PDF and PNG.
- `data/processed/`: compact analysis-ready sequence sets and search inputs.
- `analysis/`: reproducibility inputs and selected outputs for alignment, trees, ancestral reconstruction, expression support, locus history, structure prediction, and molecular-dynamics summaries.
- `scripts/`: reproducible analysis and figure scripts. Manuscript-revision utilities and reviewer-specific files are excluded.
- `config/` and `metadata/`: decision rules, commands, software versions, data manifests, and checksums.
- `MANIFEST.tsv` and `SHA256SUMS.txt`: file inventory and integrity checksums.

## GitHub and Zenodo scope

GitHub contains the manuscript-facing files, code, processed inputs, principal analysis outputs, figures, and machine-readable results. The immutable Zenodo record additionally contains the versioned third-party raw-data packages, the two large axolotl expression-search files, and full ancestral-state `.state` outputs. These bulky files are intentionally not duplicated here; their absence does not change the scientific results or the parameters recorded in this repository.

Private correspondence, review documents, manuscript-revision utilities, temporary files, caches, and machine-specific paths are excluded from both public release surfaces.

## Molecular-dynamics and structure-prediction scope

The repository includes the per-replica molecular-dynamics analysis tables, completion records, summary tables, scripts, and environment file. Raw trajectories, checkpoints, final-state XML files, prediction caches, and redundant prediction coordinates are excluded because the manuscript treats these analyses as supplementary and reports the included summary outputs. The structure-prediction release subset contains the recorded panel, confidence summaries, and analysis scripts.

## Reproducibility notes

Analyses use project-relative paths. Software versions are recorded in `metadata/`, analysis-specific environment files, and run manifests. Source accessions and provenance are recorded in `metadata/dataset_manifest.tsv` and `metadata/source_artifact_manifest.tsv`. Third-party source records remain subject to their original terms; the repository licenses apply only to rights held by the author.

`MANIFEST.tsv` and `SHA256SUMS.txt` describe the GitHub tree and therefore differ from the corresponding inventory files in the complete Zenodo archive.

## Citation

Paper:

Lu J. (2026). A muroid-specific DNMT3L G9a-association triplet among eutherians and conserved DNMT3A contacts across tetrapods. Preprint, version 1.0.1. Zenodo. https://doi.org/10.5281/zenodo.24299393

Data and code:

Lu J. (2026). A muroid-specific DNMT3L G9a-association triplet among eutherians and conserved DNMT3A contacts across tetrapods: data and code. Version 1.0.1. Zenodo. https://doi.org/10.5281/zenodo.24298873
