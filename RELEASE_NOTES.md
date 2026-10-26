# Release notes

## Version 1.1.0 - 2026-10-25

Additive revision release for **"Modular evolutionary constraint and lineage-structured switching-helix variation in DNMT3L."** The v1.0.0 analysis state remains preserved by its tag and commit.

### Added

- Ancestral-state reconstruction at human-reference positions NP349 and NP352 under taxonomy-step and unit branch lengths
- An h358-ungated placental sensitivity that restores five otherwise eligible records (213 to 218 sequences)
- Genomic-neighborhood checks for seven representative non-avian sauropsid DNMT3L loci
- A clean-start HyPhy 2.5.103 GARD analysis and partition-aware FUBAR, FEL, and MEME sensitivity for the Sauropsida MACSE alignment
- Portable revision inputs, outputs, scripts, diagnostics, and checksums under `revision_20261025/`
- Updated formatted supplementary tables and consolidated machine-readable workbook

### Updated

- Revised Figures 3 and 4 in PDF and PNG
- Repository title, findings summary, citation metadata, package metadata, commands, manifest, and checksums
- Public interpretation boundaries to keep structure and selection statements at the level supported by comparative sequence evidence

### Excluded

- Manuscript and cover-letter files
- Journal decision letters, reviewer or editor correspondence, and submission-platform records
- Local execution paths, credentials, caches, temporary files, and AppleDouble sidecars

## Version 1.0.0 - 2026-10-07

Initial public-release candidate for the broadened DNMT3L evolutionary analysis.

### Included

- Formatted supplementary tables and three machine-readable h358 tables
- Four main figures and two supplementary figures in PNG and PDF
- Three-aligner data products for four taxonomic scopes
- Orthology, phylogeny, constraint, saturation, structural-state, and selection summaries
- Analysis and figure-generation scripts
- Recorded software environment, manifests, SHA256 checksums, and public-package validator
- Dated methodological and validation provenance
- Complete raw model outputs and archived source analysis in the Zenodo bundle
- Path-specific open licensing: MIT for original code and CC BY 4.0 for other original content

### Public-release cleanup

- Excluded macOS AppleDouble files, `.DS_Store`, caches, temporary files, and Codex staging directories
- Excluded the current manuscript in Word and PDF formats
- Excluded the private journal-selection and resubmission-strategy document
- Excluded superseded manuscript and response documents from the archived source package
- Changed executable discovery in public script copies from machine-specific absolute paths to environment variables with `PATH` fallbacks
- Excluded current and archived manuscript, response-letter, and journal-strategy generation scripts because they contain nonpublic document text

No scientific parameter, result, figure, table, manuscript claim, or archived numerical output was changed during packaging.

### Remaining archive step

- Mint and add the Zenodo DOI
