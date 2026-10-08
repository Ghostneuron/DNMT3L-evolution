# Release notes

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
