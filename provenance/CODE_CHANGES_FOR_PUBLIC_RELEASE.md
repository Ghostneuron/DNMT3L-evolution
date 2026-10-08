# Code changes for the public release

The public script copies differ from the validated working-copy scripts in one packaging-only way.

1. Machine-specific executable paths were replaced by environment-variable lookups with executable-name fallbacks on `PATH`:
   - `MACSE_BIN` or `macse`
   - `MAFFT_BIN` or `mafft`
   - `PRANK_BIN` or `prank`
   - `IQTREE2_BIN` or `iqtree2`
   - `HYPHY_BIN` or `hyphy`

The current manuscript builder and archived manuscript, response-letter, and journal-strategy builders are excluded entirely because they contain nonpublic document text.

No model option, random seed, filtering threshold, taxonomic scope, correction family, result-processing rule, figure input, manuscript text, or numerical output was changed. The original commands and software versions remain recorded in historical manifests and logs in the Zenodo bundle.
