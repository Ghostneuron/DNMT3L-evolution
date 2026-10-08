# FUBAR topology-cache audit

Date identified: 2026-10-06, before selection summaries or manuscript interpretation.

## Issue

HyPhy FUBAR writes an alignment-adjacent `*.FUBAR.cache` file. The initial species-tree sensitivity and gene-tree primary run used the same MACSE alignment path. The first completed gene-tree FUBAR output was numerically identical to the earlier species-tree output despite nonidentical tree files, showing that the second run had reused topology-dependent cached likelihoods.

## Correction

- Do not interpret any FUBAR output created from the shared selection-input paths.
- Preserve the initial outputs in `analysis_v2/deprecated/fubar_shared_cache_20261006/` for provenance.
- Copy every alignment into a topology-specific directory before FUBAR execution.
- Use separate copies and therefore separate cache files for the 12 gene-tree analyses and four species-tree sensitivities.
- Overwrite the primary FUBAR JSON destinations only with results generated from those isolated inputs.
- Record alignment, tree, cache, output, command, elapsed time, and SHA256 values in `analysis_v2/fubar_isolated_inputs/fubar_topology_isolation_manifest.csv`.
- Refresh the FUBAR rows in the gene-tree and species-tree parent manifests with the final isolated commands and output hashes. Preserve each superseded shared-path command in a separate manifest field.
- Move every alignment-adjacent cache from the active `selection_inputs` directory into `analysis_v2/deprecated/fubar_shared_cache_20261006/selection_input_caches/`, record its original and archived hash, and leave no shared-path FUBAR cache beside an active selection input.

FEL, MEME, and BUSTED do not use this FUBAR cache and are unaffected.
