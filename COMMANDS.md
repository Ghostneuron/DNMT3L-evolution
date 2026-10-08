# Reproduction commands

Run commands from the root of an unpacked release bundle. The complete from-source workflow requires the Zenodo bundle; the compact GitHub bundle omits large raw model outputs and most archived source data.

## Environment

Create an isolated Python environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r environment/requirements.txt
```

Install the external executables listed in `environment/software_versions.txt`. Public script copies first use the following optional environment variables and otherwise use executable names on `PATH`:

```bash
export MACSE_BIN=/path/to/macse
export MAFFT_BIN=/path/to/mafft
export PRANK_BIN=/path/to/prank
export IQTREE2_BIN=/path/to/iqtree2
export HYPHY_BIN=/path/to/hyphy
```

PyMOL is invoked directly from the shell and must be on `PATH`.

## Validate the downloaded package

```bash
python scripts/validate_public_package.py
```

## Build inputs, alignments, and trees

The following order reproduces the broadened inputs and clade-specific alignments from the archived source analysis in the Zenodo bundle:

```bash
python scripts/prepare_broadened_inputs.py
python scripts/run_broadened_alignments.py
python scripts/filter_independent_alignments.py
python scripts/run_dnmt3c_orthology_sensitivity.py
python scripts/run_codon_trees.py
python scripts/build_ncbi_species_trees.py
python scripts/summarize_phylogeny_diagnostics.py
python scripts/build_reference_coordinate_crosswalk.py
```

`run_broadened_alignments.py` preserves the validated archived MACSE products for all four scopes and the archived placental MAFFT and PRANK products, then generates the restricted-scope MAFFT and PRANK alignments. `analysis_v2/raw_alignments/alignment_manifest.csv` records provenance, commands, and output hashes.

## Constraint and structural summaries

```bash
python scripts/summarize_constraint_landscape.py
python scripts/assess_codon_saturation.py
python scripts/assess_structural_region_window_sensitivity.py
python scripts/assess_tree_aware_constraint.py
python scripts/summarize_structural_state_survey.py
pymol -cq scripts/render_switching_helix.pml
```

## Primary selection models

```bash
python scripts/run_hyphy_selection.py
python scripts/run_species_tree_hyphy.py
python scripts/rerun_fubar_topology_isolated.py
```

Do not summarize the first-pass FUBAR results before running `rerun_fubar_topology_isolated.py`. That step isolates caches by topology and tree hash, replaces the active FUBAR JSON files with isolated-cache results, and records the operation in `analysis_v2/fubar_isolated_inputs/fubar_topology_isolation_manifest.csv`.

## Selection summaries and conditional sensitivities

```bash
python scripts/summarize_broadened_selection.py
python scripts/summarize_species_tree_sensitivity.py
python scripts/run_conditional_sensitivities.py
```

Conditional analyses are determined by the primary summaries. They qualify primary results but cannot create a primary finding. Details of correction families, focal-site filters, and retention rules are preserved in the manuscript Methods and dated records under `provenance/`.

## Figures

```bash
python scripts/build_dnmt3l_broadened_figures.py
```

The formatted supplementary-table document and its three machine-readable CSV tables are supplied as release artifacts. The manuscript files and manuscript-generation scripts are intentionally excluded from the public bundles. Their exclusion does not affect reproduction of any scientific calculation, summary, or figure.

## Computational cost

Alignment, phylogeny, and HyPhy stages can take several hours. The package already contains the reported outputs, so rerunning those stages is unnecessary for reading, validation, or secondary inspection.
