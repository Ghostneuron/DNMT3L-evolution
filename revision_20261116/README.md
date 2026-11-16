# DNMT3L falsification-audit reproducibility update

This additive repository update records the public data, figures, and code associated with the 16 November 2026 falsification audit of the current working manuscript, **“Tetrapod DNMT3L locus evidence and ancestral interface comparisons.”** It does not replace the preserved v1.1.1 release. The manuscript, audit report, private correspondence, and submission materials are not distributed in this public repository.

## Main changes represented here

- NP349 is kept separate from the 18 whole-protein substitutions that passed both topology and alignment gates.
- A leave-one-amphibian sensitivity analysis preserves K as the maximum-posterior Amniota state and Q as the sampled-Theria state on both topologies. When *Geotrypetes seraphini*, the only K-bearing amphibian, is omitted, the Amniota K posterior falls to 0.518 on the taxonomy topology and 0.499 on the gene tree, below the prespecified 0.80 call threshold. The sampled-Theria Q posterior remains 0.968 and 0.977.
- The 18 robust Amniota-to-sampled-Theria substitutions and their interface-enrichment tests are provided as separate tables.
- Selection-model null results, locus evidence, non-mammalian transcriptional support, structural-prediction seed variation, and independent molecular-dynamics replica summaries are retained in machine-readable form.
- Molecular-dynamics replicas are reported as independent observations. Replica indices are not treated as matched pairs.
- The six current figures use the audited terminology, coordinate crosswalks, and uncertainty labels.

## Directory map

- `analysis/switching_helix_asr_337/`: alignment and two fixed topologies required for the NP349 leave-one-amphibian analysis
- `analysis/therian_ff_simulation_20261104/analysis_summary/`: final independent-replica MD summaries; trajectories are intentionally omitted because of size
- `analysis/phylo_balanced_structure_pilot_20261104/results/full_monomer_analysis/`: structure-prediction confidence and seed-dispersion summaries
- `analysis/regulatory_coevolution_ptm_20261104/results/`: matched DNMT3A–DNMT3L, G9a-motif, and PTM-coordinate summaries
- `results/`: locus, expression, alignment, selection, ancestral-state, and interface result tables used in the audit
- `scripts/`: NP349 sensitivity analysis and audited MD-figure builder
- `figures/`: current Figures 1–6 in PDF and PNG

## Reproducing the new sensitivity analysis

Requirements are Python 3 with Biopython and IQ-TREE 2 available as `iqtree2` on `PATH`, or supplied through `IQTREE2_BIN`. The completed run used IQ-TREE 2.3.6.

From this directory:

```bash
python scripts/run_np349_leave_one_amphibian_sensitivity_337.py
```

The script prunes each of four amphibian tips in turn, fixes both recorded topologies, runs amino-acid ancestral reconstruction under `LG+F+R4`, verifies exact clade membership, and rewrites `results/switching_helix_337_np349_leave_one_amphibian.tsv`.

A clean-copy four-thread smoke rerun reproduced every maximum-posterior state and every 0.80-threshold decision. The largest posterior difference was 0.00008, consistent with last-decimal floating-point variation in multithreaded likelihood optimization.

The MD figure can be rebuilt with Python, pandas, NumPy, and Matplotlib:

```bash
python scripts/build_falsification_audited_md_figure_337.py
```

## Scope and interpretation

The package supports comparative and hypothesis-generating conclusions. Ancestral reconstruction is conditional on the included sequences, alignments, models, and topologies. Structural predictions and short simulations do not measure biochemical function. The simulation summary retains the contact present in the common starting template and does not claim a replica-consistent endpoint shift for the other readouts.

## Exclusions

This public update excludes the manuscript, the falsification-audit PDF, journal correspondence, temporary renders, raw molecular-dynamics trajectories, model caches, local paths, credentials, and macOS AppleDouble files.

Code under `scripts/` is licensed under the repository MIT license. Original figures, result tables, and documentation are licensed under CC BY 4.0, subject to the repository’s third-party-material exclusions.
