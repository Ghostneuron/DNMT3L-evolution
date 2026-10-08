# Clade-aware evolution of vertebrate DNMT3L

This repository accompanies the manuscript **"Clade-aware analysis of vertebrate DNMT3L reveals modular constraint and lineage-structured switching-helix residue states."** It contains the data products, scripts, figures, and supplementary material for a clade-aware analysis of vertebrate DNMT3L sequence evolution. The manuscript itself is not distributed in this release.

## Study overview

The study evaluates DNMT3L at protein scale rather than centering inference on one residue. It combines orthology and paralog checks, independently filtered codon alignments, four taxonomic scopes, three alignment strategies, phylogenetic sensitivity analyses, codon-selection models, saturation diagnostics, and structure-informed residue mapping.

The four analyzed scopes are Placentalia, Euarchontoglires, Laurasiatheria, and non-avian Sauropsida. Placentalia overlaps the two restricted placental subclades and is treated as a pooled scope, not as an independent biological replicate.

The primary selection cohort contains 269 nonhybrid records from an archived set of 321 candidates. This cohort is conditional on sequence-quality, domain-coverage, and human-reference h358-mappability filters. The h358 gate excluded eight otherwise eligible records: five placentals and three amphibians. Protein-wide estimates therefore describe the h358-gated cohort; an ungated reconstruction was not performed.

## Main findings

- The ADD domain is more conserved than the C-terminal methyltransferase-like region in placental analyses by both site entropy and a tree-aware minimum-change sensitivity.
- The switching helix and adjacent loop are more variable than the remaining C-terminal sites in sitewise placental comparisons, but neither region is exceptional against all same-length C-terminal windows after multiple-testing correction.
- Human-reference position 349, corresponding to structural Q348 in the human 9MPP comparison, shows a lineage-structured residue-state contrast: glutamine predominates among sampled therians and lysine among sampled non-avian sauropsids. This is a biochemical hypothesis, not evidence of a measured nonhuman structure or function.
- Alignment-robust support identifies 24 human-mapped diversifying candidate sites at the exploratory tier; 10 scope-site calls meet the strict tier.
- Conditional multinucleotide models retain none of seven strict episodic candidates and all three strict pervasive candidates under the prespecified cross-aligner rule.
- Human-reference position 358 is retained as a prespecified complementary analysis. It is not the organizing candidate and does not show supported episodic or cross-clade recurrence.

## Supplementary tables

The formatted tables are in:

- `documents/DNMT3L_supplementary_tables_20261007.docx`

Machine-readable versions are in:

- `supplementary/Table_S1_h358_primary_results.csv`
- `supplementary/Table_S2_h358_topology_sensitivity.csv`
- `supplementary/Table_S3_h358_multinucleotide_sensitivity.csv`

## Release bundles

Two complementary bundles are prepared from the same validated analysis state:

| Bundle | Intended use | Contents |
|---|---|---|
| GitHub | Browsing, version control, figures, documents, principal derived data, and public scripts | Selected alignments, inputs, summaries, phylogenies, structural outputs, figures, supplementary tables, and scripts. Large raw model outputs and execution logs are omitted. |
| Zenodo | Immutable reproducibility archive | The GitHub materials plus complete `analysis_v2` model outputs, logs, deprecated-output provenance, and the archived source analysis required for a full reconstruction. |

The Zenodo bundle is the complete computational record. A clean from-source rerun should use that bundle. The GitHub bundle is deliberately smaller and contains the material most useful for review and reuse.

## Directory map

- `analysis_v2/inputs`: source tables and sequence inputs for the broadened analysis
- `analysis_v2/raw_alignments`: independently generated or archived unfiltered alignments
- `analysis_v2/filtered_alignments`: alignment-specific occupancy filtering and masks
- `analysis_v2/selection_inputs`: duplicate-collapsed codon alignments used for selection analysis
- `analysis_v2/phylogenies`: inferred codon trees and diagnostics
- `analysis_v2/species_trees`: NCBI-taxonomy topology sensitivities
- `analysis_v2/hyphy`: primary gene-tree HyPhy JSON outputs (Zenodo bundle)
- `analysis_v2/species_tree_hyphy`: species-tree HyPhy JSON outputs (Zenodo bundle)
- `analysis_v2/fubar_isolated_inputs`: topology-isolated FUBAR inputs and cache record (Zenodo bundle)
- `analysis_v2/conditional_sensitivity`: conditional error-sink and multinucleotide models (Zenodo bundle)
- `analysis_v2/saturation_diagnostics`: pairwise codon diagnostics
- `analysis_v2/orthology_sensitivity`: DNMT3-family sensitivity analyses
- `analysis_v2/structure`: structural coordinate and residue-state outputs
- `analysis_v2/summaries`: machine-readable final summaries
- `figures_v2`: publication figures in PNG and PDF
- `documents`: formatted supplementary tables
- `supplementary`: machine-readable supplementary tables
- `scripts`: analysis and figure-generation scripts
- `environment`: recorded software and Python-package versions
- `provenance`: dated methodological decision records and validation evidence
- `source_package`: archived source analysis used by the broadened workflow (complete in Zenodo; minimal dependencies in GitHub)

## Quick validation

From the package root:

```bash
python scripts/validate_public_package.py
```

The validator checks required files, the package manifest, every SHA256 digest, unwanted filesystem metadata, and the GitHub 100 MB per-file limit.

Install the recorded Python dependencies in an isolated environment if you plan to rerun Python scripts:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r environment/requirements.txt
```

See `COMMANDS.md` for workflow order and external software requirements.

## Reproducibility notes

MACSE, MAFFT, PRANK, IQ-TREE, HyPhy, and PyMOL are external requirements and are not distributed in the archive. Exact versions used in the reported analysis are listed in `environment/software_versions.txt`.

Public script copies use environment-variable overrides for executable discovery and otherwise expect executable names on `PATH`: `MACSE_BIN`, `MAFFT_BIN`, `PRANK_BIN`, `IQTREE2_BIN`, and `HYPHY_BIN`. Scientific parameters are unchanged from the validated analysis. See `provenance/CODE_CHANGES_FOR_PUBLIC_RELEASE.md`.

The dated records under `provenance/` preserve the timing and rationale of sensitivity analyses. Absolute paths inside historical logs or records describe the original execution environment and are not required package locations.

## Interpretation boundaries

These data support comparative sequence and structural-context conclusions. They do not directly measure expression, methylation, protein assembly, nucleosome binding, enzyme activity, development, disease, or organismal phenotype. The human 9MPP structure supplies positional context; it does not establish the conformation or biochemical effect of a residue in a nonhuman protein.

## Citation and identifiers

Repository: <https://github.com/Ghostneuron/DNMT3L-evolution>

Citation metadata are provided in `CITATION.cff`. The Zenodo DOI remains pending and should be added after the deposit is created. Until then, cite the title, author, release version, and release date shown there.

## License status

No public-use license has yet been selected by the author. See `LICENSE_SELECTION_REQUIRED.md` before publishing the repository or Zenodo record. A license must be selected explicitly; the presence of files in these bundles does not itself grant reuse rights.
