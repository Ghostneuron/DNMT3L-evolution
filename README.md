# DNMT3L evolutionary analysis: reproducible data and code

This repository contains reproducible data, scripts, figures, and supplementary material for comparative analyses of DNMT3L evolution. The current working manuscript is **"Tetrapod DNMT3L locus evidence and ancestral interface comparisons."** The manuscript itself is not distributed in this public repository.

The additive **16 November 2026 falsification-audit update** is under `revision_20261116/`. It supplies the current six-figure set, audited summary tables, and the reproducible NP349 leave-one-amphibian sensitivity analysis. The historical formal release **v1.1.1 (25 October 2026)** remains preserved under `revision_20261025/` and in the original directories.

## Current audit update

- The 18 alignment- and topology-robust Amniota-to-sampled-Theria substitutions are separated from NP349, which failed the whole-protein alignment gate.
- Leave-one-amphibian analyses preserve the maximum-posterior NP349 endpoint states, but the Amniota K call falls below the 0.80 threshold when the only K-bearing amphibian is omitted.
- Selection-model null results and nonsignificant interface enrichment are reported alongside positive findings.
- Structural contrasts are bounded by independent-seed variation, and molecular-dynamics replicas are treated as independent rather than artificially paired.

See `revision_20261116/README.md` for exact contents, commands, and interpretation limits.

## Study overview

The study evaluates whether DNMT3L functional modules differ in evolutionary constraint and whether switching-helix residue states are structured by lineage. It combines orthology and paralog checks, independently filtered codon alignments, four taxonomic scopes, three alignment strategies, phylogenetic sensitivity analyses, codon-selection models, saturation diagnostics, structure-informed residue mapping, ancestral-state reconstruction, and representative genomic-neighborhood checks.

The four analyzed scopes are Placentalia, Euarchontoglires, Laurasiatheria, and non-avian Sauropsida. Placentalia overlaps the two restricted placental subclades and is treated as a pooled scope, not as an independent biological replicate.

The primary selection cohort contains 269 nonhybrid records from an archived set of 321 candidates. This cohort is conditional on sequence-quality, domain-coverage, and human-reference h358-mappability filters. The h358 gate excluded eight otherwise eligible records: five placentals and three amphibians. Release v1.1.0 restores the five placentals in an explicit ungated sensitivity analysis, expanding that cohort from 213 to 218 sequences.

## Main findings

- The ADD domain is more conserved than the C-terminal methyltransferase-like region in placental analyses by both site entropy and a tree-aware minimum-change sensitivity; the contrast persists after restoring the five h358-gate exclusions.
- The switching helix and adjacent loop are more variable than the remaining C-terminal sites in sitewise placental comparisons, but neither region is exceptional against all same-length C-terminal windows after multiple-testing correction.
- Human-reference position 349, corresponding to structural Q348 in the human 9MPP comparison, shows a lineage-structured residue-state contrast. The inferred therian ancestor is Q (posterior probability 0.9990) and the inferred sampled sauropsid ancestor is K (0.9969) under taxonomy-step branch lengths; both exceed 0.994 under unit branch lengths.
- At human-reference position 352, the inferred therian ancestral Q state is stable, whereas the sampled sauropsid reconstruction is branch-length sensitive.
- All seven representative sauropsid loci examined retain at least two named markers from the human DNMT3L neighborhood, supporting local locus continuity without proving function.
- A Sauropsida MACSE GARD analysis identifies two nucleotide breakpoints (165 and 894). Partition-aware site models detect no positively selected FUBAR sites and no FEL or MEME diversifying sites at q <= 0.10, while 115 sites meet the joint purifying criterion at q <= 0.05.
- Human-reference position 358 is retained as a prespecified complementary analysis. It is not the organizing candidate and does not show supported episodic or cross-clade recurrence.

All residue-state and selection results are comparative sequence inferences. They do not establish altered DNMT3L structure, nucleosome binding, methylation, or organismal phenotype.

## Supplementary tables

The current formatted tables are in:

- `documents/DNMT3L_supplementary_tables_20261025.docx`

The consolidated machine-readable workbook is:

- `supplementary/DNMT3L_machine_readable_tables_S1-S14_20261025.xlsx`

The three original v1.0.0 machine-readable h358 tables remain in:

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
- `revision_20261025`: additive v1.1.0 scripts, inputs, results, diagnostics, checksums, and a revision-specific README
- `revision_20261116`: additive falsification-audit figures, summary results, NP349 sensitivity inputs/code, and checksums

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

Citation metadata are provided in `CITATION.cff`. The Zenodo DOI remains pending and should be added after the deposit is created. Until then, cite the title, author, release version, and release date shown there. Release v1.0.0 is preserved by tag and commit `b6aa805f32471c0634519ae35b4fd2a94f3728a3`.

## License

This repository uses path-specific dual licensing:

- Original source code in `scripts/` directories is licensed under the [MIT License](LICENSE-CODE).
- Other original data, figures, supplementary material, documentation, metadata, and provenance records are licensed under [Creative Commons Attribution 4.0 International](LICENSE-DATA).

See [LICENSE](LICENSE) for scope, attribution, and third-party-material exclusions. The licenses apply only to rights held by Jie Lu; source database and structural records remain subject to their original terms.
