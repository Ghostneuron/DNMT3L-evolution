# DNMT3L revision-specific reproducibility files

This archive contains the scripts, inputs, and outputs added during the 25 October 2026 revision: ancestral-state reconstruction at NP349/NP352, the h358-ungated placental constraint sensitivity, representative sauropsid genomic-neighborhood checks, and the Sauropsida MACSE GARD/partition-aware selection sensitivity.

Within the public repository, these files constitute the additive v1.1.1 update for **"Modular evolutionary constraint and lineage-structured switching-helix variation in DNMT3L."** Version 1.1.1 supersedes v1.1.0 by restoring three checksum-listed HyPhy logs omitted by the repository ignore rule. Repository: <https://github.com/Ghostneuron/DNMT3L-evolution>.

The foundational validated analysis package remains available in DNMT3L-evolution release v1.0.0 (commit b6aa805f32471c0634519ae35b4fd2a94f3728a3). The scripts retained here document the exact revision workflow, but manuscript-building scripts also expect the frozen prior-submission DOCX templates and foundational v1.0.0 summary tables. Execution-host path strings in public text copies were normalized to portable relative labels; numerical and model content is unchanged. Manifests record both source and release hashes where normalization changes a result file.

Key software used in the revision was HyPhy 2.5.103 for the clean-start GARD analysis used for inference and the partition-aware selection reruns, MAFFT 7.526, Python 3.12/3.13, pandas, NumPy, SciPy, Biopython, Matplotlib 3.10.6, python-docx, and openpyxl. See the per-analysis manifests and script command-line arguments for model settings. The GARD manifest and diagnostics directory also record the recoverable HyPhy 2.5.94 checkpoint and the mixed-version resume that was excluded from inference after its AICc provenance was audited.

No decision letters, private correspondence, credentials, API caches, or AppleDouble sidecars are included. HyPhy's static example path in the FUBAR method-description text was also normalized in the public copy; the original working output remains unchanged.
