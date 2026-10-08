# Codon tree diagnostic audit

Date: 6 October 2026

## Purpose

This post hoc quality audit was added while the broadened HyPhy jobs were still running and before the final selection summaries were generated. It does not create or remove candidate sites. It records limitations of the inferred codon trees used by the primary HyPhy analyses.

## Checks

`scripts/summarize_phylogeny_diagnostics.py` parses all 12 IQ-TREE reports and records sequence count, codon-site count, free-parameter count, the parameter-to-site ratio, near-zero internal branches, internal-branch contribution to total tree length, and the distribution of internal SH-aLRT support. Alignment and tree SHA256 values link every diagnostic row to its inputs.

The three pooled Placentalia trees contain more estimated parameters than alignment sites because approximately 210 terminal sequences require a branch-rich tree from 373-382 retained codons. IQ-TREE explicitly warns that those parameter estimates require caution. The nonoverlapping placental subclades and Sauropsida remain below that threshold, although every inferred tree contains some near-zero internal branches.

## Interpretation

This limitation reinforces three existing rules:

1. pooled Placentalia is not treated as an independent biological replicate of its subclades;
2. a site-level claim requires concordance across independent alignments and models; and
3. every MACSE analysis is repeated on an independently constructed NCBI-taxonomy topology, with topology sensitivity allowed to weaken but not create a primary candidate.

The diagnostic does not prove that either inferred codon trees or taxonomy topologies are correct. It makes the tree uncertainty visible and keeps strongly significant gene-wide tests from being assigned automatically to a particular codon or mechanism.

