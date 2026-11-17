# DNMT3L representative protein-CDS selection rule

The 2026-11-02 NCBI Datasets archive contains 337 DNMT3L ortholog-group GeneIDs but 526 protein records and 526 CDS records. A single representative is selected per GeneID before phylogenetic analysis.

1. A protein and CDS are eligible to pair only when translation of the complete CDS, after removal of one terminal stop, exactly equals the protein sequence.
2. Candidate pairs are ranked prospectively by: sequence-QC pass; domain-QC pass; human-reference coverage (full length, C-terminal region, then ADD region); alignment score; aligned identity; curated RefSeq protein and transcript status; start methionine; proximity to the 387-aa human reference length; and accession as a deterministic final tie-breaker.
3. Sequence QC requires an in-frame CDS, no internal stop, no ambiguous codon, no protein X, and exact protein-CDS translation.
4. Domain QC uses NP_037501.2 coordinates and requires at least 0.80 full-reference coverage, at least 0.90 C-terminal-region coverage (human positions 178-379), and a mappable human position 358. ADD completeness is reported separately at at least 0.80 coverage of human positions 41-173.
5. Hybrid taxa are retained in the audit but excluded from alignment-ready sets.
6. NCBI ortholog-group membership is a candidate gate, not the final paralog-aware validation. The selected proteins are subsequently tested in a DNMT3A/DNMT3B/DNMT3L family tree before codon alignments are frozen.

These rules are applied uniformly to both retained and newly added GeneIDs; old accessions are not preferentially preserved.
