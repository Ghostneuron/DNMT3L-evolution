# Transition-informed selection decision rules

Frozen before inspecting transition-informed HyPhy outputs on 2026-11-02.

## Inputs

- Alignment: 289 validated DNMT3L protein/CDS pairs projected onto the 387 human-coordinate codons.
- NP349, NP351, NP352, and NP355 calls that disagree among available global-MAFFT, clade-MACSE, and pairwise mappings are masked as `NNN`.
- Primary topology: rooted NCBI taxonomy. Topology sensitivity: validated DNMT3L gene tree.

## Foreground definitions

- Therian-stem foreground: the single branch connecting Amniota to Theria. This branch was nominated by the NP349 reconstruction and is therefore exploratory, not an independently prespecified test.
- RELAX test set: all sampled therian branches. Reference set: all sampled sauropsid branches. Amphibian branches are unlabelled background.

## Tests and interpretation

- aBSREL tests episodic diversifying selection somewhere in the gene on the therian stem. It does not identify NP349 as causal.
- BUSTED tests gene-wide evidence that at least one site experienced episodic diversifying selection on the therian stem. It does not establish selection at NP349.
- RELAX tests whether the distribution of selection intensity differs between sampled therian and sauropsid branches. It does not test the direction or function of the NP349 substitution.
- Primary claims use the taxonomy topology. The gene tree is a sensitivity; disagreement is reported as topology sensitivity rather than combined evidence.
- P values are reported exactly. These three model families answer different questions and are not pooled into one omnibus claim. Within-test corrections performed by HyPhy are retained.
- A branch-wide positive result cannot be transferred to NP349. A negative result does not show that the residue change was neutral or that selection never occurred.

## Stop rule

- Do not run additional branch/site combinations after inspecting these outputs unless they are explicitly labeled post hoc and placed in a new testing family.
