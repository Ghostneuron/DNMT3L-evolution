# Switching-helix ancestral-history decision rules

Frozen before inspecting ancestral-state outputs on 2026-11-02.

## Estimands

- Primary site: human NP_037501.2 position 349, with Q, K, and all other amino acids retained as distinct observed states in the amino-acid reconstruction.
- Secondary site: human NP_037501.2 position 352. Its greater observed heterogeneity is analyzed separately and is not pooled with NP349.
- Supporting interval: NP341-355, used to identify neighboring state changes without assigning a functional effect.

## Inputs and topology

- Identity and orthology are inherited from the frozen 289-record selection-ready set.
- The protein alignment is the existing MAFFT DNMT3-family alignment projected onto the 387 human NP_037501.2 residue columns and restricted to the 289 validated DNMT3L records.
- At NP349, NP351, NP352, and NP355, a tip state is retained only when all available global-MAFFT, clade-MACSE, and pairwise-atlas mappings agree. A missing or discordant focal call is masked as unknown. The unmasked MAFFT projection is an alignment sensitivity.
- The primary topology is the NCBI-taxonomy topology for all 289 tips. IQ-TREE re-estimates branch lengths from the full projected protein alignment under LG+F+R4.
- The validated DNMT3L gene-tree subtree is the topology sensitivity. Agreement between topologies is sensitivity, not independent replication.

## Ancestral-state calls

- IQ-TREE empirical-Bayes marginal probabilities are reported for all amino-acid states.
- A node state is called only when its maximum posterior probability is at least 0.80. Otherwise it is `UNCERTAIN`.
- A branch transition is called only when both endpoint states meet the 0.80 threshold and differ.
- Results that differ between the taxonomy and gene-tree topologies are `topology-sensitive`.
- Ancestral reconstruction is a model-supported estimate, not an observed historical sequence.

## Recurrence

- NP349 recurrence is evaluated with stochastic character maps on the primary branch-length-optimized taxonomy tree.
- The prespecified reduced state space is Q, K, and O (all other resolved residues). Missing/gap states remain unknown.
- Equal-rates, symmetric, and all-rates-different models are compared by AICc. The lowest-AICc model is primary; models within 2 AICc units form the sensitivity set.
- At least 1,000 stochastic maps are required for final transition-count summaries. A 20-map smoke test must pass first.
- Independent evolution means at least two transitions into the same state on non-ancestral branches in a sampled map. It does not by itself establish adaptive convergence.

## Selection follow-up

- Transition-informed foreground branches are exploratory because the same sequence data identify the state changes.
- Branch-site analyses will be launched only for transitions robust to the 0.80 endpoint rule and both topologies, with exact tip/branch labeling validated first.
- All tested foreground sets enter one declared multiplicity family. A negative result is not evidence that selection never occurred.

## Claim boundary

- Sequence and ancestral analyses can establish lineage structure, model-supported timing, recurrence, and selection-model support.
- They cannot establish altered nucleosome binding, DNMT3A stimulation, methylation, fertility, or adaptation.
