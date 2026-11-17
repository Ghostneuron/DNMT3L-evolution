# Phylogenetically balanced DNMT3L structure pilot

Frozen: 2026-11-04

## Estimand

Test whether DNMT3L structured domains or its C-terminal complex with DNMT3A show reproducible, confidence-qualified geometry differences across major tetrapod lineages and the reconstructed Amniota-to-Theria contrast.

This is a structural-hypothesis screen. It cannot establish affinity, catalysis, germline deployment, genomic imprinting, or adaptation.

## Panel

The primary panel contains 22 extant orthologues plus reconstructed maximum-posterior Amniota and Theria DNMT3L sequences (24 inputs total):

- 3 amphibians: a caecilian and two salamandrid newts;
- 6 non-avian sauropsids: crocodilian, turtle, agamid, phrynosomatid lizard, gecko, and snake;
- 4 marsupials: American and Australasian representatives;
- 9 placentals spanning Euarchontoglires, Laurasiatheria, Afrotheria, and Xenarthra;
- Amniota and Theria ancestral maximum-posterior sequences from the validated taxonomy-topology IQ-TREE reconstruction.

Selection is frozen in `config/panel_species.tsv`. All extant inputs must pass the existing DNMT3L orthology/sequence gates. Complex-eligible extant inputs must also occur in the validated matched DNMT3A-DNMT3L inventory.

## Ancestral sequence rule

- Amniota is IQ-TREE node `NCBI_32524`; Theria is `NCBI_32525`.
- The prediction input uses the maximum-posterior residue at every human-coordinate site.
- Posterior probability is retained per position. Sites below 0.80 are uncertain and are excluded from residue-level or local structural comparisons.
- Domain-wide comparisons report both all mapped residues and the subset with ancestral posterior probability at least 0.80.
- A maximum-posterior sequence is a model-derived point estimate, not an observed ancestral protein.

## Prediction stages

### Monomer

- Extant inputs: complete curated DNMT3L protein selected for the orthology analysis.
- Ancestral inputs: 387-position human-coordinate maximum-posterior sequence.
- Primary comparisons: ADD domain and C-terminal domain. The N-terminal region is retained so that whole-monomer context can influence the prediction, but low-confidence N-terminal geometry is not interpreted.

### Complex

- Primary complex construct: human DNMT3A1 residues 623-912 held constant plus the DNMT3L C-terminal region corresponding to human NP_037501.2 residues 178-380.
- Holding DNMT3A constant isolates the effect of DNMT3L sequence state and avoids confounding by partner substitutions.
- Species-specific DNMT3A complexes are secondary and are admitted only if the fixed-partner screen identifies a reproducible contrast.

## Identical policy

- Predictor: Boltz Community 2.13.0 with the Boltz-2 checkpoint for every input; the community implementation is used for Apple-MPS compatibility and is not treated as an independent predictor.
- Templates: none for all inputs.
- MSA: fully local, frozen project alignments for every protein. DNMT3L uses the validated 289-orthologue human-coordinate alignment projected into each query; DNMT3A uses the validated MACSE alignment projected onto human residues 623-912. Matched species share pairing keys and unmatched sequences remain unpaired. No target sequence is sent to an external MSA service.
- Seeds: 4101, 4102, and 4103, one diffusion sample per seed.
- Recycling steps: 3; sampling steps: 200; step scale: predictor default for the pinned checkpoint.
- Output: mmCIF, confidence JSON, full PAE, and full PDE.

## Staged compute gate

1. Preflight all sequences, hashes, storage, predictor version, and hardware.
2. Run one-seed smoke predictions for human and ancestral-Amniota monomers and complexes.
3. Validate output parsing and benchmark elapsed time/memory.
4. Run three seeds for the sentinel panel: Geotrypetes, Alligator, Monodelphis, human, mouse, Amniota, and Theria.
5. Expand monomers to all 24 inputs only if the smoke and sentinel gates pass.
6. Expand complexes only for sequence contrasts whose monomer or sentinel-complex results are reproducible and confidence-qualified.

## Confidence and reproducibility gates

### Monomer domain

- At least 80% of mapped domain residues have pLDDT at least 70.
- At least two of three seeds pass the domain-confidence gate.
- Within-input domain C-alpha RMSD across passing seeds is at most 2.0 Å after domain superposition.

### Complex interface

- At least two of three seeds have ipTM at least 0.70.
- Interface residues have median pLDDT at least 70 and median inter-chain PAE at most 10 Å.
- The same chain orientation/interface is recovered in at least two seeds.

### Comparative difference

- A lineage contrast is called reproducible only when the between-input domain/interface difference exceeds two times the larger within-input seed dispersion and the same direction is present in at least two matched seed comparisons.
- Local comparisons around focal positions use only alignment-mapped residues with pLDDT at least 70, and for ancestral inputs posterior probability at least 0.80.

## MD admission

MD is not launched from one attractive prediction. A comparison must first:

1. pass the relevant confidence gate in at least two seeds;
2. show a reproducible between-input difference larger than seed variability;
3. remain supported after excluding uncertain ancestral residues and low-confidence structural neighborhoods; and
4. define one prespecified geometry/contact hypothesis that can be simulated with matched preparation and at least three independent replicas.

The already-running 2QRV-based therian-interface MD analysis is independent and is not retroactively treated as validation of this screen.

