# Alignment decision for the 337-GeneID update

## Decision

The alignments were rebuilt because the current NCBI package changes both taxon membership and representative accessions. Relative to the previous analysis, the final set increases from 269 to 289 records, 21 GeneIDs gain analysis-ready status, one loses it, and 31 of the 268 ready GeneIDs shared between versions use a different protein accession.

## Five-method placental panel

The 220 placental records were aligned with MACSE, MAFFT, PRANK, MUSCLE5, and COBALT. MACSE is the primary codon-aware workflow. PRANK is a second codon-aware sensitivity analysis. MAFFT, MUSCLE5, and COBALT align proteins and then project the exact paired CDS codons. COBALT was run in local-similarity-only mode (`-norps T`) so that the analysis is reproducible without a local Conserved Domain Database.

The five methods are technical sensitivity analyses on the same biological records. They are not biological replicates, and no majority-vote rule is used.

## Audit outcome

- All methods contain the same 220 taxa and preserve the exact source CDS content, allowing for optional removal of terminal stop codons.
- Each method maps all 387 residues of human NP_037501.2 and assigns alanine at h358.
- Raw MACSE output contains 15 frameshift markers in five records; these are alignment symbols rather than source-sequence frameshifts and are treated as unresolved during consensus filtering.
- The shared 80% occupancy gate retains 367 of 387 human-reference positions.
- Exact five-method codon consensus masks 2,460 of 80,740 retained taxon-position calls (3.05%).
- The five resulting strict-consensus files are byte-identical by design and constitute one conservative downstream alignment, not five independent selection tests.

At prespecified functional sites, all five methods agree for every fully resolved placental taxon at human positions 88, 90, 105-107, 204, 377, and 358. Switching-helix positions 349, 351, 352, and 355 pass the occupancy gate but require limited taxon-specific masking. h358 is retained with a minimum method occupancy of 85.45%; 182 taxa are resolved across all five methods and have identical amino-acid calls.

## Clade-specific alignments

Separate MACSE alignments were generated for Euarchontoglires, Laurasiatheria, Sauropsida, and Amphibia. A human coordinate anchor was added to clades that do not contain human and must be removed before clade-specific evolutionary modeling. These alignments avoid assuming that one deep-vertebrate alignment is uniformly reliable.
