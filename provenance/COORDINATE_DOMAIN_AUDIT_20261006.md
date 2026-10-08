# DNMT3L coordinate and domain-boundary audit

Date: 2026-10-06

## Finding

The structural DNMT3L construct is annotated in UniProt Q9UJW3 coordinates as residues 178-379. Direct sequence crosswalking showed that NP_037501.2 contains an additional serine at NP333 relative to Q9UJW3. The structural interval therefore maps to NP_037501.2 residues 178-380, not 178-379.

NP332-333 is a tandem `SS` run aligned to one Q9UJW3 serine, so an unconstrained global alignment can place the single-residue gap at either NP332 or NP333 with the same score. The released crosswalk uses the established rightmost convention and labels NP333 as the additional residue. This tie choice has no effect on any downstream coordinate, including NP349, NP352, NP358, or NP380.

The same offset explains why the experimentally annotated UniProt switching helix E340-K354 maps to NP341-355 and why UniProt Q348 and Q351 map to NP349 and NP352.

## Correction applied

- The C-terminal methyltransferase-like domain label was changed from NP178-379 to NP178-380.
- Entropy-based and tree-aware equal-length C-terminal background windows now end at NP380.
- Domain summaries were regenerated after the correction.
- The archived record-quality rule remains described as at least 90% coverage of NP178-379 because that was the interval actually used for the archived eligibility check.

## Analyses unaffected

The one-residue endpoint correction does not change any raw or filtered alignment sequence, occupancy decision, duplicate collapse, selection input, inferred tree, species topology, human-referenced site coordinate, switching-helix coordinate, adjacent-loop coordinate, or HyPhy model result. Selection-input hashes were compared before and after regeneration to verify this expectation.

## Sources and reproducible crosswalk

- PDBe entry 8XEE identifies the DNMT3L construct as UniProt Q9UJW3 residues 178-379: https://www.ebi.ac.uk/pdbe/entry/pdb/8xee
- The archived project contains NP_037501.2 and Q9UJW3 sequences used for direct global alignment and residue crosswalking.
- PDB 9MPP/SIFTS annotations use Q9UJW3 numbering for chain N; the PyMOL structural selection therefore correctly remains `resi 315-379` in structure coordinates.
- `analysis_v2/summaries/NP_037501_2_to_9MPP_Q9UJW3_crosswalk.csv` is generated directly from the archived human NP sequence and the chain-N/Q9UJW3 entity sequence embedded in `analysis_v2/structure/9MPP.cif`.
