#!/usr/bin/env bash
set -euo pipefail

# Run from the Vertebrate_Gene_History_Redesign_20261102 project root.
ROOT="$(pwd)"
INPUT="$ROOT/analysis/codon_alignments_337/inputs"
OUT="$ROOT/analysis/codon_alignments_337"

macse -prog alignSequences -seq "$INPUT/Placentalia_cds_with_coordinate_anchor.fasta" -out_NT "$OUT/Placentalia/MACSE/Placentalia_MACSE_NT.fasta" -out_AA "$OUT/Placentalia/MACSE/Placentalia_MACSE_AA.fasta" -gc_def 1 -score_matrix BLOSUM62 -optim 0
mafft --auto --thread 8 "$INPUT/Placentalia_protein_with_coordinate_anchor.fasta" > "$OUT/Placentalia/MAFFT/Placentalia_MAFFT_AA.fasta"
macse -prog reportGapsAA2NT -align_AA "$OUT/Placentalia/MAFFT/Placentalia_MAFFT_AA.fasta" -seq "$INPUT/Placentalia_cds_with_coordinate_anchor.fasta" -out_NT "$OUT/Placentalia/MAFFT/Placentalia_MAFFT_NT.fasta"
prank -d="$INPUT/Placentalia_cds_with_coordinate_anchor.fasta" -o="$OUT/Placentalia/PRANK/Placentalia_PRANK" -codon -f=fasta -iterate=2 -seed=20261102 -quiet

MUSCLE5="${MUSCLE5:-muscle}"
COBALT="${COBALT:-cobalt}"
"$MUSCLE5" -super5 "$INPUT/Placentalia_protein_with_coordinate_anchor.fasta" -output "$OUT/Placentalia/MUSCLE5/Placentalia_MUSCLE5_AA.fasta" -threads 8
macse -prog reportGapsAA2NT -align_AA "$OUT/Placentalia/MUSCLE5/Placentalia_MUSCLE5_AA.fasta" -seq "$INPUT/Placentalia_cds_with_coordinate_anchor.fasta" -out_NT "$OUT/Placentalia/MUSCLE5/Placentalia_MUSCLE5_NT.fasta"
"$COBALT" -i "$INPUT/Placentalia_protein_with_coordinate_anchor.fasta" -norps T -outfmt mfasta > "$OUT/Placentalia/COBALT/Placentalia_COBALT_AA.raw.fasta"
python scripts/normalize_cobalt_and_project.py \
  --protein "$INPUT/Placentalia_protein_with_coordinate_anchor.fasta" \
  --cds "$INPUT/Placentalia_cds_with_coordinate_anchor.fasta" \
  --raw-alignment "$OUT/Placentalia/COBALT/Placentalia_COBALT_AA.raw.fasta" \
  --named-alignment "$OUT/Placentalia/COBALT/Placentalia_COBALT_AA.fasta" \
  --codon-alignment "$OUT/Placentalia/COBALT/Placentalia_COBALT_NT.fasta" \
  --qc-json "$OUT/Placentalia/COBALT/Placentalia_COBALT_QC.json"

for CLADE in Euarchontoglires Laurasiatheria Sauropsida Amphibia; do
  macse -prog alignSequences -seq "$INPUT/${CLADE}_cds_with_coordinate_anchor.fasta" -out_NT "$OUT/${CLADE}/MACSE/${CLADE}_MACSE_NT.fasta" -out_AA "$OUT/${CLADE}/MACSE/${CLADE}_MACSE_AA.fasta" -gc_def 1 -score_matrix BLOSUM62 -optim 0
done
