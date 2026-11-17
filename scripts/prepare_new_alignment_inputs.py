#!/usr/bin/env python3
from __future__ import annotations

import shutil
from pathlib import Path

from Bio import SeqIO


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "data/processed/ncbi_337_representatives"
OUT = PROJECT / "analysis/codon_alignments_337/inputs"


def write_records(records, path: Path) -> None:
    with path.open("w") as handle:
        SeqIO.write(records, handle, "fasta")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_cds = list(SeqIO.parse(SOURCE / "DNMT3L_337_selection_ready_all_cds.fasta", "fasta"))
    all_protein = list(SeqIO.parse(SOURCE / "DNMT3L_337_selection_ready_all_protein.fasta", "fasta"))
    human_cds = [record for record in all_cds if record.id == "Homo_sapiens"]
    human_protein = [record for record in all_protein if record.id == "Homo_sapiens"]
    if len(human_cds) != 1 or len(human_protein) != 1:
        raise ValueError("Expected one human CDS and protein anchor")
    human_cds[0].id = "Homo_sapiens_anchor"
    human_cds[0].name = "Homo_sapiens_anchor"
    human_cds[0].description = "human DNMT3L NP_037501.2 coordinate anchor"
    human_protein[0].id = "Homo_sapiens_anchor"
    human_protein[0].name = "Homo_sapiens_anchor"
    human_protein[0].description = "human DNMT3L NP_037501.2 coordinate anchor"

    subsets = ["Placentalia", "Euarchontoglires", "Laurasiatheria", "Sauropsida", "Amphibia"]
    rows = []
    for subset in subsets:
        cds_source = SOURCE / f"DNMT3L_337_selection_ready_{subset}_cds.fasta"
        protein_source = SOURCE / f"DNMT3L_337_selection_ready_{subset}_protein.fasta"
        cds_records = list(SeqIO.parse(cds_source, "fasta"))
        protein_records = list(SeqIO.parse(protein_source, "fasta"))
        anchor_required = subset in {"Laurasiatheria", "Sauropsida", "Amphibia"}
        if anchor_required:
            cds_records.append(human_cds[0])
            protein_records.append(human_protein[0])
        cds_output = OUT / f"{subset}_cds_with_coordinate_anchor.fasta"
        protein_output = OUT / f"{subset}_protein_with_coordinate_anchor.fasta"
        write_records(cds_records, cds_output)
        write_records(protein_records, protein_output)
        # Unanchored copies are explicit inputs for analyses that do not need coordinate mapping.
        shutil.copyfile(cds_source, OUT / f"{subset}_cds.fasta")
        shutil.copyfile(protein_source, OUT / f"{subset}_protein.fasta")
        rows.append(
            {
                "subset": subset,
                "biological_records": len(cds_records) - int(anchor_required),
                "coordinate_anchor_added": anchor_required,
                "alignment_input_records": len(cds_records),
                "cds_input": str(cds_output.relative_to(PROJECT)),
                "protein_input": str(protein_output.relative_to(PROJECT)),
            }
        )

    import pandas as pd

    pd.DataFrame(rows).to_csv(
        PROJECT / "results/ncbi_337_alignment_input_manifest.tsv", sep="\t", index=False
    )
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
