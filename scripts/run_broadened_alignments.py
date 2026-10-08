from __future__ import annotations

import os

import csv
import hashlib
import shutil
import subprocess
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "source_package" / "analysis" / "orthology_rebuild_20260826"
INPUTS = PROJECT / "analysis_v2" / "inputs"
OUT = PROJECT / "analysis_v2" / "raw_alignments"
LOGS = PROJECT / "analysis_v2" / "logs" / "alignments"

SCOPES = ("Placentalia", "Euarchontoglires", "Laurasiatheria", "Sauropsida")
ALIGNERS = ("MACSE", "MAFFT", "PRANK")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def archived_alignment(scope: str, aligner: str) -> Path | None:
    base = SOURCE / "clade_alignments" / scope / aligner
    if aligner == "MACSE":
        suffix = "_with_human_anchor_NT.fasta" if scope in {"Laurasiatheria", "Sauropsida"} else "_NT.fasta"
        return base / f"{scope}_MACSE{suffix}"
    if scope == "Placentalia" and aligner == "MAFFT":
        return base / "Placentalia_MAFFT_NT.fasta"
    if scope == "Placentalia" and aligner == "PRANK":
        return base / "Placentalia_PRANK.best.fas"
    return None


def run_command(command: list[str], log_path: Path, stdout_path: Path | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log_handle:
        log_handle.write("COMMAND: " + " ".join(command) + "\n\n")
        if stdout_path is None:
            completed = subprocess.run(
                command,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
        else:
            with stdout_path.open("w") as output_handle:
                completed = subprocess.run(
                    command,
                    stdout=output_handle,
                    stderr=log_handle,
                    text=True,
                    check=False,
                )
        if completed.returncode:
            raise subprocess.CalledProcessError(completed.returncode, command)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, object]] = []

    for scope in SCOPES:
        cds = INPUTS / f"{scope}_cds_with_coordinate_anchor.fasta"
        protein = INPUTS / f"{scope}_protein_with_coordinate_anchor.fasta"
        for aligner in ALIGNERS:
            destination_dir = OUT / scope / aligner
            destination_dir.mkdir(parents=True, exist_ok=True)
            nt_output = destination_dir / f"{scope}_{aligner}_NT.fasta"
            source = archived_alignment(scope, aligner)
            provenance: str
            command_text: str

            if source is not None:
                if not source.exists():
                    raise FileNotFoundError(source)
                shutil.copy2(source, nt_output)
                provenance = "copied_from_archived_reanalysis"
                command_text = f"cp {source.relative_to(PROJECT)} {nt_output.relative_to(PROJECT)}"
            elif aligner == "MAFFT":
                aa_output = destination_dir / f"{scope}_MAFFT_AA.fasta"
                mafft_command = [os.environ.get("MAFFT_BIN", "mafft"), "--auto", "--thread", "8", str(protein)]
                run_command(
                    mafft_command,
                    LOGS / f"{scope}_MAFFT.log",
                    stdout_path=aa_output,
                )
                backtranslate_command = [
                    os.environ.get("MACSE_BIN", "macse"),
                    "-prog",
                    "reportGapsAA2NT",
                    "-align_AA",
                    str(aa_output),
                    "-seq",
                    str(cds),
                    "-out_NT",
                    str(nt_output),
                ]
                run_command(backtranslate_command, LOGS / f"{scope}_MAFFT_backtranslate.log")
                provenance = "new_clade_specific_alignment"
                command_text = " ; ".join((" ".join(mafft_command), " ".join(backtranslate_command)))
            else:
                prefix = destination_dir / f"{scope}_PRANK"
                prank_command = [
                    os.environ.get("PRANK_BIN", "prank"),
                    f"-d={cds}",
                    f"-o={prefix}",
                    "-codon",
                    "-f=fasta",
                    "-iterate=2",
                    "-seed=20261006",
                    "-quiet",
                ]
                run_command(prank_command, LOGS / f"{scope}_PRANK.log")
                best = Path(f"{prefix}.best.fas")
                if not best.exists():
                    raise FileNotFoundError(best)
                shutil.copy2(best, nt_output)
                provenance = "new_clade_specific_alignment"
                command_text = " ".join(prank_command)

            manifest_rows.append(
                {
                    "Scope": scope,
                    "Aligner": aligner,
                    "Provenance": provenance,
                    "Input_CDS": str(cds.relative_to(PROJECT)),
                    "Input_Protein": str(protein.relative_to(PROJECT)),
                    "Output_Path": str(nt_output.relative_to(PROJECT)),
                    "Output_SHA256": sha256(nt_output),
                    "Command": command_text,
                }
            )
            print(f"{scope} {aligner}: {provenance}", flush=True)

    manifest = OUT / "alignment_manifest.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]))
        writer.writeheader()
        writer.writerows(manifest_rows)


if __name__ == "__main__":
    main()
