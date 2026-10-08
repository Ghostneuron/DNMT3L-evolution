#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CHECKSUMS = ROOT / "SHA256SUMS.txt"
MANIFEST = ROOT / "MANIFEST.tsv"
PACKAGE_INFO = ROOT / "PACKAGE_INFO.json"
JUNK_PARTS = {".DS_Store", "__pycache__", ".codex_sync_staging"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fail(message: str, failures: list[str]) -> None:
    failures.append(message)
    print(f"FAIL: {message}")


def main() -> int:
    failures: list[str] = []
    required = [
        ROOT / "README.md",
        ROOT / "COMMANDS.md",
        ROOT / "CITATION.cff",
        ROOT / ".zenodo.json",
        ROOT / "RELEASE_NOTES.md",
        ROOT / "RELEASE_CHECKLIST.md",
        ROOT / "LICENSE_SELECTION_REQUIRED.md",
        ROOT / "documents" / "DNMT3L_supplementary_tables_20261007.docx",
        ROOT / "supplementary" / "Table_S1_h358_primary_results.csv",
        ROOT / "supplementary" / "Table_S2_h358_topology_sensitivity.csv",
        ROOT / "supplementary" / "Table_S3_h358_multinucleotide_sensitivity.csv",
        MANIFEST,
        CHECKSUMS,
        PACKAGE_INFO,
    ]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"required file missing or empty: {path.relative_to(ROOT)}", failures)

    excluded_public_files = [
        ROOT / "documents" / "DNMT3L_clade_aware_evolution_manuscript_20261006.docx",
        ROOT / "documents" / "DNMT3L_clade_aware_evolution_manuscript_20261006.pdf",
        ROOT / "scripts" / "build_dnmt3l_broadened_manuscript.py",
    ]
    for path in excluded_public_files:
        if path.exists():
            fail(f"nonpublic manuscript artifact is present: {path.relative_to(ROOT)}", failures)

    try:
        info = json.loads(PACKAGE_INFO.read_text())
    except Exception as exc:
        fail(f"cannot parse PACKAGE_INFO.json: {exc}", failures)
        info = {}

    files = [path for path in ROOT.rglob("*") if path.is_file()]
    for path in files:
        rel = path.relative_to(ROOT)
        if any(part in JUNK_PARTS or part.startswith("._") for part in rel.parts):
            fail(f"unwanted filesystem metadata: {rel}", failures)
        if path.name.startswith("~$") or path.suffix == ".pyc":
            fail(f"unwanted temporary/cache file: {rel}", failures)
        if info.get("bundle") == "github" and path.stat().st_size >= 100_000_000:
            fail(f"GitHub file is at least 100 MB: {rel}", failures)

    if CHECKSUMS.exists():
        listed: set[Path] = set()
        for line_number, line in enumerate(CHECKSUMS.read_text().splitlines(), 1):
            if not line.strip():
                continue
            try:
                expected, relative = line.split("  ", 1)
            except ValueError:
                fail(f"malformed checksum line {line_number}", failures)
                continue
            rel = Path(relative)
            listed.add(rel)
            path = ROOT / rel
            if not path.is_file():
                fail(f"checksum target missing: {rel}", failures)
            elif sha256(path) != expected:
                fail(f"checksum mismatch: {rel}", failures)
        expected_files = {
            path.relative_to(ROOT)
            for path in files
            if path != CHECKSUMS and not path.name.startswith("._")
        }
        if listed != expected_files:
            for rel in sorted(expected_files - listed):
                fail(f"file absent from checksum list: {rel}", failures)
            for rel in sorted(listed - expected_files):
                fail(f"checksum lists unexpected file: {rel}", failures)

    if MANIFEST.exists():
        lines = MANIFEST.read_text().splitlines()
        if not lines or lines[0] != "path\tbytes\tsha256":
            fail("MANIFEST.tsv header is invalid", failures)

    if failures:
        print(f"\nPackage validation failed: {len(failures)} issue(s).")
        return 1
    print(
        f"Package validation passed: {len(files)} files, "
        f"{sum(path.stat().st_size for path in files):,} bytes."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
