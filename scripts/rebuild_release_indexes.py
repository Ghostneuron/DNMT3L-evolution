#!/usr/bin/env python3
"""Rebuild release manifests and SHA-256 indexes from the public tree."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REVISION = ROOT / "revision_20261025"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def public_files(root: Path, excluded: set[Path]) -> list[Path]:
    return sorted(
        (
            path
            for path in root.rglob("*")
            if path.is_file()
            and ".git" not in path.relative_to(root).parts
            and path not in excluded
        ),
        key=lambda path: path.relative_to(root).as_posix(),
    )


def write_atomic(path: Path, text: str) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    temporary.replace(path)


def main() -> None:
    revision_index = REVISION / "SHA256SUMS.tsv"
    revision_files = public_files(REVISION, {revision_index})
    revision_lines = ["path\tbytes\tsha256"]
    revision_lines.extend(
        f"{path.relative_to(REVISION).as_posix()}\t{path.stat().st_size}\t{sha256(path)}"
        for path in revision_files
    )
    write_atomic(revision_index, "\n".join(revision_lines) + "\n")

    manifest = ROOT / "MANIFEST.tsv"
    checksums = ROOT / "SHA256SUMS.txt"
    manifest_files = public_files(ROOT, {manifest, checksums})
    manifest_lines = ["path\tbytes\tsha256"]
    manifest_lines.extend(
        f"{path.relative_to(ROOT).as_posix()}\t{path.stat().st_size}\t{sha256(path)}"
        for path in manifest_files
    )
    write_atomic(manifest, "\n".join(manifest_lines) + "\n")

    checksum_files = public_files(ROOT, {checksums})
    checksum_lines = [
        f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}" for path in checksum_files
    ]
    write_atomic(checksums, "\n".join(checksum_lines) + "\n")
    print(
        f"Indexed {len(manifest_files)} release files and "
        f"{len(revision_files)} revision-specific files."
    )


if __name__ == "__main__":
    main()
