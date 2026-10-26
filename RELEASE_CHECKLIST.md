# Publication checklist

Complete the unchecked items before publishing the Zenodo record.

- [x] Current manuscript Word and PDF files excluded from public bundles
- [x] Manuscript, response-letter, and journal-strategy generation scripts excluded
- [x] Formatted supplementary tables included
- [x] Machine-readable supplementary tables included
- [x] Main and supplementary figures included in PNG and PDF
- [x] Public scripts, environment records, provenance notes, manifests, and checksums included
- [x] Private journal-strategy and superseded submission documents excluded
- [x] macOS metadata, caches, lock files, and staging directories excluded
- [x] Package validator passes
- [x] Archive integrity test passes
- [x] Add MIT licensing for original code and CC BY 4.0 licensing for other original content
- [x] Create the public GitHub repository and upload the contents of the GitHub package directory
- [x] Create a GitHub release tagged `v1.0.0`
- [x] Preserve the immutable `v1.0.0` tag and stage the additive `v1.1.0` update
- [x] Publish `v1.1.1` to restore three checksum-listed logs omitted from `v1.1.0`
- [x] Add revision-specific analyses under `revision_20261025/`
- [x] Update Figures 3 and 4 and the current supplementary-table artifacts
- [x] Audit newly added v1.1.0 files for private correspondence, credentials, absolute paths, and filesystem metadata
- [ ] Reserve or mint the Zenodo DOI and upload the Zenodo ZIP archive
- [x] Add the GitHub URL to `README.md` and `CITATION.cff`
- [ ] Add the Zenodo DOI to `README.md`, `CITATION.cff`, and the manuscript Data availability statement
- [x] Update `.zenodo.json` with the path-specific mixed-license designation
- [x] Rebuild manifests and archives after license edits
- [x] Rerun `python scripts/validate_public_package.py` on the final uploaded contents
- [ ] Verify the GitHub release and Zenodo deposit are mutually linked

For journal submission, also confirm whether the journal requires city and country in the unaffiliated author address.
