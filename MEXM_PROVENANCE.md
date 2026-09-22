# Bundled mexm provenance

## Historical identity

- Component: `mexm` (“Materials Ex Machina”)
- Upstream repository: <https://github.com/eragasa/mexm-base>
- Historical branch: `dev`
- Source commit: `3963a1a30ba1595fba69cd0f791dc31fe2e08f7c`
- Source tree: `a19cd86088aad5ecb0b3bbf8a4c73e0159a22fd3`
- Commit date: `2020-07-20T15:38:43-04:00`
- License: MIT; retained in `MEXM_LICENSE`

The source was copied from the byte-verified Project Koios preservation
snapshot for the exact commit above. `MEXM_SOURCE_SHA256SUMS` records the
original SHA-256 identity of each active bundled file, and
`MEXM_SOURCE_GIT_MODES` records its historical Git mode.

## Why mexm is bundled

The recovered `pymatmc2` implementation imports `mexm` directly for crystal
structures, VASP file I/O, simulation representations, filesystem parsing,
exceptions, and job-submission objects. Bundling the historical dependency at
`src/mexm/` preserves that source relationship without relying on the continued
availability or mutability of another repository.

This is source preservation, not a claim that the combined tree is a modern,
supported, or independently validated distribution.

## Included boundary

The active bundle contains 133 checksum-verified Python files from the
historical `src/mexm/` package. Repository-relative package paths are retained.
The upstream test suite, scripts, editor configuration, generated calculator
outputs, pseudopotentials, databases, and Git metadata are not copied into this
repository.

Five package files that do not compile are preserved instead under
`references/mexm-broken/`. They are historical evidence, not active package
source.

## Historical dependencies

`MEXM_REQUIREMENTS.historical.txt` preserves the exact upstream
`requirements.txt` from the selected commit. It is an archival dependency
record, not a lockfile or a statement of compatibility with current package
versions. The `mexm` paths imported by `pymatmc2` require NumPy and PyYAML;
other bundled modules reference additional packages from the historical list,
including SciPy, pandas, ASE, scikit-learn, and database adapters.

No dependency was installed or executed as part of this recovery.

## Evidentiary limits

Checksum agreement establishes source identity only. Successful syntax
compilation establishes only that the active Python files parse under the
interpreter used for recovery. Neither establishes runtime compatibility,
scientific correctness, numerical verification, calculator availability, or
authorization to submit external jobs.
