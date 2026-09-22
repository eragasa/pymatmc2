# Source-recovery provenance

## Historical source identity

- Repository: `https://github.com/eragasa/pymatmc2`
- Historical branch: `develop`
- Source commit: `3063f031301c04e395cb9d7bc5064f2a52e1110c`
- Source tree: `54866bd536c8131b1932354d61f1ff33ae881a69`
- Commit date: `2020-10-09T17:22:25-04:00`
- Recovery branch: `master`
- Recovery date: `2026-09-23`

The source was recovered from the byte-verified Project Koios preservation
snapshot corresponding to that exact Git commit. `SOURCE_SHA256SUMS` records
the SHA-256 identity of each recovered file before this recovery commit;
`SOURCE_GIT_MODES` records the corresponding historical Git modes.

## Included boundary

The active recovery includes 143 checksum-verified source files comprising:

- the `pymatmc2` Python package;
- Python and MATLAB test source;
- source-bearing scripts;
- Sphinx source and diagrams; and
- compact authored configuration or test inputs selected by the preservation
  snapshot.

Repository-relative paths are retained.

## Excluded boundary

The recovery excludes historical calculator data and generated artifacts,
including large VASP output trees, wavefunctions, restart data, generated
results, caches, and databases. The historical `develop` Git object remains the
identity for those omitted files; omission does not assert that they are
unimportant or independently reproducible.

The initial `master` implementation is retained separately under
`references/master-initial/` rather than mixed into the later `develop`
package. Its tracked `POTCAR` is removed from the recovered tip because
pseudopotential redistribution authority was not established; the reference
record retains its Git-blob and SHA-256 identities without duplicating the
bytes.

`src/pymatmc2/multicellmc_spawn.py` is also excluded from the MIT package. It
was an unsuccessful experiment substantially derived from the GPL-3.0-or-later
abICS `run_base_mpi.py` implementation. A clearly separated historical
reference and its licensing record are retained under `references/mc2/`.

Nine additional historical files are excluded from active source because they
are incomplete, contain invalid Python syntax, or depend directly on those
incomplete experiments. They remain byte-preserved under `references/broken/`,
with their dispositions and identities recorded there.

## Evidentiary limits

Source recovery and checksum agreement establish file identity only. They do
not establish runtime compatibility, correctness, numerical verification,
scientific validation, or authorization to execute external calculators.
