# Incomplete historical source references

This directory preserves historical source that was present at upstream commit
`3063f031301c04e395cb9d7bc5064f2a52e1110c` but is not admitted into the
recovered package or maintained test tree.

The files are retained as source evidence only. They are not imported,
collected as tests, or represented as working software.

## Dispositions

- `src/pymatmc2/widom.py`: incomplete implementation; it contains an empty
  method body and does not compile.
- `tests/unittests/WidomTest/`: development inputs for the incomplete Widom
  implementation.
- `tests/unittests/MultiCell/test__after_first_iteration/`: experimental local
  implementation containing invalid Python syntax.
- `tests/unittests/Pymatmc2Results/data.py`: incomplete experimental result
  writer containing an empty method body.
- `tests/unittests/Pymatmc2Results/results.py`: incomplete experimental result
  reader containing multiple syntax errors.

`REFERENCE_SHA256SUMS` records the retained reference bytes. Their original
repository-relative paths are preserved beneath this directory.
