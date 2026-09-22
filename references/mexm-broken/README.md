# Incomplete historical mexm references

This directory preserves five `mexm` package files from upstream commit
`3963a1a30ba1595fba69cd0f791dc31fe2e08f7c` that do not compile. They are
retained at their original repository-relative paths but are not admitted into
the active package under `src/mexm/`.

## Dispositions

- `src/mexm/io/hpcutil.py`: empty class body causing `IndentationError`.
- `src/mexm/io/phonopy/__init__.py`: malformed function call causing
  `SyntaxError`.
- `src/mexm/io/slurm/daemon.py`: empty conditional body causing
  `IndentationError`.
- `src/mexm/io/torque.py`: malformed expression causing `SyntaxError` and
  references to modules absent from the selected source tree.
- `src/mexm/simulation/multicell.py`: empty loop body causing
  `IndentationError`.

`REFERENCE_SHA256SUMS` records the exact retained bytes. These files remain
covered by the MIT terms in `MEXM_LICENSE` at the repository root.
