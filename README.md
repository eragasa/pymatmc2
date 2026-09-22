# pymatmc2

`pymatmc2` is a historical Python implementation of a multi-cell Monte Carlo
workflow for searching compositional space among competing phases.

## Recovery status

This branch is a source-only recovery from historical `develop` commit
`3063f031301c04e395cb9d7bc5064f2a52e1110c`. The recovered boundary contains
package source, source-bearing tests, documentation, and scripts whose exact
pre-recovery identities are listed in `SOURCE_SHA256SUMS`.

The historical `develop` tree also contained gigabytes of generated calculator
and simulation output. Those outputs are deliberately not restored here. In
particular, this recovery excludes recorded `OUTCAR`, `vasprun.xml`, restart,
wavefunction, generated-result, and other calculation-artifact trees. See
`PROVENANCE.md` for the exact recovery and exclusion boundary.

Known incomplete or syntactically invalid historical files are preserved under
`references/broken/` rather than admitted into the package or maintained test
tree. The compact implementation from the initial `master` commit is retained
under `references/master-initial/`; its pseudopotential file is omitted because
redistribution authority was not established. The restored code remains
historical research software: recovery preserves source but does not establish
current executability, numerical verification, scientific validation, or
support for production calculations.

## Licensing

The recovered `pymatmc2` source is distributed under the MIT License; see
`LICENSE`.

The non-operational historical reference under `references/mc2/` is separate.
It records an unsuccessful MPI-spawn experiment derived from GPL-licensed
abICS code and is distributed under GPL-3.0-or-later. It is not part of the
`pymatmc2` package and is not imported by the recovered source. Historical
references under `references/broken/` remain covered by the MIT License but are
explicitly excluded from supported source.

## Scientific attribution

The multi-cell Monte Carlo method is attributed to:

C. Niu, Y. Rao, W. Windl, and M. Ghazisaeidi, “Multi-cell Monte Carlo method
for phase prediction,” *npj Computational Materials* **5**, 120 (2019).
<https://doi.org/10.1038/s41524-019-0259-z>

This citation identifies the scientific method. It does not imply endorsement
of this historical implementation by the article's authors.
