# pymatmc2

`pymatmc2` is my Python implementation of the Multi-Cell Monte Carlo
(MC²) approach developed by Maryam Ghazisaeidi and collaborators [1–3].
The method searches the compositional space of competing crystalline phases to
study phase stability and coexistence in alloys.

## Project context

I was a postdoctoral researcher at The Ohio State University from 2020 to 2021,
during the COVID-19 pandemic. I have restored and updated this repository in
preparation for returning to this problem and exploring additional mechanisms
within the MC² framework. That future development is planned work; the
current repository remains a recovered historical implementation and should
not be interpreted as evidence that those mechanisms are already implemented
or validated.

## Installation

Installation is currently intended for source-based development rather than a
stable release:

```console
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install --editable .
```

Package and dependency metadata are maintained in `pyproject.toml`. See
[`docs/installation.md`](docs/installation.md) for optional dependencies,
verification, external-software requirements, and known limitations.

The current source baseline is **v0.1**. This version identifies the recovered,
installable starting point for renewed development; it does not represent
scientific validation of the implementation.

## Documentation

- [Documentation index](docs/index.md)
- [Installation](docs/installation.md)
- [Repository structure](docs/repository-structure.md)
- [Development guide](docs/development.md)
- [Scientific background](docs/scientific-background.md)

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

## Bundled historical dependency

`pymatmc2` directly depends on the historical `mexm` project for structure
representations, VASP I/O, simulation objects, configuration parsing, and job
submission. A source-only snapshot of `mexm` commit
`3963a1a30ba1595fba69cd0f791dc31fe2e08f7c` is therefore bundled at
`src/mexm/`.

The bundle contains 133 checksum-verified package files. Five syntactically
incomplete upstream files are preserved separately under
`references/mexm-broken/`. See `MEXM_PROVENANCE.md` and
`MEXM_SOURCE_SHA256SUMS` for the exact boundary and identities. Historical
package requirements are retained in `MEXM_REQUIREMENTS.historical.txt`; they
have not been resolved or validated against current Python versions.

## Licensing

The recovered `pymatmc2` source is distributed under the MIT License; see
`LICENSE`. The bundled historical `mexm` source is also MIT-licensed under its
own retained notice in `MEXM_LICENSE`.

The non-operational historical reference under `references/mc2/` is separate.
It records an unsuccessful MPI-spawn experiment derived from GPL-licensed
abICS code and is distributed under GPL-3.0-or-later. It is not part of the
`pymatmc2` package and is not imported by the recovered source. Historical
references under `references/broken/` remain covered by the MIT License but are
explicitly excluded from supported source.

## References

1. C. Niu, W. Windl, and M. Ghazisaeidi, “Multi-Cell Monte Carlo Relaxation
   method for predicting phase stability of alloys,” *Scripta Materialia*
   **132**, 9–12 (2017).
   <https://doi.org/10.1016/j.scriptamat.2017.01.001>
2. C. Niu, Y. Rao, W. Windl, and M. Ghazisaeidi, “Multi-cell Monte Carlo method
   for phase prediction,” *npj Computational Materials* **5**, 120 (2019).
   <https://doi.org/10.1038/s41524-019-0259-z>
3. E. Antillon and M. Ghazisaeidi, “Efficient determination of solid-state
   phase equilibrium with the multicell Monte Carlo method,” *Physical Review
   E* **101**, 063306 (2020).
   <https://doi.org/10.1103/PhysRevE.101.063306>

These papers define and develop the scientific method that motivated this
implementation. This repository does not imply endorsement by their authors.
