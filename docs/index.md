# pymatmc2 documentation

`pymatmc2` is a recovered Python implementation of the Multi-Cell Monte Carlo
(MC²) method for studying phase stability and coexistence in alloys. The
repository is being prepared for renewed development, but the recovered code is
not yet presented as a validated production package.

## Start here

- [Installation](installation.md)
- [Repository structure](repository-structure.md)
- [Development guide](development.md)
- [Scientific background](scientific-background.md)

## Status

The active tree contains:

- the historical `pymatmc2` implementation under `src/pymatmc2/`;
- the historical `mexm` dependency under `src/mexm/`;
- source-oriented legacy tests under `tests/`; and
- provenance and separately retained incomplete material.

Installation makes the Python packages importable. It does not install or
license external calculators, supply pseudopotentials, configure an HPC
scheduler, reproduce a published calculation, or establish scientific
validation.

For exact recovery boundaries, see [`PROVENANCE.md`](../PROVENANCE.md) and
[`MEXM_PROVENANCE.md`](../MEXM_PROVENANCE.md).
