# Installation

## Status and supported use

The installation metadata is intended for source-based development and
inspection of the recovered implementation. The project has not yet established
a tested compatibility matrix or a stable release. Python 3.10 or newer is
declared provisionally in `pyproject.toml`.

## Create an isolated environment

From the repository root:

```console
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install --editable .
```

This installs both packages from the shared `src/` tree:

- `pymatmc2`, the MC² implementation; and
- `mexm`, its bundled historical materials-simulation dependency.

For plotting support:

```console
python -m pip install --editable '.[plot]'
```

For development tools:

```console
python -m pip install --editable '.[dev]'
```

Both extras can be installed together:

```console
python -m pip install --editable '.[dev,plot]'
```

## Verify imports

```console
python -c "import pymatmc2; import mexm; print('imports succeeded')"
```

A successful import is a packaging smoke check only. It does not establish that
a simulation workflow is correct or ready to run.

## Dependencies

The active `pymatmc2` import path requires:

- ASE;
- NumPy;
- SciPy;
- pandas;
- statsmodels; and
- PyYAML.

Matplotlib is optional and used by plotting behavior. PostgreSQL support
belongs to a broader historical `mexm` surface and is available through the
`legacy-postgresql` extra. Other historical `mexm` modules may have additional
unmaintained dependencies.

`MEXM_REQUIREMENTS.historical.txt` is an archival copy of the upstream 2020
requirements. It is not the installation authority and should not be installed
as a lockfile. Current installation dependencies are declared only in
`pyproject.toml`.

## External software is not installed

Python installation does not provide:

- VASP or a VASP license;
- VASP pseudopotentials;
- LAMMPS;
- GULP;
- PhonTS;
- MPI launchers;
- Torque/PBS, Slurm, or another scheduler; or
- site-specific modules, executable paths, accounts, and credentials.

Do not attempt an external-calculator or scheduler workflow until its inputs,
software license, computational cost, and execution environment have been
reviewed and explicitly authorized.

## Known limitations

- The package is currently marked `0.1`; this identifies the recovered and
  installable source baseline, not a scientifically validated release.
- The recovered test tree is historical and is not yet a clean modern test
  suite.
- Some historical `mexm` modules are incomplete or rely on unavailable legacy
  modules. Five files that do not parse are retained under
  `references/mexm-broken/` instead of the active package.
- Successful installation does not imply numerical verification or scientific
  validation.
