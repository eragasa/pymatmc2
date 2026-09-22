# Development guide

## Set up the project

Create an isolated environment and install the editable package with development
and plotting dependencies:

```console
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install --editable '.[dev,plot]'
```

All maintained package and dependency metadata belongs in `pyproject.toml`.

## Current baseline

This repository began as a source recovery. Before changing behavior, determine
whether the relevant code is:

- active `pymatmc2` source;
- bundled historical `mexm` source;
- an incomplete reference retained under `references/`; or
- a development script in the recovered test tree.

The historical checksum manifests record original identities. Do not regenerate
them to make an edited file appear unchanged.

## Checks

Run syntax compilation without writing bytecode into the repository:

```console
cache_dir="$(mktemp -d)"
PYTHONPYCACHEPREFIX="$cache_dir" python -m compileall -q src
rm -rf "$cache_dir"
```

Run a packaging/import smoke check after installation:

```console
python -c "import pymatmc2; import mexm"
```

The recovered tests can be explored with:

```console
python -m pytest --collect-only
```

The legacy test tree is not yet guaranteed to collect or pass as a complete
suite. Start with focused tests for the code being changed. Never run a test
that launches an external calculator or submits a scheduler job without
explicit authorization.

Build distribution metadata with:

```console
python -m build
```

Do not publish the resulting artifacts unless publication is explicitly
authorized. Remove local `dist/`, `build/`, and `*.egg-info` outputs after
inspection.

## Adding new MC² mechanisms

For each proposed mechanism:

1. state the ensemble, state variables, invariants, and acceptance rule;
2. identify how mass balance, phase fractions, and detailed balance are
   maintained;
3. distinguish structural relaxation from Monte Carlo sampling;
4. add deterministic software tests before external calculations;
5. define numerical-verification evidence separately from scientific
   validation; and
6. document assumptions and limitations.

Do not infer scientific correctness from historical behavior or from a passing
unit test.

## External calculations

The source contains interfaces for VASP and historical scheduler workflows.
Those interfaces do not confer permission to run them. Before any authorized
calculation, record the executable, input system, numerical settings, expected
outputs, estimated resources, and data-retention plan.
