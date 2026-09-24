# Repository structure

```text
pymatmc2/
├── src/
│   ├── pymatmc2/       Active recovered MC² implementation
│   └── mexm/           Bundled historical dependency
├── tests/              Recovered source-oriented tests and resources
├── docs/               Maintained project documentation
├── doc/                Historical Sphinx configuration fragment
├── references/         Non-operational historical material
├── sbin/               Recovered command-line scripts
├── devtools/           Historical offline-artifact selection policy
├── pyproject.toml      Installation and dependency metadata
├── OFFLINE_ARTIFACTS.md  Private preservation evidence and limits
├── PROVENANCE.md       pymatmc2 recovery boundary
└── MEXM_PROVENANCE.md  mexm recovery boundary
```

## Active implementation

`src/pymatmc2/` owns the MC² workflow, concentration and phase-fraction
representations, mutation behavior, results handling, and orchestration logic.
Future MC² development should normally occur here.

## Bundled mexm dependency

`src/mexm/` supplies the historical structure, VASP-I/O, simulation, filesystem,
and job-submission objects imported by `pymatmc2`. It is bundled to preserve the
original source relationship. It is not a separately modernized dependency.

Changes to `src/mexm/` should be narrowly justified. The original identities in
`MEXM_SOURCE_SHA256SUMS` are immutable provenance records, not checksums to be
regenerated after future edits.

## Historical references

Material under `references/` is retained for provenance and investigation, not
normal imports:

- `broken/` contains incomplete `pymatmc2` source;
- `master-initial/` contains the compact initial implementation;
- `mc2/` contains a GPL-licensed, non-operational abICS-derived experiment; and
- `mexm-broken/` contains incomplete historical `mexm` source.

The GPL reference must remain separated from the MIT-licensed active package.

## Tests and resources

The recovered `tests/` tree mixes unit tests, development scripts, and authored
input generators. It should be modernized incrementally. Tests added for new
work must be deterministic and must not launch external calculators or submit
scheduler jobs.

Large generated calculation trees, pseudopotentials, wavefunctions, restart
files, and databases do not belong in this repository. Selected historical
bytes remain privately recoverable under the policy and limitations documented
in `OFFLINE_ARTIFACTS.md`; they are never package or test inputs.
