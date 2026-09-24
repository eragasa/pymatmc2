# AGENTS.md

## Project

`pymatmc2` is research software for the Multi-Cell Monte Carlo (MC²) method for
alloy phase stability and coexistence. The repository currently combines a
recovered historical implementation with preparation for renewed development.
Treat recovered behavior and new work as provisional unless supported by
explicit tests and scientific evidence.

## Authoritative locations

| Subject | Location |
|---|---|
| Active MC² implementation | `src/pymatmc2/` |
| Bundled historical dependency | `src/mexm/` |
| Installation metadata | `pyproject.toml` |
| Maintained documentation | `README.md` and `docs/` |
| pymatmc2 recovery provenance | `PROVENANCE.md` |
| mexm recovery provenance | `MEXM_PROVENANCE.md` |
| Historical and excluded material | `references/` |
| Offline artifact policy and evidence | `devtools/offline-artifacts.toml` and `OFFLINE_ARTIFACTS.md` |
| Historical Sphinx fragment | `doc/` |

Read the relevant provenance record before changing recovered source or moving
files across these boundaries.

## Preservation rules

- `SOURCE_SHA256SUMS` and `SOURCE_GIT_MODES` record the original recovered
  `pymatmc2` snapshot. Do not rewrite them to disguise later source changes.
- `MEXM_SOURCE_SHA256SUMS` and `MEXM_SOURCE_GIT_MODES` serve the same purpose
  for the bundled `mexm` snapshot.
- Material beneath `references/` is non-operational historical evidence. Do not
  import it into active packages or silently repair it in place.
- `OFFLINE_ARTIFACT_SHA256SUMS` records private historical artifact identities.
  Do not restore those bytes into the maintained tree or use them as package or
  test inputs.
- Keep the GPL-licensed `references/mc2/` material separate from the
  MIT-licensed active packages.
- Do not add pseudopotentials, wavefunctions, calculator outputs, restart data,
  databases, generated results, credentials, or scheduler secrets to Git.

A checksum mismatch after an intentional source edit is not fixed by changing
the historical checksum. Document the new change and preserve the historical
identity record.

## Scientific and execution boundaries

- Distinguish expected behavior, software verification, numerical verification,
  and scientific validation.
- Do not describe passing tests as scientific validation.
- Do not claim convergence, phase stability, or agreement with a publication
  without retained inputs, settings, outputs, and an identified comparison.
- Do not execute VASP, LAMMPS, GULP, PhonTS, MPI jobs, scheduler submissions, or
  other external calculations without explicit human authorization.
- Never add or redistribute a `POTCAR` unless its redistribution rights and
  exact intended use have been established by the human owner.
- New Monte Carlo mechanisms are planned work until implemented, tested, and
  documented. Do not present roadmap items as current capabilities.

## Python development

- Support installation through `pyproject.toml`; do not duplicate dependency
  authority in `requirements.txt`.
- Keep `pymatmc2` and the bundled `mexm` package importable from the shared
  `src/` layout.
- Prefer new implementation work in `src/pymatmc2/`. Modify bundled `mexm`
  only when the dependency boundary genuinely requires it, and document the
  departure from the preserved snapshot.
- Add focused tests for changed behavior. Tests must not launch external
  calculators, scheduler jobs, or network services.
- Use compact authored fixtures and temporary runtime directories. Do not add
  generated simulation trees as fixtures.
- Do not weaken assertions, tolerances, or test selection merely to obtain a
  pass.

## Documentation

- Keep installation and development instructions in `docs/` and link them from
  `README.md`.
- Label historical behavior, current behavior, and proposed work distinctly.
- Cite the primary MC² papers when describing the scientific method.
- Keep provenance details in the provenance records instead of duplicating the
  complete recovery narrative across every page.

## Working procedure

1. Inspect the current branch and working tree.
2. Read the relevant source, tests, documentation, and provenance records.
3. Make the smallest coherent change.
4. Run the cheapest relevant checks first, then broader checks where feasible.
5. Report changed paths, commands run, warnings, and remaining limitations.
6. Commit or push only when explicitly requested.
