# Release-maintenance policy

`offline-artifacts.toml` is the reviewed pymatmc2 selection policy for generated
and restricted files retained in historical commit
`3063f031301c04e395cb9d7bc5064f2a52e1110c`. The maintained v0.1 source tree
already excludes these paths; do not apply this policy to the current checkout.

The generic staging, verification, archive, and recovery implementation is
pinned to `projectkoios-bootstrap` commit
`c0f406685a49838a690c1519a6174564bbaa057b` (tree
`3be0101c8a582e9cdcaa943197176770d52c7f78`). Apply it only to a clean,
disposable checkout detached at the historical commit. The checkout must ignore
`.offline/` locally.

```bash
PYTHONPATH=/explicit/projectkoios-bootstrap/python python3.14 -m \
  projectkoios.bootstrap.harness.offline_artifacts stage \
  --repository-root /disposable/historical/checkout \
  --policy /maintained/pymatmc2/devtools/offline-artifacts.toml \
  --destination .offline/release-readiness \
  --checksum-manifest OFFLINE_ARTIFACT_SHA256SUMS
```

Do not use `--remove-originals` in the disposable historical checkout. Archive
creation independently verifies every selected byte against the exact Git
object. Do not execute retained Python, calculators, schedulers, shell scripts,
or recovered workflow content.
