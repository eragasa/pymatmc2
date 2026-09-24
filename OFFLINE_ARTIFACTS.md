# Offline artifact preservation

The maintained v0.1 tree is source-only. Selected generated calculator output,
restricted VASP pseudopotentials, workflow state, results, and generated
documentation from historical commit
`3063f031301c04e395cb9d7bc5064f2a52e1110c` remain privately recoverable
without restoring them to the maintained source tree.

The repository-owned historical selection policy is
`devtools/offline-artifacts.toml`. Generic staging, verification, deterministic
archive, and recovery mechanics are provided by `projectkoios-bootstrap` commit
`c0f406685a49838a690c1519a6174564bbaa057b`.

The policy selected 3,191 files totaling 3,166,798,699 bytes. Every staged byte
was verified against the historical Git commit and its exact commit tree
`fd9788560365e93faae596ac91f2530aa6bd92a8`. A complete archive recovery into a
new directory reverified all paths, sizes, checksums, and Git identities.

The private preservation bundle is stored locally under:

```text
~/Library/CloudStorage/Dropbox/pymatmc2/v0.1/data/
```

Its archive is:

```text
pymatmc2-offline-artifacts-git-3063f031301c.tar
SHA-256 c26ba00feff668323b5d02aab12ed96b6b9270eb2d1677f52b581c7ebb2df0c5
```

`OFFLINE_ARTIFACT_SHA256SUMS` records conventional per-file checksums without
including private bytes. The bundle contains its canonical manifest, source
identity, archive checksum, and recovery instructions. Dropbox remote
synchronization and an independent second backup must be confirmed separately.

The archive is private historical evidence, not a package input, public release
asset, or reproducibility claim. Do not execute recovered content. In
particular, no redistribution authority is claimed for the 30 retained
`POTCAR` files.
