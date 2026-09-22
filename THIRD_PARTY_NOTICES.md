# Third-party notices

## Bundled mexm component

The historical `mexm` package under `src/mexm/` originated as the separate
Materials Ex Machina project and is bundled because `pymatmc2` imports it
directly.

- Upstream project: <https://github.com/eragasa/mexm-base>
- Upstream commit: `3963a1a30ba1595fba69cd0f791dc31fe2e08f7c`
- Copyright: Copyright (c) 2019 Eugene Ragasa
- License: MIT
- Exact license text: `MEXM_LICENSE`
- Detailed provenance: `MEXM_PROVENANCE.md`

## abICS-derived historical reference

`references/mc2/multicellmc_spawn.py` is retained only as a non-operational
historical reference. It derives substantially from
`abics/applications/latgas_abinitio_interface/run_base_mpi.py` in abICS commit
`3023806c3cdc613a38f135a9278d54bd08942b7e`.

- Upstream project: <https://github.com/issp-center-dev/abICS>
- Upstream copyright: Copyright (C) 2019–, The University of Tokyo
- License: GNU General Public License, version 3 or later
- Exact license text: `references/mc2/LICENSE`

The reference is not imported by or included in the MIT-licensed `pymatmc2`
package. Its former historical MIT header does not relicense the abICS-derived
portions.
