# Historical MPI-spawn reference

This directory preserves the unsuccessful historical
`MultiCellMonteCarloMpiSpawn` experiment as a reference for later analysis. It
is not package source, is not imported by `pymatmc2`, and is not represented as
working software.

The historical file was created in `pymatmc2` commit
`5aa620857fbdcde5f8d76ea65711d356e82c667e` and cites the abICS
`run_base_mpi.py` implementation. Comparison with the contemporaneous abICS
revision shows substantial shared structure and implementation text.

- abICS source commit: `3023806c3cdc613a38f135a9278d54bd08942b7e`
- Upstream file:
  `abics/applications/latgas_abinitio_interface/run_base_mpi.py`
- Upstream license: GPL-3.0-or-later
- Historical `pymatmc2` file SHA-256:
  `8d827d53a2c8bebcf4739dedf665f28bd87364f453c929f5f63b661046dfeda3`

The reference copy adds a corrective GPL notice before the historical content.
Its full license is retained in `LICENSE`. The original bytes remain preserved
in the separate Project Koios legacy archive.
