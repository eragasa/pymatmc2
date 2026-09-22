# Scientific background

The Multi-Cell Monte Carlo method, written MC², represents possible coexisting
solid phases in separate simulation cells. This avoids explicitly modeling the
interfaces between bulk phases. Composition changes and phase fractions are
coupled so that the overall composition remains constrained.

The three primary papers associated with development of the method are:

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

The 2017 work introduced a multi-cell relaxation approach with fixed cell sizes.
The 2019 work introduced lever-rule coupling to permit variable phase fractions
and applied the method with density-functional-theory energies. The 2020 work
developed more general acceptance criteria and a predictor-corrector consistency
check associated with chemical-potential equilibrium, and demonstrated an
approach using classical interatomic potentials.

A later overview is:

- M. Ghazisaeidi, “Alloy thermodynamics via the Multi-cell Monte Carlo (MC)2
  method,” *Computational Materials Science* **193**, 110322 (2021).
  <https://doi.org/10.1016/j.commatsci.2021.110322>

## Relationship to this repository

`pymatmc2` is Eugene J. Ragasa's Python implementation motivated by this body of
work. The repository is not the source archive for the papers and does not imply
endorsement by their authors.

The current code is a recovered historical implementation. Statements in the
papers are not automatically capabilities of this software. Each mechanism,
ensemble, and acceptance rule must be traced to source behavior and verified
before it is described as implemented.
