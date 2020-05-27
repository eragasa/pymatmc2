import numpy as np
from numpy import linalg
from pymatmc2.phase_molar_fraction import PhaseMolarFraction

expected_concentration_matrix = \
    [
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4]
    ]
expected_total_concentration = [0.25, 0.25, 0.25, 0.25]
expected_cell_energies = [ -543.09136824, -543.29211709, -540.62577658, -543.60646897]

is_debug = True

A = np.array(expected_concentration_matrix)
E = np.array(expected_cell_energies)
b = np.array(expected_total_concentration)

PhaseMolarFraction.calculate(A=A, b=b, E=E, is_debug=True)
