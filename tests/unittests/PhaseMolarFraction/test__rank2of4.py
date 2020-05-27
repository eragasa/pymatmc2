import pytest
import numpy as np
from numpy import linalg
from phase_molar_fraction import PhaseMolarFraction
from pymatmc2.cell_concentration_matrix import CellConcentrationMatrix

expected_concentration_matrix = \
    [
        [3/4,   3/4, 1/12, 1/12],
        [1/12, 1/12,  3/4, 3/4],
        [1/12, 1/12, 1/12, 1/12],
        [1/12, 1/12, 1/12, 1/12]
    ]
expected_total_concentration = [0.25, 0.25, 0.25, 0.25]
expected_cell_energies = [ -543.09136824, -543.29211709, -540.62577658, -543.60646897]


def test__inputs_are_valid():

is_debug = True

A = np.array(expected_concentration_matrix)
E = np.array(expected_cell_energies)
b = np.array(expected_total_concentration)

PhaseMolarFraction.calculate(A=A, b=b, E=E, is_debug=True)
