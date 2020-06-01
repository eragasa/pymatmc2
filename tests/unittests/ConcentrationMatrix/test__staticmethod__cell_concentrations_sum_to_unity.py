import pytest

import numpy as np
from numpy import linalg
from typing import List
from typing import Tuple

from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import ConcentrationMatrix
from pymatmc2.error import Pymatmc2ConcentrationMatrixError

case__good_concentration_matrix = {}
case__good_concentration_matrix['A'] = np.array(
    [
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4]
    ]
)
case__good_concentration_matrix['all_sum_to_unity'] = True
case__good_concentration_matrix['phase_sums_to_unity'] = [True, True, True, True]

case__bad_concentration_matrix = {}
case__bad_concentration_matrix['A'] = np.array(
    [
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4],
        [1/5, 1/4, 1/4, 1/4]
    ]
)
case__bad_concentration_matrix['all_sum_to_unity']= False
case__bad_concentration_matrix['phase_sums_to_unity'] = [False, True, True, True]

def test__w_good_concentration_matrix():
    A = case__good_concentration_matrix['A']
    expected_all_sum_to_unity = case__good_concentration_matrix['all_sum_to_unity']
    expected_phase_sums_to_unity = case__good_concentration_matrix['phase_sums_to_unity']

    all_sum_to_unity, phase_sums_to_unity = ConcentrationMatrix.cell_concentrations_sum_to_unity(A=A)

    assert all_sum_to_unity == expected_all_sum_to_unity
    for v, expected_v in zip(phase_sums_to_unity, expected_phase_sums_to_unity):
        assert v == expected_v

def test__w_bad_concentration_matrix():
    A = case__bad_concentration_matrix['A']
    expected_all_sum_to_unity = case__bad_concentration_matrix['all_sum_to_unity']
    expected_phase_sums_to_unity = case__bad_concentration_matrix['phase_sums_to_unity']

    all_sum_to_unity, phase_sums_to_unity = ConcentrationMatrix.cell_concentrations_sum_to_unity(A=A)

    assert all_sum_to_unity == expected_all_sum_to_unity
    for v, expected_v in zip(phase_sums_to_unity, expected_phase_sums_to_unity):
        assert v == expected_v
