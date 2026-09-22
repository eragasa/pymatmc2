import pytest
import numpy as np
from numpy import linalg
from typing import List
from typing import Tuple
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import ConcentrationMatrix
from pymatmc2.error import Pymatmc2ConcentrationMatrixError

concentration_matrix_ = np.array(
    [
        [1/2, 1/2], 
        [1/2, 1/2]
    ]
)

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


@pytest.fixture
def concentration_matrix() -> np.ndarray:
    return concentration_matrix_

def test__gettter__not_set():
    o = ConcentrationMatrix()
    assert o.A is None

def test__setter__as_numpy_array(concentration_matrix):
    A = concentration_matrix
    assert isinstance(A, np.ndarray)

    o = ConcentrationMatrix()
    o.A = A
    assert isinstance(o.A, np.ndarray)

def test__setter__as_list(concentration_matrix):
    A = concentration_matrix.tolist()
    assert isinstance(A, list)

    o = ConcentrationMatrix()
    o.A = A
    assert isinstance(o.A, np.ndarray)
    assert np.allclose(o.A, concentration_matrix)

def test__setter__w_bad_concentration_matrix():
    A = case__bad_concentration_matrix['A']

    o = ConcentrationMatrix()

    # testing that correct exception gets raised
    with pytest.raises(Pymatmc2ConcentrationMatrixError) as e:
        o.A = A

    # checking attributes
    try:
        o.A = A
    except Pymatmc2ConcentrationMatrixError as e:
        assert e.kwargs['phase_sums_to_unity'] == case__bad_concentration_matrix['phase_sums_to_unity']

def dev__setter__w_bad_concentration_matrix():
    A = case__bad_concentration_matrix['A']   
    o = ConcentrationMatrix()
    try:
        o.A = A
    except Pymatmc2ConcentrationMatrixError as e:
        print(type(e))
        print(e.kwargs['phase_sums_to_unity'])
        print(e.kwargs['phase_sums_to_unity'] == case__bad_concentration_matrix['phase_sums_to_unity'])

if __name__ == "__main__":
    dev__setter__w_bad_concentration_matrix()
