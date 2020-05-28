import pytest
import numpy as np
from numpy import linalg
from typing import List
from typing import Tuple
from pymatmc2 import Pymatmc2Configuration
from pymatmc2.concentration_matrix import ConcentrationMatrix

concentration_matrix_ = np.array(
    [
        [1/2, 1/2], 
        [1/2], [1/2]
    ]
)

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
    #assert np.allclose(o.A, concentration_matrix)

