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

# DEFINE TESTING FIXTURES

@pytest.fixture
def configuration_path():
    configuration_path_ = 'pymatmc2.config'
    return configuration_path_

@pytest.fixture
def concentration_matrix() -> np.ndarray:
    return concentration_matrix_

@pytest.fixture
def configuration(configuration_path):
    configuration_ = Pymatmc2Configuration()
    configuration_.read(path=configuration_path)
    return configuration_

@pytest.fixture
def o_concentration_matrix():
    o = ConcentrationMatrix()
    return o

# DEFINE TESTS
def test__constructor__no_args():
    o = ConcentrationMatrix()
    assert o.configuration is None
    assert o.A is None
    assert o.U is None
    assert o.S is None
    assert o.Vt is None
    assert o.AInv is None
    assert o.SVD_method is None

def test__constructor__w_configuration_as_Pymatmc2Configuration(configuration):
    o = ConcentrationMatrix(configuration=configuration)
    assert isinstance(o.configuration, Pymatmc2Configuration)
    assert o.A is None
    assert o.U is None
    assert o.S is None
    assert o.Vt is None
    assert o.AInv is None
    assert o.SVD_method is None

def test__constructor__w_configuration_as_invalid_class():
    bad_configuration = {}
    with pytest.raises(TypeError) as e:
        o = ConcentrationMatrix(configuration=bad_configuration)
    expected_error_msg = "configuration needs to be an instance of Pymatmc2Configuration"
    assert str(e.value) == expected_error_msg

def test__constructor__w_A_as_numpy_array(concentration_matrix):
    A = concentration_matrix
    o = ConcentrationMatrix(A=A)
    assert isinstance(o.A, np.ndarray)
    assert o.U is None
    assert o.S is None
    assert o.Vt is None
    assert o.AInv is None
    assert o.SVD_method is None

def test__constructor__w_A_as_list(concentration_matrix):
    A = concentration_matrix.tolist()
    assert isinstance(A, list)

    o = ConcentrationMatrix(A=A)
    assert isinstance(o.A, np.ndarray)
    assert o.U is None
    assert o.S is None
    assert o.Vt is None
    assert o.AInv is None
    assert o.SVD_method is None

def test__property__configuration(configuration):
    o = ConcentrationMatrix()
    o.configuration = configuration
    assert isinstance(o.configuration, Pymatmc2Configuration)
