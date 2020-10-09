import pytest
import os
import platform
import numpy as np
from numpy import linalg
from typing import List
from typing import Tuple
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import ConcentrationMatrix

concentration_matrix_ = np.array(
    [
        [1/2, 1/2], 
        [1/2, 1/2]
    ]
)

def get_resource_path():
    resource_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__))
    )

def get_configuration_path():
    configuration_path_dict = {
        'Windows':'pymatmc2.config.windows',
        'Mac':'pymatmc2.config.osx',
        'Linux':'pymatmc2.config.linux'
    }
    os_type = platform.system()
    configuration_path_ = os.path.join(
        get_resource_path(),
        configuration_path_dict[os_type]
    )
    return configuration_path_

# DEFINE TESTING FIXTURES
@pytest.fixture
def resource_path():
    resource_path_ = os.path.join(
        os.path.dirname(os.path.abspath(__file__))
    )
    return resource_path_

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
    assert o.X is None
    assert o.U is None
    assert o.S is None
    assert o.Vt is None
    assert o.XInv is None
    assert o.SVD_method is None

def test__constructor__w_configuration_as_Pymatmc2Configuration(configuration):
    o = ConcentrationMatrix(configuration=configuration)
    assert isinstance(o.configuration, Pymatmc2Configuration)
    assert o.X is None
    assert o.U is None
    assert o.S is None
    assert o.Vt is None
    assert o.XInv is None
    assert o.SVD_method is None

def test__constructor__w_configuration_as_invalid_class():
    bad_configuration = {}
    with pytest.raises(TypeError) as e:
        o = ConcentrationMatrix(configuration=bad_configuration)
    expected_error_msg = "configuration needs to be an instance of Pymatmc2Configuration"
    assert str(e.value) == expected_error_msg

def test__constructor__w_X_as_numpy_array(concentration_matrix):
    X = concentration_matrix
    o = ConcentrationMatrix(X=X)
    assert isinstance(o.X, np.ndarray)
    assert o.U is None
    assert o.S is None
    assert o.Vt is None
    assert o.XInv is None
    assert o.SVD_method is None

def test__constructor__w_X_as_list(concentration_matrix):
    X = concentration_matrix.tolist()
    assert isinstance(X, list)

    o = ConcentrationMatrix(X=X)
    assert isinstance(o.X, np.ndarray)
    assert o.U is None
    assert o.S is None
    assert o.Vt is None
    assert o.XInv is None
    assert o.SVD_method is None

def test__property__configuration(configuration):
    o = ConcentrationMatrix()
    o.configuration = configuration
    assert isinstance(o.configuration, Pymatmc2Configuration)
