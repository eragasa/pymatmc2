import pytest

import numpy as np
from numpy import linalg
from typing import List
from typing import Tuple

from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import ConcentrationMatrix
from pymatmc2.error import Pymatmc2ConcentrationMatrixError 
# DEFINE SOME MODULE LEVEL VARIABLES
# these variables ensure that the same variables used in development routines
# and testing routines are the same

concentration_matrix_ = np.array(
    [
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4]
    ]
)

total_concentration_ = np.array(
    [1/4, 1/4, 1/4]
)

# DEFINE TESTING FIXTURES
@pytest.fixture
def configuration_path():
    configuration_path_ = 'pymatmc2.config'
    return configuration_path_

@pytest.fixture
def configuration(configuration_path_):
    configuration_ = Pymatmc2Configuration()
    configuration_.read(path=configuration_path)

@pytest.fixture
def o_concentration_matrix():
    o = ConcentrationMatrix()
    return o

@pytest.fixture
def concentration_matrix():
    A = concentration_matrix_
    return A

# DEFINE TESTS
def test__input_matrices_satisfy_testing_condictions(concentration_matrix):
    all_sum_to_unity, phase_sums_to_unity = ConcentrationMatrix.cell_concentrations_sum_to_unity(A=concentration_matrix)
    assert all_sum_to_unity
    assert linalg.matrix_rank(concentration_matrix) < min(concentration_matrix.shape)

def test__is_rank_deficient(concentration_matrix):
    assert ConcentrationMatrix.is_rank_deficient(A=concentration_matrix_)

  
# DEFINE DEVELOPMENT ROUTINES
def dev__input_matrices_satisfy_testing_condictions():
    all_sum_to_unity, phase_sums_to_unity = ConcentrationMatrix.cell_concentrations_sum_to_unity(A=concentration_matrix_)
    print(all_sum_to_unity)
    print(phase_sums_to_unity)

def dev__is_rank_deficient():
    is_rank_deficient = ConcentrationMatrix(A=concentration_matrix_)

def dev__staticmethod__SVD():
    A = concentration_matrix_
    ConcentrationMatrix.SVD(A=A, is_debug=True)


def dev__do_SVD__no_args():
    A = concentration_matrix_
    o = ConcentrationMatrix()
    o.do_SVD()
    
def dev__do_SVD__w_A():
    A = concentration_matrix_

    o = ConcentrationMatrix()
    o.do_SVD(A=A)

def dev__SVD():
    is_debug = True
    # Ax = b
    A = cell_concentration_matrix_
    b = total_concentration_

    if is_debug:
        print('A:\n{}'.format(A))

    U, S, Vt = linalg.svd(A, full_matrices=True)
    if is_debug:
        print('U:\n{}'.format(U))
        print('S:\n{}'.format(S))
        print('Vt\n:{}'.format(Vt))

    # reconstruction test
    m, n = A.shape
    if m == n:
        USV = U @ np.diag(S) @ Vt
    else:
        USV = U[:,:n] @ np.diag(S) @ Vt[:m,:]
    passes_reconstruction_test = np.allclose(A, USV)

    if is_debug:    
        print('USV:\n{}'.format(USV))
        print('A - USV:\n{}'.format(A-USV))
        print('passes_reconstruction_test:{}'.format(passes_reconstruction_test))

    
    # determine singular solutions
    tolerance = np.max(np.size(A)) * np.spacing(np.max(np.diag(S)))
    p = np.sum(S > tolerance)
    # reduced space
    Up = np.matrix(U[:, :p])
    Vp = np.matrix(Vt[:p, :])
    Sp = np.diag(S)[:p, :p]
    if is_debug:
        print('tolerance:{}'.format(tolerance))
        print('p:{}'.format(p))
        print('Up:{}'.format(Up))
        print('Vp:{}'.format(Vp))
        print('S_diag:{}'.format(Sp))

    # singular solutions reconstruction test   
    USVp = Up @ Sp @ Vp
    print('USVp:{}'.format(USVp)) 
    print('Vp.shape:{}'.format(Vp.shape))
    print('Up.shape:{}'.format(Up.shape))

    # invert using SVD
    print('invert using SVD')
    SpInv_diag = np.diag(1/S)[:p,:p]
    print('SpInv_diag:{}'.format(SpInv_diag))
    print('SpInv_diag.shape:{}'.format(SpInv_diag.shape))
    AInv_p = Vp.T @ SpInv_diag @ Up.T
    print('AInv_p:{}'.format(AInv_p))
    print('AInt_p:{}'.format(AInv_p))
    print(linalg.eig(AInv_p)[0])
    print(linalg.eig(AInv_p)[1])
    b = b.reshape(m, 1)
    phase_molar_fraction = np.matrix(AInv_p) @ np.matrix(b)
    print('phase_molar_fraction={}'.format(phase_molar_fraction))



if __name__ == "__main__":
    dev__input_matrices_satisfy_testing_condictions()
    dev__is_rank_deficient()
    dev__staticmethod__SVD()
    dev__do_SVD__no_args()
