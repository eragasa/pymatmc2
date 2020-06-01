import pytest
import numpy as np
from numpy import linalg
from pymatmc2 import ConcentrationMatrix
case__4by4__rank1 = {}
# concentration matrix
case__4by4__rank1['X'] = np.array(
    [
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4]
    ]
)
case__4by4__rank1['X1.shape'] = [4,1]
case__4by4__rank1['c'] = np.array(
    # total concentration vector
    [0.25, 0.25, 0.25, 0.25]
)
case__4by4__rank1['E'] = np.array(
    # cell energies
    [ -543.09136824, -543.29211709, -540.62577658, -543.60646897]
)

def test__4x4__rank1__SVD():
    X0 = case__4by4__rank1['X']
    U, S, Vt = ConcentrationMatrix.SVD(X0, is_debug=True)
    assert isinstance(U, np.ndarray)
    assert isinstance(S, np.ndarray)
    assert isinstance(Vt, np.ndarray)
    assert np.allclose(U @ np.diag(S) @ Vt, X0)

def test__4x4__rank1__SVD_w_X1():
    X0 = case__4by4__rank1['X']
    X1 = ConcentrationMatrix.remove_degenerate_compositions(X=X0, is_debug=True)
    assert isinstance(X1, np.ndarray)
    assert X1.shape[0] == case__4by4__rank1['X1.shape'][0]
    assert X1.shape[1] == case__4by4__rank1['X1.shape'][1]

    U, S, Vt = ConcentrationMatrix.SVD(X1, is_debug=True)
    assert isinstance(U, np.ndarray)
    assert isinstance(S, np.ndarray)
    assert isinstance(Vt, np.ndarray)
