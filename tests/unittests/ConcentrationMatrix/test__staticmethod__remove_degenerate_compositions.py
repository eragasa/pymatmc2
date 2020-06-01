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

def test__4x4__rank1():
    X0 = case__4by4__rank1['X']
    X1 = ConcentrationMatrix.remove_degenerate_compositions(X=X0, is_debug=True)
    assert isinstance(X1, np.ndarray)
    assert X1.shape[0] == case__4by4__rank1['X1.shape'][0]
    assert X1.shape[1] == case__4by4__rank1['X1.shape'][1]

def dev__remove_degenerate_compositions():
    X0 = case__4by4__rank1['X']
    X1 = ConcentrationMatrix.remove_degenerate_compositions(X=X0, is_debug=True)
    print(X1.shape)

if __name__ == "__main__":
    dev__remove_degenerate_compositions()