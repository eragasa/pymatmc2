import numpy as np
from numpy import linalg
from pymatmc2 import ConcentrationMatrix
from pymatmc2 import PhaseMolarFraction

case__4by4__rank1 = {}
case__4by4__rank1['X'] = np.array(
    # concentration matrix
    [
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4],
        [1/4, 1/4, 1/4, 1/4]
    ]
)
case__4by4__rank1['c'] = np.array(
    # total concentration vector
    [0.25, 0.25, 0.25, 0.25]
)
case__4by4__rank1['E'] = np.array(
    # cell energies
    [ -543.09136824, -543.29211709, -540.62577658, -543.60646897]
)

def test__calculate():
    is_debug=True
    X0 = case__4by4__rank1['X']
    c = case__4by4__rank1['c']
    E = case__4by4__rank1['E']

    f = PhaseMolarFraction.calculate(X=X0, c=c, E=E, is_debug=is_debug)
    assert isinstance(f, np.ndarray)

def dev__calculate():
    is_debug = True
    X0 = case__4by4__rank1['X']
    c = case__4by4__rank1['c']
    E = case__4by4__rank1['E']

    f = PhaseMolarFraction.calculate(X=X0, c=c, E=E, is_debug=is_debug)

if __name__ == "__main__":
    dev__calculate()
