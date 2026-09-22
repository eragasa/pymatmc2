import pytest
import numpy as np
from numpy import linalg
from pymatmc2 import ConcentrationMatrix
from pymatmc2 import PhaseMolarFraction

case__4by4__rank2 = {}
case__4by4__rank2['X'] = np.array(
    [
        [3/4,   3/4, 1/12, 1/12],
        [1/12, 1/12,  3/4,  3/4],
        [1/12, 1/12, 1/12, 1/12],
        [1/12, 1/12, 1/12, 1/12]
    ]
)
case__4by4__rank2['c'] = np.array([0.25, 0.25, 0.25, 0.25])
case__4by4__rank2['E'] = np.array([ -543.09136824, -543.29211709, -540.62577658, -543.60646897])

def test__calculate():
    is_debug=True
    X0 = case__4by4__rank2['X']
    c = case__4by4__rank2['c']
    E = case__4by4__rank2['E']

    f = PhaseMolarFraction.calculate(X=X0, c=c, E=E, is_debug=is_debug)
    assert isinstance(f, np.ndarray)

def dev__calculate():
    is_debug = True
    X0 = case__4by4__rank2['X']
    c = case__4by4__rank2['c']
    E = case__4by4__rank2['E']

    f = PhaseMolarFraction.calculate(X=X0, c=c, E=E, is_debug=is_debug)

if __name__ == "__main__":
    dev__calculate()
