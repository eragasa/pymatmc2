import numpy as np
from numpy import linalg
from pymatmc2 import ConcentrationMatrix
from pymatmc2 import PhaseMolarFraction

case__4by4__rank3 = {}
case__4by4__rank3['X'] = np.array(
    [
        [0.25925926, 0.24074074, 0.24074074, 0.25925926],
        [0.25925926, 0.24074074, 0.25925926, 0.24074074],
        [0.24074074, 0.25925926, 0.24074074, 0.25925926],
        [0.24074074, 0.25925926, 0.25925926, 0.24074074]
    ]
)
case__4by4__rank3['c'] = np.array(
    [0.25, 0.25, 0.25, 0.25]
)
case__4by4__rank3['E'] = np.array(
    [ -543.09136824, -543.29211709, -540.62577658, -543.60646897]
)

is_debug = True
X0 = case__4by4__rank3['X']
c = case__4by4__rank3['c']
X1 = ConcentrationMatrix.remove_degenerate_compositions(X=X0, is_debug=True)
X1_matrix_rank = linalg.matrix_rank(X1)
print('X1_matrix_rank:{}'.format(X1_matrix_rank))
U, S, Vt = ConcentrationMatrix.SVD(X1, is_debug=True)
U, S, Vt = ConcentrationMatrix.reduced_SVD(X1, is_debug=True)

# invert using SVD
m, n = X1.shape
S_inv = np.diag(1/S)
C_inv = Vt.T @ S_inv @ U.T
c = c.reshape(m, 1)

f = np.matrix(C_inv) @ np.matrix(c)
if is_debug:
    print(80*'-')
    print('svd_solution:\n{}'.format(f))

phase_molar_fraction_is_valid = PhaseMolarFraction.is_valid(f=f)
print('phase_molar_fraction_is_valid:{}'.format(phase_molar_fraction_is_valid))

#PhaseMolarFraction.calculate(A=A, b=b, E=E, is_debug=True)
