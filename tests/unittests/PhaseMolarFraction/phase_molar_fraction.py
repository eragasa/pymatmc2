import numpy as np
from numpy import linalg
from typing import Optional

class PhaseMolarFraction():

    @staticmethod
    def sums_to_unity(f: np.ndarray) -> bool:
        if np.sum(f) == 1:
            sums_to_unity = True
        else:
            sums_to_unity = False
        return sums_to_unity

    @staticmethod
    def has_no_negative_components(f: np.ndarray) -> bool:
        return not np.any(f < 0)

    @staticmethod
    def is_valid(f):
        tests = [
            PhaseMolarFraction.sums_to_unity,
            PhaseMolarFraction.has_no_negative_components
        ]

        return all([test(f=f) for test in tests])

    @staticmethod
    def calculate_from_svd(A: np.ndarray, b: np.ndarray, is_debug: Optional[bool] = False):
        assert isinstance(A, np.ndarray)
        assert isinstance(b, np.ndarray)
        if is_debug:
            print('A:\n{}'.format(A))
            print('b:\n{}'.format(b))

        U, S, Vt = linalg.svd(A, full_matrices=True)
        if is_debug:
            print(80*'-')
            print('DECOMPOSITION OF A')
            print('U:\n{}'.format(U))
            print('S:\n{}'.format(S))
            print('Vt\n:{}'.format(Vt))

        # test decomposition A
        m, n = A.shape
        if m == n:
            USV = U @ np.diag(S) @ Vt
        else:
            USV = U[:,:n] @ np.diag(S) @ Vt[:m,:]
        is_USV_equals_A = np.allclose(A, USV)

        if is_debug:
            print(80*'-')
            print('TEST DECOMPOSITION OF A')
            print('USV:\n{}'.format(USV))
            print('A - USV:\n{}'.format(A-USV))
            print('is_USV_equals_A:{}'.format(is_USV_equals_A))

        # invert using SVD
        S_inv = np.diag(1/S)
        A_inv = Vt.T @ S_inv @ U.T
        b = b.reshape(m, 1)

        x = np.matrix(A_inv) @ np.matrix(b)
        if is_debug:
            print(80*'-')
            print('svd_solution:\n{}'.format(x))

        return x

    @staticmethod
    def calculate_from_reduced_svd(A: np.ndarray, b: np.ndarray, is_debug: Optional[bool] = False):
        assert isinstance(A, np.ndarray)
        assert isinstance(b, np.ndarray)
        if is_debug:
            print(80*'=')
            print('Solve using reduced SVD')
            print(80*'=')
            print('A:\n{}'.format(A))
            print('b:\n{}'.format(b))

        U, S, Vt = linalg.svd(A, full_matrices=True)
        if is_debug:
            print(80*'-')
            print('DECOMPOSITION OF A')
            print('U:\n{}'.format(U))
            print('S:\n{}'.format(S))
            print('Vt\n:{}'.format(Vt))
    
  
        # determine singular solutions
        tolerance = np.max(np.size(A)) * np.spacing(np.max(np.diag(S)))
        p = np.sum(S > tolerance)
        if is_debug:
            print('tolerance:{}'.format(tolerance))
            print('p:{}'.format(p))

        # reduced space
        Up = np.matrix(U[:, :p])
        Vp = np.matrix(Vt[:p, :])
        Sp = np.diag(S)[:p, :p]
        if is_debug:
            print('Up:{}'.format(Up))
            print('Vp:{}'.format(Vp))
            print('Sp:{}'.format(Sp))
   
        USVp = Up @ Sp @ Vp
        if is_debug:
            print('USVp:{}'.format(USVp)) 
            print('Vp.shape:{}'.format(Vp.shape))
            print('Up.shape:{}'.format(Up.shape))

        # invert using SVD
        Sp_inv = np.diag(1/S)[:p, :p]
        Ap_inv = Vp.T @ Sp_inv @ Up.T
        if is_debug:
           print('invert using SVD')
           print('SpInv_diag:{}'.format(Sp_inv))
           print('AInv_p:{}'.format(Ap_inv))

        m, n = A.shape
        b = b.reshape(m, 1)
        x = np.matrix(Ap_inv) @ np.matrix(b)
        if is_debug:
            print(80*'-')
            print('reduced_svd_solution:\n{}'.format(x))

        return x

    @staticmethod
    def calculate(A: np.ndarray, b: np.ndarray, E: np.ndarray, is_debug: Optional[bool]=False):
        m, n = A.shape
        if linalg.matrix_rank(A)  ==min(m, n):
            x = PhaseMolarFraction.calculate_from_svd(A=A, b=b, is_debug=is_debug)
        else:
            x = PhaseMolarFraction.calculate_from_reduced_svd(A=A, b=b, is_debug=is_debug)
        
        A0 = A
        A1 = np.unique(A, axis=1)
        if is_debug:
            print('A0:\n{}'.format(A0))
            print('A1:\n{}'.format(A1))

        m0, n0 = A0.shape
        m1, n1 = A1.shape

        if n0 != n1:
            if is_debug:
               print('duplicate compositions!')
               print('using lowest energy phase for degenerate composition')
               print('E:\n{}'.format(E))
            for idx_n1 in range(n1):
               lowest_energy = None
               for idx_n0 in range(n0):
                   if np.allclose(A0[:,idx_n0], A1[:,idx_n1]):
                       if lowest_energy is None:
                           lowest_energy = E[idx_n0]
                           idx_lowest_energy = idx_n0
                       elif E[idx_n0] < lowest_energy:
                           lowest_energy = E[idx_n0]
                           idx_lowest_energy = idx_n0
                       else:
                           pass
               sum_phase_fraction = 0
               for idx_n0 in range(n0):
                   if np.allclose(A0[:,idx_n0], A1[:, idx_n1]):
                       sum_phase_fraction += x[idx_n0]
                       x[idx_n0] = 0
               x[idx_lowest_energy] = sum_phase_fraction
        if is_debug:
            print('X1:{}'.format(x))         
        return x
