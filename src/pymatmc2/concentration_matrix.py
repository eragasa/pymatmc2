from copy import deepcopy
import numpy as np
from numpy import linalg
from typing import List
from typing import Tuple
from typing import Optional

from mexm.exception import MexmException
from pymatmc2 import Pymatmc2Configuration
from pymatmc2.error import Pymatmc2ConcentrationMatrixError

class ConcentrationMatrix():

    supported_svd_methods = ['default']

    def __init__(
        self, 
        X: Optional[np.ndarray]=None,
        configuration: Optional[Pymatmc2Configuration]=None
    ):
        self._configuration = None
        self._X = None
        self._U = None
        self._S = None
        self._Vt = None
        self._XInv = None
        self._SVD_method = None

        if configuration is not None: 
            self.configuration = configuration
        if X is not None: 
            self.X = X

    @property
    def configuration(self):
        return self._configuration

    @configuration.setter
    def configuration(self, configuration):
        if not isinstance(configuration, Pymatmc2Configuration):
            raise TypeError("configuration needs to be an instance of Pymatmc2Configuration")
        self._configuration = configuration


    @property
    def X(self) -> np.ndarray:
        return self._X

    @X.setter
    def X(self, X: np.ndarray):
        if isinstance(X, np.ndarray):
            X_ = deepcopy(X)
        elif isinstance(X, list):
            X_ = np.array(X)
        else:
            raise TypeError('cannot convert argument into a numpy array')

        all_sum_to_unity, phase_sums_to_unity = ConcentrationMatrix.cell_concentrations_sum_to_unity(X=X_)
        if not all_sum_to_unity:
            msg = "The columns of the concentration matrix do not sum to unity"
            kwargs = {
                'phase_sums_to_unity':phase_sums_to_unity
            }
            raise Pymatmc2ConcentrationMatrixError(msg, **kwargs)
        
        # assign to private member variable
        self._X = X_

    @property
    def U(self) -> np.ndarray:
        return self._U

    def _set_U(self, U):
        if isinstance(U, np.ndarray):
            self._U = U
        elif isinstance(U, list):
            self._U = np.array(U)
        else:
            raise TypeError('cannot convert arguement into a numpy array')

    @property
    def S(self) -> np.ndarray:
        return self._S

    def _set_U(self, S: np.ndarray):
        if isinstance(S, np.ndarray):
            self._S = S
        elif isinstance(S, list):
            self._S = np.array(S)
        else:
            raise TypeError('cannot convert argument into a numpy array')

    @property
    def Vt(self):
        return self._Vt

    def _set_Vt(self, Vt: np.ndarray):
        if isinstance(Vt, np.ndarray):
            self._Vt = Vt
        elif isinstance(Vt, list):
            self._Vt = np.array(Vt)
        else:
            raise TypeError('cannot convert argument into a numpy array')

    @property
    def XInv(self):
        return self._XInv

    def _set_XInv(self, XInv: np.ndarray):
        if isinstance(XInv, np.ndarray):
            self._XInv = XInv
        elif isinstance(AInv, list):
            self._XInv = np.array(XInv)
        else:
            raise TypeError('cannot onvert argument into a numpy array')

    @property
    def SVD_method(self) -> str:
        return self._SVD_method

    def _set_SVD_method(self, s: str):
        if not isinstance(s, str):
            raise TypeError('SVD must be a string type')

        if s not in self.supported_svd_method:
            raise ValueError('{} is not a supported SVD type'.format(s))

        self._SVD_method = SVD_method

    @staticmethod
    def cell_concentrations_sum_to_unity(X: np.ndarray, is_debug = False) -> Tuple[bool, List[bool]]:
        """
        Arguments:
            A (numpy.ndarray): concentration matrix
        Returns:
            Tuple[bool, List[bool]
        """

        m, n = X.shape
        # m is the number of rows, correspodning with the number of elements
        # n is the number of columns, corresponding with the number of simulation cells

        column_sums = X.sum(axis=0)
        if is_debug: 
            print(column_sums)

        sums_to_unity = []
        for column_sum in column_sums.tolist():
            sums_to_unity.append(column_sum == 1.)
        
        if is_debug:
            print(sums_to_unity)
        return all(sums_to_unity), sums_to_unity

    @staticmethod
    def is_rank_deficient(X: np.ndarray, is_debug: Optional[bool] = False) -> bool:
        """

        A matrix is rank deficient if it does not have full rank.  A matrix is full rank
        if its rank equals the largest possible for a matrix lesser of the numbers of rows and columns.

        Arguments:
            X (numpy.ndarray): concentration matrix
        Returns:
            Tuple[bool, List[
        """
        m, n = X.shape 
        # m is the number of rows, correspodning with the number of elements
        # n is the number of columns, corresponding with the number of simulation cells

        if linalg.matrix_rank(X) < min(m, n):
           is_rank_deficient_ = True
        else:
           is_rank_deficient_ = False

        return is_rank_deficient_



    @staticmethod
    def SVD(
        X: np.ndarray, 
        is_debug: Optional[bool] = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Arguments:
            X (np.ndarray): The matrix to be decomposed.
            is_debug (Optional[bool]):  
        Returns:
            (Tuple[np.ndarray, np.ndarray, np.ndarray]): U S Vt portions of the 
                SVD decomposition of a matrix.
        """

        if is_debug:
            print('A:\n{}'.format(X))

        U, S, Vt = linalg.svd(X, full_matrices=True)

        m, n = X.shape
        if m == n:
            S = np.diag(S)
            USVt = U @ S @ Vt
        else:
            U = U[:, :n]
            Vt = Vt[:m, :]
            sigma_shape = [U.shape[1], Vt.shape[0]]
            sigma = np.zeros(sigma_shape)
            for i, v in enumerate(S.tolist()):
                sigma[i,i] = v

            USVt =  U[:, :n] @ sigma @ Vt[:m, :]

        if is_debug:
            print('U:\n{}'.format(U))
            print('S:\n{}'.format(S))
            print('Vt:\n{}'.format(Vt))
            print('USVt:\n{}'.format(USVt))
            print('reconstruction_passed:{}'.format(np.allclose(USVt, X)))
        return U, S, Vt

    @staticmethod
    def SVD_reconstruction_test(
        X: np.ndarray,
        U: np.ndarray,
        S: np.ndarray,
        Vt: np.ndarray,
        is_debug=False
    ) -> Tuple[bool, np.ndarray]:

        m, n = X.shape
        if m == n:
            USVt = U @ np.diag(S) @ Vt
        else:
            USVt =  U[:, :n] @ np.diag(S) @ Vt[:m, :]

        is_passed = np.allclose(X, USVt)
        if is_debug:
            print('RECONSTRUCTION TEST')
            print('X:\n{}'.format(X))
            print('X - USVt:\n{}'.format(X-USVt))
            print('is_passed:{}'.format(is_passed))
        return is_passed, X - USVt

    @staticmethod
    def invert_using_SVD_components(U: np.ndarray, S: np.ndarray, Vt: np.ndarray) -> np.ndarray:
        S_inv = np.diag(1/S)
        X_inv = Vt.T @ S_inv @ U.T
        return X_inv



    @staticmethod
    def reduced_SVD(X: np.ndarray,is_debug: Optional[bool] = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        X (np.ndarray): cell concentration matrix
        is_debug (bool): Default value is false

        Returns:
        ========
        (np.ndarray): U
        (np.ndarray): S
        (np.ndarray): V
        """
        assert isinstance(X, np.ndarray)
        if is_debug:
            print('X:\n{}'.format(X))

        if is_debug:
            print(80*'=')
            print('reduced space SVD')
            print(80*'=')
            print('X:\n{}'.format(X))

        U, S, Vt = linalg.svd(X, full_matrices=True)
        if is_debug:
            print(80*'-')
            print('DECOMPOSITION OF X')
            print('U:\n{}'.format(U))
            print('S:\n{}'.format(S))
            print('Vt\n:{}'.format(Vt))
    
          # determine singular solutions
        tolerance = np.max(np.size(X)) * np.spacing(np.max(np.diag(S)))
        p = np.sum(S > tolerance)
        if is_debug:
            print('tolerance:{}'.format(tolerance))
            print('p:{}'.format(p))

        # reduced space
        Up = np.matrix(U[:, :p])
        Vp = np.matrix(Vt[:p, :])
        Sp = S[:p]
        Sp_diag = np.diag(S)[:p, :p]
        if is_debug:
            print('Up:\n{}'.format(Up))
            print('Vp:\n{}'.format(Vp))
            print('Sp:\n{}'.format(Sp))
   
        USVp = Up @ Sp_diag @ Vp
        if is_debug:
            print('USVp:\n{}'.format(USVp)) 
            print('Vp.shape:{}'.format(Vp.shape))
            print('Up.shape:{}'.format(Up.shape))
            print('Sp.shape:{}'.format(Sp.shape))

        return Up, Sp, Vp

    def do_SVD(self, A: Optional[np.ndarray]=None):
        if A is not None:
            self.A = A

        U, S, Vt = SVD(A)
        reconstruction_test_results = ConcentrationMatrix.SVD_reconstruction_test(A, U, Vt)
        if reconstruction_test_results[0]:
            self._U = U
            self._S = S
            self._Vt = Vt
        else:
            msg = "SVD fails reconstruction matrix test"
            kwargs = {
                'U': U,
                'S': S,
                'Vt': Vt,
                'A_less_USV': reconstruction_test_results[1]
            }
            raise Pymatmc2ConcentrationMatrixError(m, **kwargs)

        return U, S, Vt

    def remove_degenerate_compositions(X: np.ndarray, is_debug: Optional[bool]=False) -> np.ndarray:
        if isinstance(X, np.ndarray):
            X0 = X
        elif isinstance(X, list):
            X0 = np.array(X)
        else:
            raise TypeError('unable to cast X into numpy array')

        X1 = np.unique(X0, axis=1)
        if is_debug:
            print('remove_degenerate_compositions()')
            print('X0:\n{}'.format(X0))
            print('X1:\n{}'.format(X1))
        return X1