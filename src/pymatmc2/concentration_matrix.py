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
        A: Optional[np.ndarray]=None,
        configuration: Optional[Pymatmc2Configuration]=None
    ):
        self._configuration = None
        self._A = None
        self._U = None
        self._S = None
        self._Vt = None
        self._AInv = None
        self._SVD_method = None

        if configuration is not None: 
            self.configuration = configuration
        if A is not None: 
            self.A = A

    @property
    def configuration(self):
        return self._configuration

    @configuration.setter
    def configuration(self, configuration):
        if not isinstance(configuration, Pymatmc2Configuration):
            raise TypeError("configuration needs to be an instance of Pymatmc2Configuration")
        self._configuration = configuration


    @property
    def A(self) -> np.ndarray:
        return self._A

    @A.setter
    def A(self, A: np.ndarray):
        if isinstance(A, np.ndarray):
            self._A = A
        elif isinstance(A, list):
            self._A = np.array(A)
        else:
            raise TypeError('cannot convert argument into a numpy array')

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
    def AInv(self):
        return self._AInv

    def _set_AInv(self, AInv: np.ndarray):
        if isinstance(AInv, np.ndarray):
            self._AInv = AInv
        elif isinstance(AInv, list):
            self._AInv = np.array(AInv)
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
    def cell_concentrations_sum_to_unity(A: np.ndarray, is_debug = False) -> Tuple[bool, List[bool]]:
        """
        Arguments:
            A (numpy.ndarray): concentration matrix
        Returns:
            Tuple[bool, List[bool]
        """

        m, n = A.shape
        # m is the number of rows, correspodning with the number of elements
        # n is the number of columns, corresponding with the number of simulation cells

        column_sums = A.sum(axis=0)
        if is_debug: 
            print(column_sums)

        sums_to_unity = []
        for column_sum in column_sums.tolist():
            sums_to_unity.append(column_sum == 1.)
        
        if is_debug:
            print(sums_to_unity)
        return all(sums_to_unity), sums_to_unity

    @staticmethod
    def is_rank_deficient(A: np.ndarray, is_debug: Optional[bool] = False) -> bool:
        """

        A matrix is rank deficient if it does not have full rank.  A matrix is full rank
        if its rank equals the largest possible for a matrix lesser of the numbers of rows and columns.

        Arguments:
            A (numpy.ndarray): concentration matrix
        Returns:
            Tuple[bool, List[
        """
        m, n = A.shape 
        # m is the number of rows, correspodning with the number of elements
        # n is the number of columns, corresponding with the number of simulation cells

        if linalg.matrix_rank(A) < min(m, n):
           is_rank_deficient_ = True
        else:
           is_rank_deficient_ = False

        return is_rank_deficient_

    @staticmethod
    def SVD(A: np.ndarray, is_debug = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if is_debug:
            print('A:\n{}'.format(A))

        U, S, Vt = linalg.svd(A, full_matrices=True)
        if is_debug:
            print('U:\n{}'.format(U))
            print('S:\n{}'.format(S))
            print('Vt:\n{}'.format(Vt))

        return U, S, Vt

    @staticmethod
    def SVD_reconstruction_test(
        A: np.ndarray,
        U: np.ndarray,
        S: np.ndarray,
        Vt: np.ndarray,
        is_debug=False
    ) -> Tuple[bool, np.ndarray]:

        m, n = A.shape
        if m == n:
            USV = U @ np.diag(S) @ Vt
        else:
            USV =  U[:, :n] @ np.diag(S) @ Vt[:m, :]

        if is_debug:
            print('A:\n{}'.format(A))
            print('A - USV:\n{}'.format(A-USV))

        return np.allclose(A, USV), A - USV

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

