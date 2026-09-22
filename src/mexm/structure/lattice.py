from copy import deepcopy
import numpy as np

class Lattice(object):
    unit_lattice = [[1,0,0],[0,1,0],[0,0,1]]
    primitive_lattice = [[1,0,0],[0,1,0],[0,0,1]]
    """ a crystal lattice

    Args:
        a0 (float): the length of the a1 lattice vector
        H (list of list): a 3x3 matrix containing the lattice vectors as ranks
    """

    def __init__(self, a0=1.0, H='unit'):
        self.a0 = a0
        self.H_ = None

        self.initialize_H_matrix(H=H)

    @staticmethod
    def initialize_from_dict(obj_dict):
        o = Lattice(
            a0 = obj_dict['a0'],
            H = obj_dict['H']
        )
        return o

    def initialize_H_matrix(self, H):
        if isinstance(H, str):
            if H == 'unit':
                self.H = np.array(self.unit_lattice, dtype=float)
            elif H == 'primitive':
                self.H = np.array(self.primitive_lattice, dtype=float)
            else:
                raise ValueError

        elif isinstance(H, list):
            assert len(H) == 3
            assert all([len(k) == 3 for k in H])
            self.H = np.array(H)

        elif isinstance(H, np.ndarray):
            self.H = np.copy(H)

        else:
            raise TypeError

    @property
    def volume(self):
        return np.dot(self.a1, np.cross(self.a2, self.a3))

    @property
    def reciprocal_volume(self):
        return np.dot(self.b1, np.cross(self.b2, self.b3))

    @property
    def H(self):
        return self.H_

    @H.setter
    def H(self, Hmatrix):
        if isinstance(Hmatrix, list):
            self.H_ = np.array(Hmatrix)
        elif isinstance(Hmatrix, np.ndarray):
            self.H_ = np.copy(Hmatrix)

    @property
    def h1(self):
        return self.H[:,0]

    @h1.setter
    def h1(self, h1):
        if isinstance(h1, list):
            self.H[:,0] = np.array(h1, dtype=float)
        else:
            self.H[:,0] = h1

    @property
    def h2(self):
        return self.H[:,1]

    @h2.setter
    def h2(self, h2):
        if isinstance(h2, list):
            self.H[:,1] = np.array(h2, dtype=float)
        else:
            self.H[:,1] = h2
    @property
    def h3(self):
        return self.H[:,2]

    @h3.setter
    def h3(self, h3):
        if isinstance(h3, list):
            self.H[:,2] = np.array(h3, dtype=float)
        else:
            self.H[:,2] = h3

    @property
    def a1(self):
        return self.a0 * self.h1

    @property
    def a2(self):
        return self.a0 * self.h2

    @property
    def a3(self):
        return self.a0 * self.h3

    @property
    def a1_length(self):
        """float: length of the a1 lattice vector"""
        a0 = self.a0
        h1 = self.H[:,0]
        a1_length = a0 * np.dot(h1,h1)**2
        return a1_length

    @property
    def a2_length(self):
        """float: length of the a1 lattice vector"""
        a0 = self.a0
        h2 = self.H[:,1]
        a2_length = a0 * np.dot(h2,h2)**2
        return a2_length

    @property
    def a3_length(self):
        """float: length of the a1 lattice vector"""
        a0 = self.a0
        h3 = self.H[:,2]
        a3_length = a0 * np.dot(h3,h3)**2
        return a3_length

    @property
    def b1(self):
        """numpy.array: this is a 3x1 numpy array, in reciprocal space"""
        a1 = self.a1
        a2 = self.a2
        a3 = self.a3

        b1 = 2 * np.pi * np.cross(a2,a3) / np.dot(a1, np.cross(a2,a3))
        return b1

    @property
    def b2(self):
        """numpy.array: this is a 3x1 numpy array, in reciprocal space"""
        a1 = self.a1
        a2 = self.a2
        a3 = self.a3

        b2 = 2 * np.pi * np.cross(a3,a1) / np.dot(a2, np.cross(a3,a1))
        return b2

    @property
    def b3(self):
        """numpy.array: this is a 3x1 numpy array, in reciprocal space"""
        a1 = self.a1
        a2 = self.a2
        a3 = self.a3

        b3 = 2 * np.pi * np.cross(a1,a2) / np.dot(a3, np.cross(a1,a2))
        return b3

    def to_dict(self):
        return_dict = {
            'a0':self.a0,
            'H':self.H
        }

        return return_dict

    def __str__(self):
        return_rows = [
            "a0:{:10.6f}".format(self.a0),
            "H-MATRIX",
            "\th1:{:10.6f} {:10.6f} {:10.6f}".format(*self.h1),
            "\th2:{:10.6f} {:10.6f} {:10.6f}".format(*self.h2),
            "\th3:{:10.6f} {:10.6f} {:10.6f}".format(*self.h3),
            "REAL SPACE LATTICE VECTORS",
            "\ta1:{:10.6f} {:10.6f} {:10.6f}".format(*self.a1),
            "\ta2:{:10.6f} {:10.6f} {:10.6f}".format(*self.a2),
            "\ta3:{:10.6f} {:10.6f} {:10.6f}".format(*self.a3),
            "RECIPROCAL SPACE LATTICE VECTORS",
            "\tb1:{:10.6f} {:10.6f} {:10.6f}".format(*self.b1),
            "\tb2:{:10.6f} {:10.6f} {:10.6f}".format(*self.b2),
            "\tb3:{:10.6f} {:10.6f} {:10.6f}".format(*self.b3)
        ]
        return_str = "\n".join(return_rows)

        return return_str
