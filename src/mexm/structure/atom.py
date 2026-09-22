import numpy as np
from mexm.util import is_number

class Atom(object):
    """description of an atom

    This position is a data structure which contains information about an
    individual atom

    Args:
        symbol (str): the standard ISO symbol for an element
        position (list of float): the position of the atom the units
           are dependent upon use.
        magmom (float): the magnetic moment of the atom.

    Attributes:
        symbol (str): the standard ISO symbol for an element
        position (numpy.ndarray): the position of the atom, usually in direct
            coordinates
        magnetic_moment (float): the magnetic moment of the atom
    """
    def __init__(self, symbol, position, magmom=0, atom_id=None):
        assert isinstance(symbol, str)

        self.symbol = symbol
        self.position = None
        self.magnetic_moment = None
        self.atom_id = None

        self._initialize_position(position=position)
        self._initialize_magnetic_moment(magmom=magmom)
        self._initialize_atom_id(atom_id=None)

    def _initialize_position(self, position):
        assert isinstance(position, list) or isinstance(position, np.ndarray)

        if isinstance(position, list):
            self.position = np.array(position)
        else:
            self.position = position.copy()

    def _initialize_magnetic_moment(self, magmom):
        assert is_number(magmom)
        self.magnetic_moment = float(magmom)

    def _initialize_atom_id(self, atom_id):
        if atom_id is None:
            self.atom_id = None
        elif isinstance(atom_id, str):
            self.atom_id = atom_id
        else:
            self.atom_id = str(atom_id)

    @property
    def magmom(self):
        return self.magnetic_moment

    @magmom.setter
    def magmom(self, magmom):
        self.magnetic_moment = magmom
