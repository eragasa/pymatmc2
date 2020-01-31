from collections import OrderedDict
from mexm.io.vasp import Poscar

class SimulationCells:

    def __init__(self, n_cells):
        """

        Arguments:
            n_cells (int): number of cells
        """
        assert isinstance(n_cells, int)
        self.n_cells = n_cells
        self.cells = [Poscar() for i in range(self.n_cells)]

    @property
    def symbols(self):
        symbols = []
        for cell in self.cells:
            for symbol in cell.symbols:
                if symbol not in symbols:
                    symbols.append(symbol)
        return symbols

    def get_number_of_atoms(self, symbol=None):
        n_atoms = 0
        for cell in self.cells:
            n_atoms += cell.get_number_of_atoms(symbol)
        return n_atoms
