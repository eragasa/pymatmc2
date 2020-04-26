""" This class implements the Widom Test """
from copy import deepcopy
from typing import Tuple
import numpy as np
from mexm.structure import SimulationCell
from pymatmc2 import Pymatmc2Configuration
from pymatmc2.mutator import CellFlip
class WidomTest(object):

    def __init__(self):
        self._configuration = None

    @property
    def configuration(self) -> Pymatmc2Configuration:
        return self._configuration

    @configuration.setter
    def configuration(self, configuration: Pymatmc2Configuration):
        self._configuration = configuration

    def get_widom_test_cells(self):
        N_iterations = self.configuration.widom_n_iterations
        for i in range(N_iterations):
            

    def get_widom_simulation_results(self):
        pass

    def get_widom_statistic(self):
        N_iterations = self.configuration.widom_n_iterations
        for i in range(N_iterations):

class WidomTestCell(self):

    def __init__(self):
        self._configuration = None

    @property
    def configuration(self) -> Pymatmc2Configuration:
        return self._configuration

    @configuration.setter
    def configuration(self, configuration: Pymatmc2Configuration):
        self._configuration = configuration

    def mutate_cell(self, cell: SimulationCell) -> Tuple[str, str, SimulationCell]:
        new_cell = deepcopy(cell)

        # 1. select an atom
        
        n_atoms = cell.get_n_atoms('all')
        idx_atm = np.choose(range(n_atoms), 1)[0]

        
        # 2. choose symbol to flip to
        s1 = cell.atomic_basis[idx_atm].symbol
        symbols = deepcopy(configuration.symbols)
        symbols.remove(s1)
        s2 = np.choose(range(symbols), 1)[0]

        # 3. flip the symbol
        new_cell.atomic_basis[idx_atm].symbol = s2

        return s1, s2, new_cell