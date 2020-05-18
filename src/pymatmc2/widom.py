""" This class implements the Widom Test """
import os
from copy import deepcopy
from typing import Tuple
import numpy as np
from mexm.structure import SimulationCell
from mexm.io.vasp import Poscar

from pymatmc2 import constants
from pymatmc2 import Pymatmc2Configuration

class WidomTest(object):

    def __init__(self):
        self._configuration = None
        self._cell_initial = None
        self._n_iterations = None
        self._symbol1 = None
        self._symbol2 = None
        self._widom_results_path = None


    @property
    def configuration(self) -> Pymatmc2Configuration:
        return self._configuration

    @configuration.setter
    def configuration(self, configuration: Pymatmc2Configuration):
        self._configuration = configuration

    @property
    def cell_initial(self) -> SimulationCell:
        return self._cell_initial

    @cell_initial.setter
    def cell_initial(self, cell: SimulationCell):
        self._cell_initial = cell

    @property
    def n_iterations(self) -> int:
        return self._n_iterations

    @n_iterations.setter
    def n_iterations(self, n: int):
        if not isinstance(n, int):
            msg = 'n_iterations must be integer type'
            raise TypeError(msg)

        self._n_iterations = n

    @property
    def symbol1(self) -> str:
        return self._symbol1

    @symbol1.setter
    def symbol1(self, symbol):
        if not isinstance(symbol, str):
            msg = 'symbol1 must be a string type'
            raise TypeError(msg)
        self._symbol1 = symbol

    @property
    def symbol2(self) -> str:
        return self._symbol2

    @symbol2.setter
    def symbol2(self, symbol):
        if not isinstance(symbol, str):
            msg = 'symbol2 must be a string type'
            raise TypeError(msg)
        self._symbol2 = symbol

    def generate_structure_files(self, n_iterations: int):
        self.n_iterations = n_iterations

        assert issubclass(self.cell_initial, SimulationCell)
        assert self.symbol1 != self.symbol2
        assert self.symbol1 in self.cell_initial.symbols
        assert self.symbol2 in self.configuration.symbols
        assert self.phase in self.configuration.cell_names

        for i_iteration in range(self.n_iterations + 1):

            if i_iteration == 0:
                self.create_widom_simulation_path()
                self.create_widom_results_path()
                cell = deepcopy(self.cell_initial)
            else:
                cell = deepcopy(self.cell_initial)
                idxs_symbol1 = []
                for i, atom in enumerate(self.cell_initial.atomic_basis):
                    if atom.symbol == self.symbol1:
                        idxs_symbol1.append(i)
                idx = np.random.choice(idxs_symbol1, 1)[0]
                # create new cell
                cell = deepcopy(self.cell_initial)
                cell.atomic_basis[idx].symbol = self.symbol2

            calculator_type = 'vasp'

            simulation_generators = {}
            simulation_generators['vasp'] = self.generate_mutated_cell_w_calculator_vasp

            simulation = simulation_generators[calculator_type](cell)
            cell_generators[calculator_type](cell=self.initial_cell)

            # write new cell
            if isinstance(cell):
                simulation_path = self.configuration.simulation_path
                phase = 'fcc1'
                cell_name = '{}_{}_{:05}.vasp'.format(self.symbol1, self.symbol2, i_iteration)
                cell_name.write(path=cell_name)

    def generate_simulation_w_calculator_vasp(self, cell: Poscar) -> VaspSimulation:

        incar_path = ""
        kpoints_path = ""
        o_vasp = VaspSimulation()
        o_vasp.incar = Incar()
        o_vasp.poscar = cell
        o_vasp.kpoints = Kpoints()

        return o_vasp

    def create_simulation(self, simulation_path: str, phase: 'str', cell: SimulationCell):

        if phase not in self.configuration.cell_names:
            msg = 'cannot identify phase, {}, in the configuration file'
            raise ValueError(msg)






    def get_widom_simulation_results(self):
        

    def get_widom_statistic(self, n_iterations: int):
        N_iterations = self.configuration.widom_n_iterations
        for i in range(N_iterations):
            pass