import os
import shutil
from copy import deepcopy
from collections import OrderedDict
from typing import List, Dict
import numpy as np
from numpy import linalg
from mexm.io.vasp import Poscar
from mexm.simulation import VaspSimulation

from pymatmc2 import Pymatmc2Configuration

class MultiCellError(BaseException): pass

class MultiCell:
    """

    Attributes:
        configuration (Pymatmc2Configuration)
        molar_fraction_total (List[float])
        simulations (Dict[str, Simulation])
    """

    def __init__(self):
        """
        """
        self.configuration = None
        self._molar_fraction_total = None
        self.simulations = None

    @staticmethod
    def initialize_from_obj(multicell):
        """

        Args:
            multicell (MultiCell)
        Returns:
            MultiCell:
        """
        obj = MultiCell()
        obj.configuration = deepcopy(multicell.configuration)
        obj.simulations = deepcopy(multicell.simulations)
        return obj

    @staticmethod
    def initialize_from_pymatmc2_configuration(
        configuration: Pymatmc2Configuration
    ):
        """
            configuration (Pymatmc2Configuration)
        """
        obj = MultiCell()
        obj.configuration = configuration

        simulation_cells = configuration.simulation_cells
        order_simulation_cells = OrderedDict(
            sorted(simulation_cells.items())
        )
        obj.simulations = OrderedDict()
        if configuration.calculator_type == 'vasp':
            for k, v in order_simulation_cells.items():
                obj.simulations[k] = VaspSimulation()
                obj.simulations[k].poscar.read(path=v['poscar'])
                obj.simulations[k].incar.read(path=v['incar'])
                obj.simulations[k].potcar.read(path=v['potcar'])
                obj.simulations[k].kpoints.read(path=v['kpoints']) 

        return obj
    
    @property
    def total_energy(self):
        # see Niu, et al "Multi-cell Monte Carlo method for phase prediction",
        # npj Computational Materials, Eq (3) with P=0.

        # number of phases
        m = len(self.simulations)
        
        # 
        sum_U = 0
        for phase in self.simulations:
            U = self.simulations[phase].total_energy
            f = self.phase_molar_fraction[phase]
            sum_U += U * f

        return m * sum_U

    @property
    def concentration(self) -> Dict[str, float]:

        assert isinstance(self.configuration, Pymatmc2Configuration)

        sum_concentration = sum(self.configuration.concentration.values())
        
        concentration = OrderedDict()
        for s in self.symbols:
            concentration_symbol =  self.configuration.concentration[s]
            concentration[s] = concentration_symbol/sum_concentration

        return concentration
    
    @property
    def cell_concentration(self) -> Dict[str, float]:

        cell_concentration = OrderedDict()
        for c in self.cell_names:
            cell_concentration[c] = OrderedDict()

        for c in self.cell_names:
            for s in self.symbols:
                n_atoms_symbol = self.cells[c].get_number_of_atoms(s)
                n_atoms_total = self.cells[c].get_number_of_atoms()
                cell_concentration[c][s] = n_atoms_symbol/n_atoms_total

        return cell_concentration
    
    @property
    def cell_concentration_matrix(self) -> np.ndarray:

        X = []
        for c in self.cell_names:
            X.append(list(self.cell_concentration[c].values()))

        return np.array(X)

    @property
    def molar_fraction_total(self) -> Dict[str, float]:
        assert isinstance(self.configuration, Pymatmc2Configuration)

        molar_fraction_total = self.configuration.molar_fraction_total
        molar_fraction_total = OrderedDict()

        sum_molar_fraction_total = sum(self.configuration.molar_fraction_total.values())
        for k, v in self.configuration.molar_fraction_total.items():
            molar_fraction_total[k] = v/sum_molar_fraction_total

        return molar_fraction_total



    def get_simulation_paths(self, path):
        return [os.path.join(path, k) for k in self.simulations]

    def read(self, path: str):
        """ read simulations from disk

        Arguments:
            path (str)
        """

        assert os.path.isdir(path)
        self.simulations = OrderedDict()
        
        for cell_name in self.configuration.simulation_cells:
            if self.configuration.calculator_type == 'vasp':
                vasp_simulation_path = os.path.join(path, cell_name)
                self.simulations[cell_name] = VaspSimulation()
                self.simulations[cell_name].read(vasp_simulation_path)


    def write(self, path: str):
        """ write simulations to disk

        Writes out the simulations to path.  Currently only tested with VASP simulations, but
        should work with any subclass of the mexm.simulation.Simulation.

        Arguments:
            path (str): the path of the directory which to write the simulations
        """
        if os.path.isdir(path):
            shutil.rmtree(path)
        
        os.mkdir(path)
        for simulation_name, simulation_obj in self.simulations.items():
            simulation_path = os.path.join(path, simulation_name)
            os.mkdir(simulation_path)
            simulation_obj.write(simulation_path=simulation_path)

    def configure(self, configuration: Pymatmc2Configuration):
        """ configure class from a Pymatmc2Configuration
        
        this method sets up the cells attribute

        Arguments:
            configuration (Pymatmc2Configuration): instance of the 
                Pymatmc2Configuration in which to setup this class.
         """

        assert isinstance(configuration, Pymatmc2Configuration)
        # set argument to attribute value
        self.configuration = configuration

        # configure simulations
        ordered_simulation_cells = OrderedDict(
            sorted(self.configuration.simulation_cells.items())
        )
        self.simulations = OrderedDict()
        if self.configuration.calculator_type == 'vasp':
            for k, v in ordered_simulation_cells.items():
                self.simulations[k] = VaspSimulation()
                self.simulations[k].poscar.read(path=v['poscar'])
                self.simulations[k].incar.read(path=v['incar'])
                self.simulations[k].potcar.read(path=v['potcar'])
                self.simulations[k].kpoints.read(path=v['kpoints'])

        
        # configure molar fraction total
        sum_molar_fraction_total \
            = sum(self.configuration.molar_fraction_total.values())
        self.molar_fraction_total = OrderedDict()
        for k, v in self.configuration.molar_fraction_total.items():
            self.molar_fraction_total[k] = v/sum_molar_fraction_total
            

    def add_cells_from_dict(self, simulation_cells):
        """ add cells from dictionary object

        simulation_cells = {
            '<cell_name>' = {'<file_type'>:'<path_to_file>'}
        }
        Arguments:
            simulation_cells (dict) = dictionary of simulation cells
        """

    @property
    def cell_names(self):
        ordered_cell_names = sorted([k for k in self.simulations])
        return ordered_cell_names
    
    @property
    def cells(self):
        cells = OrderedDict()
        for k in self.cell_names:
            if isinstance(self.simulations[k], VaspSimulation):
                cells[k] = self.simulations[k].poscar
        return cells

    @property
    def cell_molar_fraction(self) -> List[List[float]]:
        X = OrderedDict()
        for k in self.cell_names:

            X[k] = OrderedDict()
            sum_n_atoms = self.cells[k].n_atoms
            for s in self.symbols:
                n_atoms = self.cells[k].get_number_of_atoms(s)
                X[k][s] = n_atoms/sum_n_atoms
        return X
    
    @property
    def symbols(self):
        symbols = []
        for cell in self.cells.values():
            for s in cell.symbols:
                if s not in symbols:
                    symbols.append(s)
        symbols.sort()
        return symbols

    @property
    def total_molar_fraction(self):
        c = OrderedDict(
            [(s, self.molar_fraction_total[s]) for s in self.symbols]
        )
        return c
    
    @property
    def phase_molar_fraction(self):

        # transform cell molar fraction in a matrix
        X = []
        for _, v in self.cell_molar_fraction.items():
            X.append([v[s]for s in self.symbols])
        X = np.array(X)

        # transform total molar fraction into a column vector
        c = [v for _, v in self.total_molar_fraction.items()]            
        c = np.array(c)

        # solve for phase molar fraction
        f = np.dot(linalg.inv(X), c)
        for k in f:
            if k < 0:
                # msg = "phase molar fraction cannot be negative\n"
                # msg += "X:{}\n".format(X)
                # msg += "c:{}\n".format(c)
                # msg += "f:{}\n".format(f)
                # raise MultiCellError(msg)
                from scipy.optimize import nnls
                f, residuals = nnls(X,c)
                break 
        
        return {v:f[i] for i,v in enumerate(self.cell_names)}
        
    def get_number_of_atoms(self, symbol=None):
        n_atoms = 0
        for cell in self.cells.values():
            n_atoms += cell.get_number_of_atoms(symbol)
        return n_atoms
