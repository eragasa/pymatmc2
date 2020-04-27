import os
import shutil
from copy import deepcopy
from collections import OrderedDict
from typing import List, Dict

import numpy as np
from numpy import linalg

from mexm.structure import SimulationCell
from mexm.io.vasp import Poscar
from mexm.simulation import Simulation
from mexm.simulation import VaspSimulation

from pymatmc2 import Pymatmc2Configuration

class MultiCellError(BaseException): pass

class MultiCell:
    """

    Attributes:
        src_path (str)
        dst_path (str)
        configuration (Pymatmc2Configuration)
        molar_fraction_total (List[float])
        simulations (Dict[str, Simulation])
    """

    def __init__(self):
        """
        """
        self.src_path = None
        self.dst_path = None
        self.configuration = None
        self._molar_fraction_total = None
        self._simulations = None

    @property
    def simulations(self) -> Dict[str, Simulation]:
        return self._simulations

    @simulations.setter
    def simulations(self, simulations: Dict[str, Simulation]):
        for k, v in simulations.items():
            assert isinstance(k, str)
            assert issubclass(v, Simulation)
        self._simulations = simulations

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
    def initialize_from_configuration(
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
            U_per_atom = U /self.simulations[phase].n_atoms
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

        X = np.array(X)
        return X.T

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
        Raises:
            FileNotFoundError: if path directory does not exist.
        """

        if not os.path.isdir(path):
            msg = 'cannot read read the directory:{}'.format(path)
            raise FileNotFoundError(msg)
        
        self.src_path = path
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

    def archive(self, dst_path: str):

        if os.path.isdir(dst_path):
            shutil.rmtree(dst_path)
        os.mkdir(dst_path)

        for simulation_name, simulation_obj in self.simulations.items():
            simulation_path = os.path.join(dst_path, simulation_name)
            simulation_obj.archive(dst_path=simulation_path)

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
        X = self.cell_concentration_matrix

        # transform total molar fraction into a column vector
        c = [self.concentration[s] for s in self.symbols]

        # Since X c = f
        # f = X^{-1} . c
        # solve for phase molar fraction
        f = np.dot(linalg.inv(X), c)
        #for k in f:
        #    if k < 0:
                # msg = "phase molar fraction cannot be negative\n"
                # msg += "X:{}\n".format(X)
                # msg += "c:{}\n".format(c)
                # msg += "f:{}\n".format(f)
                # raise MultiCellError(msg)
        #        from scipy.optimize import nnls
        #        f, residuals = nnls(X,c)
        #        break 
        
        return {v:f[i] for i,v in enumerate(self.cell_names)}
        
    def get_number_of_atoms(self, symbol=None):
        n_atoms = 0
        for cell in self.cells.values():
            n_atoms += cell.get_number_of_atoms(symbol)
        return n_atoms

import math
import numpy as np
from typing import Dict, List
from mexm.elements import ELEMENTS
from mexm.structure import SimulationCell
from mexm.structure import make_super_cell

def get_bcc_cell(r: float) -> SimulationCell:
    #TODO:

    a0 = 4 / math.sqrt(3) * r
    simulation_cell = SimulationCell()
    simulation_cell.a0 = a0
    simulation_cell.add_atom('Ni',[0.0, 0.0, 0.0])
    simulation_cell.add_atom('Ni',[0.5, 0.5, 0.5])
    return simulation_cell

def get_fcc_cell(r: float) -> SimulationCell:
    #TODO:
    
    a0 = math.sqrt(8) * r
    simulation_cell = SimulationCell()
    simulation_cell.a0 = a0
    simulation_cell.add_atom('Ni',[0.0, 0.0, 0.0])
    simulation_cell.add_atom('Ni',[0.5, 0.5, 0.0])
    simulation_cell.add_atom('Ni',[0.5, 0.0, 0.5])
    simulation_cell.add_atom('Ni',[0.0, 0.5, 0.5])
    return simulation_cell

def get_hcp_cell():
    #TODO:
    simulation_cell = SimulationCell()
    return simulation_cell

def create_random_cell(
    cell_type: str,
    composition: Dict[str, float],
    supercell: List[int]
) -> SimulationCell:
    """
    Arguments:
        cell_type (str): the type of base cell structure to use.
        composition (List[str, float]): the composition of the alloy
        supercell (List[int]): the supercell 
    """

    assert isinstance(cell_type, str)
    assert isinstance(composition, dict)
    assert isinstance(supercell, list)
    assert all([isinstance(k, int) for k in supercell])

    sum_composition = sum(composition.values())
    composition_ = {
        k:v/sum_composition for k, v in composition.items()
    }
    max_atmrad = max([ELEMENTS[k].vdwrad for k in composition_])
    cell_type_to_unit_cell_map = {
        'bcc': get_bcc_cell,
        'fcc': get_fcc_cell,
        'hcp': get_hcp_cell
    } 
    cell = cell_type_to_unit_cell_map[cell_type](r=max_atmrad)
    cell = make_super_cell(structure=cell, sc=supercell)

    idx_atoms_all = [k for k in range(cell.n_atoms)]
    idx_atoms = {}

    atoms_processed = []
    for symbol in composition_:
        idx_atoms[symbol] = []
        
        sum_composition = sum([
            v for k, v in composition_.items() if k not in atoms_processed
        ])
        probability = composition_[symbol]/sum_composition

        n_atoms = int(cell.n_atoms * probability) \
            - sum([len(v) for v in idx_atoms.values()])

        for i_atom in range(n_atoms):
            idx = np.random.choice(idx_atoms_all, 1).tolist()[0]
            idx_atoms[symbol].append(idx)
            idx_atoms_all.remove(idx)
            atoms_processed.append(symbol)
    
    for symbol in idx_atoms:
        for idx_atom in idx_atoms[symbol]:
            cell.atomic_basis[idx_atom].symbol = symbol

    return cell
