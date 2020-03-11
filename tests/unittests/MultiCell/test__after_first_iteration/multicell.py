# coding: utf-8
# Copyright (c) Eugene J. Ragasa
# Distributed under the terms of the MIT License

""" multicell montecarlo class

This module implements holding simulation classes for multiple cells

TODO:
[3/5/2020] - Generalize this class to use the SimulationFactory.
"""

__author__ = "Eugene J. Ragasa"
__email__ = "ragasa.2@osu.edu"
__copyright__ = "Copyright 2020, Eugene J. Ragasa"
__maintainer__ = "Eugene J. Ragasa"
__date__ = "2020/02/22"

from collections import OrderedDict
from typing import List
import numpy as np
from numpy import linalg
from mexm.io.vasp import Poscar
from mexm.simulation import VaspSimulation
from pymatmc2 import Pymatmc2Configuration

class MultiCell:
    """

    Attributes:
        cells (Dict[str, Poscar])
        molar_fraction_total (List[(float)])
    """

    def __init__(self):
        """
        """
        self.molar_fraction_total = None
        self.n_cells = None
        self.cells = None
        self.simulations = None

    @staticmethod
    def initialize_from_pymatmc2_configuration(
        configuration: Pymatmc2Configuration
    ):
        """
            configuration (Pymatmc2Configuration)
        """
        obj = MultiCell()
        obj.cells = {}

        simulation_cells = configuration.simulation_cells
        order_simulation_cells = OrderedDict(
            sorted(simulation_cells.items())
        )
        for k, v in order_simulation_cells.items():
            obj.cells[k] = Poscar()
            obj.cells[k].read(path=v['poscar'])

        molar_fraction_total = configuration.molar_fraction_total
        obj.molar_fraction_total = {
            k:v/sum(molar_fraction_total.values()) for k, v in molar_fraction_total.items()
        }
        return obj
    
    def add_simulation(
        self,
        simulation_type: str, 
        simulation_name: str,
        simulation_path: str
    ):
        if self.simulations is None:
            self.simulations = {}

        if self.cells is None:
            self.cells = {}

        if simulation_type == 'vasp':
            self.simulations[simulation_name] = VaspSimulation()
            self.simulations[simulation_name].read(simulation_path)
            self.cells[simulation_name] = self.simulations[simulation_name].poscar

    def configure(self, configuration: Pymatmc2Configuration):
        """ configure class from a Pymatmc2Configuration
        
        this method sets up the cells attribute

        Arguments:
            configuration (Pymatmc2Configuration): instance of the 
                Pymatmc2Configuration in which to setup this class.
         """
        self.cells = {}
        for k,v in configuration.simulation_cells.items():
            self.cells[k] = Poscar()
            self.cells[k].read(path=v['poscar'])            

    def add_cells_from_dict(self, simulation_cells):
        """ add cells from dictionary object

        simulation_cells = {
            '<cell_name>' = {'<file_type'>:'<path_to_file>'}
        }
        Arguments:
            simulation_cells (dict) = dictionary of simulation cells
        """

    @property
    def cell_names(self) => List[str]:
        ordered_cells = OrderedDict(sorted(self.cells.items()))
        return [k for k in ordered_cells.keys()]
        
    @property
    def cell_molar_fraction(self) -> List[List[float]]:
        X = []
        ordered_cells = OrderedDict(sorted(self.cells.items()))
        for v in ordered_cells.values():
            n_atoms = [v.get_number_of_atoms(s) for s in self.symbols]
            X.append([k/sum(n_atoms) for k in n_atoms])
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
        c = [self.molar_fraction_total[s] for s in self.symbols]
        return c
    
    @property
    def phase_molar_fraction(self):
        X = np.array(self.cell_molar_fraction)
        c = np.array(self.total_molar_fraction)

        f = np.dot(linalg.inv(X), c)
        for k in f:
            pass
            # assert k > 0
        return f.tolist()
        
    def get_number_of_atoms(self, symbol=None):
        n_atoms = 0
        for cell in self.cells.values():
            n_atoms += cell.get_number_of_atoms(symbol)
        return n_atoms

    def modify_atoms(self, modify_type: str, modify_options):
        """

        Arguments:
            modify_type (str):  Supported options are 'interphase_swap', 'species_flip', and 'intraphase_swap'
        """
        modify_dict = {
            'interphase_swap': self.modify_atoms_interphase_swap,
            'species_flip':  self.modify_atoms_species_flip,
            'intraphase_swap': self.modify_atoms_species_flip
        }
