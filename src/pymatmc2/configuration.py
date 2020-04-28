# coding: utf-8
# Copyright (c) Eugene J. Ragasa
# Distributed under the terms of the MIT License

""" configuration class for pymatmc2

This module implmements an class for convenient access to configure
the pymatmc2 classes
"""

__author__ = "Eugene J. Ragasa"
__email__ = "ragasa.2@osu.edu"
__copyright__ = "Copyright 2020, Eugene J. Ragasa"
__maintainer__ = "Eugene J. Ragasa"
__date__ = "2020/02/22"

import os
import copy
import yaml
from collections import OrderedDict
from typing import Dict, List
from mexm.io.filesystem import OrderedDictYAMLLoader

class Pymatmc2Configuration():
    """ configuration file for pymatmc2

    Attributes:
        path (str): the path to the directory
        configuration (dict): dictionary of the configuration
    """
    def __init__(self):
        self._path = None
        self.configuration = None
        pass

    @staticmethod
    def initialize_from_dict(config_dict: dict):
        """
        Arguments:
            config_dict (dict): configuration dictionary
        """
        obj = Pymatmc2Configuration()
        obj.configure_from_dict(config_dict=config_dict)

    def configure_from_dict(self, config_dict: dict):
        self.configuration = copy.deepcopy(config_dict)
        
    @property
    def path(self):
        """(str): absolute path of the configuration file"""
        
        return self._path

    @path.setter
    def path(self, path):

        # convert to absolute
        if isinstance(path, str):
            if os.path.isabs(path):
                self._path = path
            else:
                self._path = os.path.abspath(path)
        else:
            msg = "path must be a string"
            raise TypeError(msg)

    @property
    def calculator_type(self) -> str:
        """:(str) the energy calculator code being code """
        return self.configuration['calculator']['calculator']

    @property
    def simulation_type(self) -> str:
        """:(str) the simulation type of the code being used """
        return self.configuration['calculator']['simulation_type']

    @property
    def n_cells(self) -> int:
        """:(int) number of simulation_cells """     
        return len(self.configuration['atomic_configuration']['simulation_cells'])

    @property
    def simulation_cells(self) -> Dict[str, Dict[str, str]]:
        """:dict simulations cells"""
        return self.configuration['atomic_configuration']['simulation_cells']
    
    @property
    def cell_names(self) -> List[str]:
        """:List(str) names of the cells"""
        cell_names = [k for k in self.simulation_cells] 
        return cell_names

    @property
    def phase_point_string(self) -> str:
        T = int(self.temperature)
        P = int(self.pressure)

        fmt = "{}K_{}GPa"
        return  fmt.format(T, P)

    def get_iteration_string(self, i: int) -> str:
        fmt = '{:05}'

        return fmt.format(i)

    @property
    def temperature(self) -> float:
        temperature = self.configuration['environment_variables']['temperature']
        return temperature

    @property
    def pressure(self) -> float:
        return self.configuration['environment_variables']['pressure']

    @property
    def concentration(self) -> Dict[str, float]:

        # make local copy of the composition
        concentration = OrderedDict()
        for k, v in self.configuration['atomic_configuration']['molar_fraction_total'].items():
            concentration[k] = v

        # sum of composition, should be equal to one.  force normalization.
        sum_concentration = sum(concentration.values())
        for k, v in concentration.items():
            concentration[k] = v/sum_concentration

        return concentration

    @property
    def symbols(self) -> List[str]:
        
        symbols = []
        for k in self.configuration['atomic_configuration']['molar_fraction_total']:
            symbols.append(k)

        return symbols

    @property
    def hpc_manager(self) -> dict:
        return self.configuration['hpc_manager']
        
    # deprecate
    @property
    def molar_fraction_total(self) -> Dict[str, float]:
        
        # make local copy of the composition
        composition = OrderedDict()
        for k, v in self.configuration['atomic_configuration']['molar_fraction_total'].items():
            composition[k] = v

        # sum of composition, should be equal to one.  force normalization.
        sum_composition = sum(composition.values())
        for k, v in composition.items():
            composition[k] = v/sum_composition

        return composition

    @property
    def mutation_weights(self):
        mutation_weights = OrderedDict()
        for k, v in self.configuration['mutation_weights'].items():
            mutation_weights[k] = v

        return mutation_weights

    @mutation_weights.setter
    def mutation_weights(self, mutation_weights: Dict[str, float]):
        self.configuration['mutation_weights'] = OrderedDict()
        for k, v in mutation_weights.items():
            self.configuration['mutation_weights'][k] = v
    
    @property
    def max_iterations(self):
        return self.configuration['max_iterations']

    @property
    def results_path(self):
        return self.configuration['results']['dir']

    @property
    def results_tar_path(self):
        return self.configuration['results']['tar_path']

    @property
    def results_file_path(self):
        return self.configuration['results']['file_path']

    def read(self, path):
        """ read the configuration files

        Arguments:
            path (str): path to the configuration file
        """
        self.path = path

        try:
            with open(self.path, 'r') as f:
                self.configuration= yaml.load(f, OrderedDictYAMLLoader)
        except FileNotFoundError:
            raise

    def write(self, path):
        """ write the configuration files

        Arguments:
            path (str): path to the configuration file
        """
        self.path = path

        with open(self.path, 'w') as f:
            yaml.dump(
                self.configuration, 
                f,
                default_flow_style=False
            )

