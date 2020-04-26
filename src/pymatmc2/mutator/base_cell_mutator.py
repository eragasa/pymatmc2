# coding: utf-8
# Copyright (c) Eugene J. Ragasa
# Distributed under the terms of the MIT License

""" CellClip class

the module implements the CellFlip class
"""

__all__ = ["CellFlip"]
__author__ = "Eugene J. Ragasa"
__email__ = "ragasa.2@osu.edu"
__copyright__ = "Copyright 2020, Eugene J. Ragasa"
__maintainer__ = "Eugene J. Ragasa"
__date__ = "2020/04/20"

from abc import ABC, abstractmethod
from typing import List
from typing import Tuple
from copy import deepcopy
import numpy as np
from mexm.structure import SimulationCell
from pymatmc2 import Pymatmc2Configuration

class BaseCellMutator(ABC):
    
    def __init__(self):
        self._configuration = None
        self._cell_initial = None
        self._cell_candidate = None
        self._cell_final = None
    
    @property
    def configuration(self) -> Pymatmc2Configuration:
        """:Pymatmc2Configuration """
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
    def cell_candidate(self) -> SimulationCell:
        return self._cell_candidate

    @cell_candidate.setter
    def cell_candidate(self, cell: SimulationCell):
        self._cell_candidate = cell
    
    @property
    def cell_final(self) -> SimulationCell:
        return self._cell_final
        
    @cell_final.setter
    def cell_final(self, cell: SimulationCell):
        self._cell_final = cell

    @abstractmethod
    def get_cell_differences(self):
        """ return the differences between the initial and candidate cells

        Returns:
            List[Tuple[int, str, str]]: the first entry is the index of the atom 
                which was flipped,  the second entry the initial chemical
                symbol of the atom which was flipped, and the third entry
                is candidate chemical symbol of the atom which was flipped.
        """
        raise NotImplementedError

    @abstractmethod
    def mutate(self, force_different_symbol=False) -> SimulationCell:
        """
        Raises:
            ValueError: if configuration attribute has not been set
        """
        raise NotImplementedError   

if __name__ == "__main__":
    obj = CellFlip()
