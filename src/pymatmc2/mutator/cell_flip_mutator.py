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

from typing import List
from copy import deepcopy
import numpy as np
from mexm.structure import SimulationCell
from pymatmc2 import Pymatmc2Configuration

class CellFlipMutator(object):
    
    def __init__(self):
        super().__init__(self)
    
    def get_cell_differences(self):
        """ return the differences between the initial and candidate cells

        Returns:
            List[int, str, str]: the first entry is the index of the atom 
                which was flipped,  the second entry the initial chemical
                symbol of the atom which was flipped, and the third entry
                is candidate chemical symbol of the atom which was flipped.
        """
        diffs = []
        n_atoms = self.cell_initial.n_atoms
        for i in range(n_atoms):
            symbol_initial = self.cell_initial.atomic_basis[i].symbol
            symbol_candidate = self.cell_candidate.atomic_basis[i].symbol
            if symbol_initial != symbol_candidate:
                diffs.append[i, symbol_initial, symbol_candidate]
        return diffs

    def mutate(self, force_different_symbol=False) -> SimulationCell:
        """
        Raises:
            ValueError: if configuration attribute has not been set
        """ 

        if self.configuration is None:
            msg = "the configuration attribute has not been set"
            raise ValueError(msg)

        if not issubclass(self.configuration, Pymatmc2Configuration):
            msg = "the configuration attribute must be a subclass of Pymatmc2Configuration"
            raise TypeError(msg)

        if self.cell_initial is None:
            msg = "cell_initial attribute has not been set"
            raise ValueError(msg)
        
        if not issubclass(self.cell_initial, SimulationCell):
            msg = "cell_initial must be a subclass of mexm.structure.SimulationCell"
            raise TypeError(msg)

        n_atoms = self.cell_initial.n_atoms
        idx_atoms = list(range(n_atoms))
        idx_atom = np.random.choice(idx_atoms, 1)[0]

        if force_different_symbol:
            old_symbol = self.cell_initial.atomic_basis[idx_atom].symbol
            symbols = self.configuration.symbols
            new_symbol = np.random.choice(symbols, 1)[0]
        else:
            symbols = self.cell_initial.symbols
            new_symbol = np.random.choice(symbols, 1)[0]
        
        self.cell_candidate = deepcopy(self.cell_initial)
        self.cell_candidate.atomic_basis[idx_atom].symbol = new_symbol

        return self.cell_candidate
    

if __name__ == "__main__":
    obj = CellFlip()
