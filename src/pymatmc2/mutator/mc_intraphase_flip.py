from copy import deepcopy

# typing libraries
from typing import Tuple
from typing import Optional

# linear algebra imports
import numpy as np
from numpy import linalg
import warnings

# Materials Ex Machina imports
from mexm.structure import SimulationCell
from mexm.io.vasp import Poscar
from mexm.io.vasp import Contcar
from mexm.simulation import VaspSimulation

# Pymatmc2 import
from pymatmc2 import constants
from pymatmc2 import MultiCell
from pymatmc2.mutator import BaseMultiCellMutator

class IntraphaseFlipMutator(BaseMultiCellMutator):
    mutate_type = 'intraphase_flip'

    @property
    def n_atoms_to_flip_max(self):
        return self.n_atoms_to_flip_max_

    @n_atoms_to_flip_max.setter
    def n_atoms_to_flip_max(self, n:int):
        if not isinstance(n, int):
            raise TypeError('n_atoms_to_flip_max must be an integer')
        self.n_atoms_to_flip_max_ = n

    @property
    def n_atoms_to_flip_min(self):
        return self.n_atoms_to_flip_min_

    @n_atoms_to_flip_min.setter
    def n_atoms_to_flip_min(self, n:int):
        if not isinstance(n, int):
            raise TypeError('n_atoms_to_flip_min must be an integer')
        self.n_atoms_to_flip_min_ = n


    
    def mutate_cell(
        self, 
        cell: SimulationCell,
        is_new_symbol_different= False
    ) -> SimulationCell:
        """ mutate a single cell

        Arguments:
            cell (SimulationCell)
            is_new_symbol_different (bool)
        Returns:
            SimulationCell: the mutated simulation cell
        """

        if not issubclass(cell, SimulationCell):
            msg = "cell arguments must be mexm.str ture.SimulationCell type"
            raise TypeError(msg)

        rtn_cell = deepcopy(cell)

        n_atoms = cell.n_atoms
        idx_atoms = list(range(cell.n_atoms))
        idx_atom = np.random.choice(idx_atoms, 1)[0]

        if is_new_symbol_different:
            # force the symbol to be different
            old_symbol = cell.atomic_basis[idx_atom].symbol
            symbols = list(self.configuration.symbols)
            symbols.remove(old_symbol)
            new_symbol = np.random.choice(symbols, 1)[0]
        else:
            # the symbol can be the same
            old_symbol = cell.atomic_basis[idx_atom].symbol
            symbols = list(self.configuration.symbols)
            new_symbol = np.random.choice(symbols, 1)[0]

        # assign new symbol
        rtn_cell.atomic_basis[idx_atom].symbol = new_symbol
    
        return rtn_cell

    def mutate_multicell(self, 
        multicell: MultiCell, 
        n_atoms_to_flip_min:Optional[int] = None,
        n_atoms_to_flip_max:Optional[int] = None, 
        is_debug:Optional[bool] = False
    ) -> MultiCell:
        """ create a new candidate multicell

        Create a deepcopy of multicell to the attribute multicell_initial.
        The number of atoms chosen is determined by an equal probability
        mass function over the open interval of integers from 
        n_atoms_to_flip_min to n_atoms_to_flip_max

        Arguments:
            multicell (MultiCell): the initial multicell 
            is_debug (Optional[bool] = False): if set to True, create debug output to stdout

        Returns:
            (MultiCell): the mutated multicell.  This results is also stored
                as the attribute 'mc_candidate'
        """
        self.multicell_initial = deepcopy(multicell)
        self.multicell_candidate = deepcopy(multicell)
        
        if n_atoms_to_flip_min is None:
            n_atoms_to_flip_min = self.n_atoms_to_flip_min
        else:
            self.n_atoms_to_flip_min = n_atoms_to_flip_min

        if n_atoms_to_flip_max is None:
            n_atoms_to_flip_max = self.n_atoms_to_flip_max
        else:
            self.n_atoms_to_flip_max = n_atoms_to_flip_max
                    
        # determine phases to alter
        
        for phase in multicell.simulations:
            simulation = multicell.simulations[phase]

            # extract the structure file
            if isinstance(simulation, VaspSimulation):
                if isinstance(simulation.contcar, Poscar):
                    cell = deepcopy(simulation.contcar)
                else:
                    cell = deepcopy(simulation.poscar)
            else:
                raise ValueError('unknown simulation type')
            assert isinstance(cell, SimulationCell)

            # mutate the cell
            # determine number of atoms to flip
            n_atoms_to_flip = 0
            if n_atoms_to_flip_min == n_atoms_to_flip_max:
                n_atoms_to_flip == n_atoms_to_flip_min
            elif n_atoms_to_flip_max > n_atoms_to_flip_min:
                n_atoms_to_flip_options = np.linspace(
                    start=n_atoms_to_flip_min, 
                    stop=n_atoms_to_flip_max,
                    num=n_atoms_to_flip_max - n_atoms_to_flip_min
                )
                n_atoms_to_flip = np.random.choice(n_atoms_to_flip_options, 1)
            else:
                msg = 'n_atoms_to_flip ({}) must be less than n_atoms_flip_max ({})'
                raise ValueError(
                    msg.format(
                        n_atoms_to_flip_min,
                        n_atoms_to_flip_max
                    )
                )

            # select atom at random
            idx_atoms = list(range(cell.n_atoms))
            idx_atoms_to_flip = np.random.choice(idx_atoms, n_atoms_to_flip)[0]
            old_symbols = cell.atomic_basis[idx_atoms_to_flip].symbol
                    
            # flip the symbol
            for old_symbol in old_symbols:
                symbols = cell.symbols
                symbols.remove(old_symbol)
                new_symbol = np.random.choice(symbols, 1)[0]
                    
                # assign new symbol
                cell.atomic_basis[idx_atom].symbol = new_symbol
                if isinstance(simulation, VaspSimulation):
                    self.multicell_candidate.simulations[phase].poscar = Poscar.initialize_from_object(obj=cell)
                else:
                    raise ValueError("unknown simulation type")
       
        assert isinstance(self.multicell_candidate, MultiCell) 
        return self.multicell_candidate

    def acceptance_probability(
        self, 
        E_initial: float, 
        E_candidate: float, 
        temperature: float,
        pressure: float
    ) -> float:
        """ Calculate the acceptance/rejection probability

        This calculates the probability of accepting the the changes of a 
        multicell flip. 
        """

        kB = constants.BOLTZMANN
        T = temperature
        beta = 1/kB/T
        n_phases = len(self.configuration.cell_names)
 
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')

            E1 = E_candidate
            E0 = E_initial

            p_accept = min(1,np.exp(-n_phases*(E1-E0)*beta))

        self.p_accept = {}
        self.p_accept['total'] = p_accept
        return p_accept

    def accept_or_reject(
        self,
        multicell_initial: MultiCell, 
        multicell_candidate: MultiCell,
        temperature: float,
        pressure: float,
        is_debug: Optional[bool] = False
    ) -> Tuple[bool, MultiCell]:
        
        if not isinstance(temperature, float):
           msg = 'temperature must be a numeric type'
           raise TypeError(msg)

        if not isinstance(pressure, float):
            msg = 'pressure must be a numeric type'
            raise TypeError(msg)

        self.multicell_initial = deepcopy(multicell_initial)
        self.multicell_candidate = deepcopy(multicell_candidate)

 
        E0 = self.multicell_initial.total_energy
        E1 = self.multicell_candidate.total_energy

        
        self.is_accept = {}
        for k in self.configuration.cell_names:
            self.is_accept[k] = None
        self.is_accept['total'] = None

        self.p_accept = {}
        for k in self.configuration.cell_names:
            self.p_accept[k] = None
        self.p_accept['total'] = None

        self.p_accept_rnd = {}
        for k in self.configuration.cell_names:
            self.p_accept_rnd[k] = None
        self.p_accept_rnd['total'] = None        

        if E1 < E0:
            self.is_accept['total'] = True
        else:
            self.p_accept['total'] = self.acceptance_probability(
                E_initial = E0, 
                E_candidate = E1, 
                temperature = temperature,
                pressure = pressure
            )
 
            self.p_accept_rnd['total'] = np.random.random()
            if self.p_accept_rnd['total'] < self.p_accept['total']:
                self.is_accept['total'] = True
            else:
                self.is_accept['total'] = False

        self.multicell_final = deepcopy(self.multicell_candidate)

        return self.is_accept['total'], self.multicell_final
