from copy import deepcopy

# typing libraries
from typing import Tuple

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
from pymatmc2.mutator import BaseMultiCellMutator
from pymatmc2 import MultiCell

class IntraphaseFlipMutator(BaseMultiCellMutator):
    mutate_type = 'intraphase_flip'

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

    def mutate_multicell(self, multicell: MultiCell) -> MultiCell:
        """ create a new candidate multicell

        Create a deepcopy of multicell to the attribute multicell_initial.


        Arguments:
            multicell (MultiCell): the initial multicell 

        Returns:
            (MultiCell): the mutated multicell.  This results is also stored
                as the attribute 'mc_candidate'
        """
        self.multicell_initial = deepcopy(multicell)
        self.multicell_candidate = deepcopy(multicell)
        while True:
            try:
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
                    # select atom at random
                    idx_atoms = list(range(cell.n_atoms))
                    idx_atom = np.random.choice(idx_atoms, 1)[0]
                    old_symbol = cell.atomic_basis[idx_atom].symbol
                    
                    # flip the symbol
                    symbols = cell.symbols
                    symbols.remove(old_symbol)
                    new_symbol = np.random.choice(symbols, 1)[0]
                    
                    # assign new symbol
                    cell.atomic_basis[idx_atom].symbol = new_symbol
                
                    if isinstance(simulation, VaspSimulation):
                        rtn_multicell.simulations[phase].poscar \
                            = Poscar.initialize_from_object(obj=cell)
                    else:
                        raise ValueError("unknown simulation type")
                    
                # this will raise a LinAlg error if rank deficient
                self.multicell_candidate.phase_molar_fraction
                break
            except linalg.LinAlgError:
                pass
        
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
        pressure: float
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
