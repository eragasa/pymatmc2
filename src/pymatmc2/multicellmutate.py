from abc import ABC, abstractmethod
import warnings
from typing import Tuple
import numpy as np

# imports from Materials Ex Machina
from mexm.structure import SimulationCell
from mexm.io.vasp import Poscar
from mexm.simulation import VaspSimulation

# import from pymatmc2
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell
from pymatmc2 import constants

class MultiCellMutateAlgorithm(ABC):
    def __init__(self):
        self.multicell = None

    @abstractmethod
    def mutate_multicell(
        self,
        multicell: MultiCell
    ) -> MultiCell:
        raise NotImplementedError

    @abstractmethod
    def accept_or_reject(
        self,
        multicell_initial: MultiCell, 
        multicell_candidate: MultiCell,
        temperature, float
    ) -> Tuple[bool, MultiCell]:
        raise NotImplementedError

class IntraphaseSwap(MultiCellMutateAlgorithm):
    mutate_type = 'intraphase_swap'
    """ implements the Interphase swap algorithm 
    
    The interphase swap algorithm swaps atomic positions within each 
    phase
    """

    def mutate_multicell(
        self,
        multicell: MultiCell
    ) -> MultiCell:
        """

        This algorithm mutates each phase by an swapping the positions
        of two atoms of different species.

        Arguments:
            multicell (MultiCell): the cell to be mutated

        Returns:
            MultiCell: the mutated cell
        """

        rtn_multicell = MultiCell.initialize_from_obj(multicell=multicell)
        
        for phase in multicell.simulations:
            simulation = multicell.simulations[phase]

            if isinstance(simulation, VaspSimulation):
                cell = SimulationCell.initialize_from_object(
                    obj = simulation.contcar
                )
            atom_types = [v.symbol for v in cell.atomic_basis]
            
            symbols = list(cell.symbols)
            symbol1 = np.random.choice(symbols, 1)[0]
            symbols.remove(symbol1)
            symbol2 = np.random.choice(symbols, 1)[0]
            
            idx_symbol_1 = [
                i for i, v in enumerate(cell.atomic_basis) if v.symbol == symbol1
            ]
            idx_symbol_2 = [
                i for i, v in enumerate(cell.atomic_basis) if v.symbol == symbol2
            ]

            # the index of the atoms i am swapping chosen randomly
            idx_atm_1 = np.random.choice(idx_symbol_1, 1)[0]
            idx_atm_2 = np.random.choice(idx_symbol_2, 1)[0]

            cell.atomic_basis[idx_atm_1].symbol = symbol2
            cell.atomic_basis[idx_atm_2].symbol = symbol1
            
            if isinstance(simulation, VaspSimulation):
                rtn_multicell.simulations[phase].poscar \
                    = Poscar.initialize_from_object(obj=cell)
        return rtn_multicell

    def acceptance_probability(
        self, 
        E0: float,
        E1: float,
        temperature: float
    ) -> float:
        """
        
        Arguments:
            E0 (float): energy of the initial phase in eV
            E1 (float): energy of the candidate phase in eV
            temperature (float): temperature of the simulation in K

        Note:
        
        There is a possible runtime error which is suppressed due to
        ..math:
          exp ( frac{1}{k_B T} (E_1 - E_0) )


        """

        kB = constants.BOLTZMANN
        
        # metropolis probability of acceptance
        with warnings.catch_warnings():
            # turned off warnings here due to possible overflow 
            
            warnings.filterwarnings('ignore')

            p_accept = min(
                1,
                np.exp(-(E1 - E0)/(kB*temperature))
            )

        return p_accept

    def accept_or_reject(
        self,
        multicell_initial: MultiCell, 
        multicell_candidate: MultiCell,
        temperature: float
    ) -> Tuple[bool, MultiCell]:

        # object to return, since each phase has an acceptance/rejection 
        # criteria      
        rtn_multicell = MultiCell()
        rtn_multicell.simulations = {}

        phases = [k for k in multicell_initial.simulations]
        for phase in phases:
            E0 = multicell_initial.simulations[phase].total_energy
            E1 = multicell_candidate.simulations[phase].total_energy
            if E1 < E0:
                is_accept_phase = True
            else:
                p_accept = self.acceptance_probability(
                    E0=E0, 
                    E1=E1, 
                    temperature=temperature
                )
                if np.random.random() < p_accept:
                    is_accept_phase = True
                else:
                    is_accept_phase = False

            if is_accept_phase:
                rtn_multicell.simulations[phase] \
                    = multicell_candidate.simulations[phase]
            else:
                rtn_multicell.simulations[phase] \
                    = multicell_initial.simulations[phase]
            
            return rtn_multicell

class InterphaseSwap(MultiCellMutateAlgorithm):
    mutate_type = 'interphase_swap'

    def mutate_multicell(
        self,
        multicell: MultiCell
    ) -> MultiCell:
        raise NotImplementedError


    def accept_or_reject(
        self,
        multicell_initial: MultiCell, 
        multicell_candidate: MultiCell,
    ) -> Tuple[bool, MultiCell]:
        raise NotImplementedError


class IntraphaseFlip(MultiCellMutateAlgorithm):
    mutate_type = 'intraphase_flip'

    def mutate_multicell(
        self,
        multicell: MultiCell
    ) -> MultiCell:

        rtn_multicell = MultiCell.initialize_from_obj(multicell=multicell)
        for phase in multicell.simulations:
            simulation = multicell.simulations[phase]

            if isinstance(simulation, VaspSimulation):
                cell = SimulationCell.initialize_from_object(
                obj = simulation.contcar
            )

            idx_atoms = list(range(cell.n_atoms))
            idx_atom = np.random.choice(idx_atoms, 1)[0]

            old_symbol = cell.atomic_basis[idx_atom].symbol
            symbols = cell.symbols
            symbols.remove(old_symbol)
            new_symbol = np.random.choice(symbols, 1)[0]

            cell.atomic_basis[idx_atom].symbol = new_symbol
        
            if isinstance(simulation, VaspSimulation):
                 rtn_multicell.simulations[phase].poscar \
                    = Poscar.initialize_from_object(obj=cell)
               
        return rtn_multicell

    def acceptance_probability(
        self, 
        E0: float, 
        E1: float, 
        temperature: float
    ):
        kB = constants.BOLTZMANN
        T = temperature
        
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            p_accept = min(
                1,
                np.exp(- (1/(kB*T)) * (E1 - E0))
            )

        return p_accept

    def accept_or_reject(
        self,
        multicell_initial: MultiCell, 
        multicell_candidate: MultiCell,
        temperature: float
    ) -> Tuple[bool, MultiCell]:
        

        E0 = multicell_initial.total_energy
        E1 = multicell_candidate.total_energy

        if E1 < E0:
            is_accept = True
        else:
            p_accept = self.acceptance_probability(
                E0 = E0, 
                E1 = E1, 
                temperature = temperature
            )
            
            if np.random.random() < p_accept:
                is_accept = True
            else:
                is_accept = False

        if is_accept:
            return multicell_candidate
        else:
            return multicell_initial

class MultiCellMutateAlgorithmFactory(ABC):
    factories = {
        'interphase_swap': InterphaseSwap,
        'intraphase_swap': IntraphaseSwap,
        'intraphase_flip': IntraphaseFlip
    }

    def __init__(self):
        self.configuration = None
        self.mutation_types = None
        self.mutation_weights = None
        self.cumulative_weights = None

    def configure(self, configuration: Pymatmc2Configuration):
        self.configuration = configuration
        self.mutation_types = []
        self.mutation_weights = []
        for k, v in configuration.mutation_weights.items():
            self.mutation_types.append(k)
            self.mutation_weights.append(v)
        
        self.cumulative_weights = []
        for k in range(len(self.mutation_weights)):
            self.cumulative_weights.append(
                sum(self.mutation_weights[:k+1])
            )

        self.cell_names = self.configuration.cell_names

    def read_simulations(self, i_iteration):
        self.multicell_initial = MultiCell()
        self.multicell_initial.configuration = self.configuration
        self.multicell_candidate = MultiCell()
        self.multicell_candidate.conditation = self.configuration
        

    def accept_or_reject(
        self, 
        multicell_initial: MultiCell, 
        multicell_candidate: MultiCell,
        temperature: float,
        mutate_type: str
    ) -> Tuple[bool, MultiCell]:
        """

        Returns:
            bool: True if the candidate structure is accepted
            MultiCell: returns the new structure
        """

        mutator = self.factories[mutate_type]()
        mutator.configuration = self.configuration
        is_accept, multicell = mutator.accept_or_reject(
            multicell_initial = multicell_initial,
            multicell_candidate = multicell_candidate,
            temperature = temperature
        )

        return is_accept, multicell

    def mutate_cells(self, multicell: MultiCell) -> MultiCell:
        self.mutate_type = self.determine_mutate_algorithm()
        mutator = self.factories[self.mutate_type]()
        return mutator.mutate_multicell(multicell = multicell)

    def determine_mutate_algorithm(self):
        probability = np.random.random()

        for i, p in enumerate(self.cumulative_weights):
            if probability < p:
                mutate_type = self.mutation_types[i]
                break
        self.mutate_type = mutate_type
        return mutate_type

    
