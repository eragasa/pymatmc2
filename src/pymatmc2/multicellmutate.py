from copy import deepcopy
from abc import ABC, abstractmethod
import warnings
from typing import Tuple, List, Dict
import numpy as np
import numpy.linalg as linalg

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
        self._configuration = None
        self.multicell = None
        self._multicell_initial = None
        self._multicell_candidate = None
        self._multicell_final = None

    @property
    def configuration(self) -> Pymatmc2Configuration:
        return self._configuration

    @configuration.setter
    def configuration(self, configuration: Pymatmc2Configuration):
        self._configuration = configuration

    @property
    def multicell_initial(self) -> MultiCell:
        return self._multicell_initial

    @multicell_initial.setter
    def multicell_initial(self, mc: MultiCell):
        assert isinstance(mc, MultiCell)
        self._multicell_initial = mc

    @property
    def multicell_candidate(self) -> MultiCell:
        return self._multicell_candidate

    @multicell_candidate.setter
    def multicell_candidate(self, mc: MultiCell):
        assert isinstance(mc, MultiCell)
        self._multicell_candidate = mc

    @property
    def multicell_final(self) -> MultiCell:
        return self._multicell_final

    @multicell_final.setter
    def multicell_final(self, mc: MultiCell):
        assert isinstance(mc, MultiCell)
        self._multicell_final = mc
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
        temperature: float
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
        is_accept_phase_array = []
        for phase in phases:
            n_atoms_0 = multicell_initial.simulations[phase].poscar.n_atoms
            E0 = multicell_initial.simulations[phase].total_energy/n_atoms_0

            n_atoms_1 = multicell_candidate.simulations[phase].poscar.n_atoms
            E1 = multicell_candidate.simulations[phase].total_energy/n_atoms_1
            assert isinstance(E0, float)
            assert isinstance(E1, float)
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
            is_accept_phase_array.append(is_accept_phase)
            self.multicell_candidate = rtn_multicell 
        return any(is_accept_phase_array), rtn_multicell

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

    def mutate_cell(self, phase_name: str) -> SimulationCell:
        simulation = self.multicell_initial.simulations[phase_name]

        if isinstance(simulation, VaspSimulation):
            if isinstance(simulation.contcar, Poscar):
                cell = deepcopy(simulation.contcar)
            else:
                cell = deepcopy(simulation.poscar)
        else:
            msg = "unknown simulation type: {}"
            msg = msg.format(type(simulation))
            raise TypeError(msg)

        idx_atoms = list(range(cell.n_atoms))
        idx_atom = np.random.choice(idx_atoms, 1)[0]

        #old_symbol = cell.atomic_basis[idx_atom].symbol 
        #symbols = cell.symbols
        #symbols.remove(old_symbol)
        symbols = cell.symbols
        new_symbol = np.random.choice(symbols, 1)[0]

        cell.atomic_basis[idx_atom].symbol = new_symbol

        if isinstance(simulation, VaspSimulation):
            self.multicell_candidate.simulations[phase_name].poscar \
                = Poscar.initialize_from_object(obj=cell)
        else:
            msg = "unknown simulation type: {}"
            msg = msg.format(type(simulation))
            raise TypeError(msg)
        
    def mutate_multicell(self, multicell: MultiCell) -> MultiCell:
        """ create a new candidate multicell

        Create a deepcopy of multicell to the attribute multicell_initial.


        Arguments:
            multicell (MultiCell): the initial multicell from which to get a    candidate
        """
        self.multicell_initial = deepcopy(multicell)
        self.multicell_candidate = deepcopy(multicell)
        
        is_good_multicell = False
        while not is_good_multicell:
            for phase in multicell.simulations:
                self.mutate_cell(phase_name=phase)
            try:
                self.multicell_candidate.phase_molar_fraction
                is_full_rank = True
            except linalg.LinAlgError:
                is_full_rank = False
            
            if is_full_rank:
                is_phase_fraction_good_array = []
                for f in self.multicell_candidate.phase_molar_fraction.values():
                    if f > 1 or f < 0:
                        is_phase_fraction_good_array.append(False)
                    else:
                        is_phase_fraction_good_array.append(True)
                is_valid_phase_fraction = all(is_phase_fraction_good_array)
            else:
                is_valid_phase_fraction = False
            
            is_different_concentration_array = []
            for cn in self.multicell_initial.cell_names:
                cn0 = self.multicell_initial.cell_concentration[cn]
                cn1 = self.multicell_candidate.cell_concentration[cn]
                if cn0 == cn1:
                    is_different_concentration_array.append(False)
                else:
                    is_different_concentration_array.append(True)
            is_different_concentration = any(is_different_concentration_array)

            if all([
                is_full_rank, 
                is_valid_phase_fraction,
                is_different_concentration
            ]):
                is_good_multicell = True
                
        return self.multicell_candidate

    def acceptance_probability(
        self, 
        E0: float, 
        E1: float, 
        temperature: float
    ):
        kB = constants.BOLTZMANN
        T = temperature
        beta = 1/kB/T
        n_phases = len(self.multicell_initial.cell_names) 
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')

            p_accept = min(1,np.exp(-n_phases*(E1-E0)*beta))

        self.p_accept = {}
        self.p_accept['total'] = p_accept
        return p_accept

    def accept_or_reject(
        self,
        multicell_initial: MultiCell, 
        multicell_candidate: MultiCell,
        temperature: float
    ) -> Tuple[bool, MultiCell]:
       

        self.multicell_initial = multicell_initial
        self.multicell_candidate = multicell_candidate

        E0 = self.multicell_initial.total_energy
        E1 = self.multicell_candidate.total_energy
        
        if E1 < E0:
            self.is_accept = True
        else:
            p_accept = self.acceptance_probability(
                E0 = E0, 
                E1 = E1, 
                temperature = temperature
            )
            
            self.p = np.random.random()
            if self.p < p_accept:
                self.is_accept = True
            else:
                self.is_accept = False

        if self.is_accept:
            self.multicell_final = self.multicell_candidate
        else:
            self.multicell_final = self.multicell_initial

        return self.is_accept, self.multicell_final

class MultiCellMutateAlgorithmFactory(ABC):
    factories = {
        'interphase_swap': InterphaseSwap,
        'intraphase_swap': IntraphaseSwap,
        'intraphase_flip': IntraphaseFlip
    }

    """
    Attributes:
        configuration (Pymatmc2Configuration)
    """

    @property
    def mutation_types(self) -> List[str]:
        assert isinstance(self.configuration, Pymatmc2Configuration)
        mutation_types = []
        for k in self.configuration.mutation_weights.keys():
            mutation_types.append(k)
        return mutation_types

    @property
    def mutation_weights(self) -> List[float]:
        assert isinstance(self.configuration, Pymatmc2Configuration)
        mutation_weights = []
        for v in self.configuration.mutation_weights.values():
            mutation_weights.append(v)
        return mutation_weights

    @property
    def cumulative_weights(self) -> List[float]:
        assert isinstance(self.configuration, Pymatmc2Configuration)
        mutation_weights = self.mutation_weights
        cumulative_weights = []
        for k in range(len(self.mutation_weights)):
            cum_w = sum(mutation_weights[:k+1]) 
            cumulative_weights.append(cum_w)
        return cumulative_weights

    @property
    def cell_names(self) -> List[str]:
        return self.configuration.cell_names

    def configure(self, configuration: Pymatmc2Configuration):
        self.configuration = configuration

    def read_simulations(self, i_iteration:int):
        self.multicell_initial = MultiCell()
        self.multicell_initial.configuration = self.configuration
        self.multicell_candidate = MultiCell()
        self.multicell_candidate.configuration = self.configuration

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

        assert isinstance(multicell_initial, MultiCell)
        assert isinstance(multicell_candidate, MultiCell)
        assert isinstance(temperature, float)
        assert isinstance(mutate_type, str)

        mutator = self.factories[mutate_type]()
        mutator.configuration = self.configuration
        is_accept, multicell = mutator.accept_or_reject(
            multicell_initial = multicell_initial,
            multicell_candidate = multicell_candidate,
            temperature = temperature
        )

        return is_accept, multicell

    def mutate_cells(self, multicell: MultiCell) -> Tuple[str, MultiCell]:
        
        self.mutate_type = self.determine_mutate_algorithm()
        
        mutator = self.factories[self.mutate_type]()

        self.multicell_initial = multicell
        self.multicell_candidate = mutator.mutate_multicell(multicell = multicell)


        return self.mutate_type, mutator.mutate_multicell(multicell = multicell)

    def determine_mutate_algorithm(self):
        probability = np.random.random()

        for i, p in enumerate(self.cumulative_weights):
            if probability < p:
                mutate_type = self.mutation_types[i]
                break
        self.mutate_type = mutate_type
        return mutate_type

    
