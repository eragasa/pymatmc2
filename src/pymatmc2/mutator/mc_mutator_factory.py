from typing import List
from typing import Tuple
import numpy as np
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell
from pymatmc2.mutator import IntraphaseFlipMutator
from pymatmc2.mutator import IntraphaseSwapMutator
from pymatmc2.mutator import BaseMultiCellMutator

class MultiCellMutatorFactory(BaseMultiCellMutator):

    factories = {
        v.mutate_type: v for v in [IntraphaseFlipMutator, IntraphaseSwapMutator]
    }

    def __init__(self):
        super().__init__()
        self._mutate_type = None

    @property
    def mutate_type(self) -> str:
        return self._mutate_type

    @mutate_type.setter
    def mutate_type(self, mutate_type: str):
        self._mutate_type = mutate_type

    @property
    def multicell_mutators(self):
        factory_dict = {
            'intraphase_swap': IntraphaseSwapMutator,
            'intraphase_flip': IntraphaseFlipMutator
        }
        return factory_dict

    @property
    def mutation_types(self) -> List[float]:
        if not isinstance(self.configuration, Pymatmc2Configuration):
            raise TypeError('the configuration attribute has not been set')
        mutation_types = [k for k in self.configuration.mutation_weights.keys()]
        return mutation_types

    @property
    def mutation_weights(self) -> List[float]:
        if not isinstance(self.configuration, Pymatmc2Configuration):
            raise TypeError('the configuration attribute has not been set')
        mutation_weights = [v for v in self.configuration.mutation_weights.values()]
        return mutation_weights

    @property
    def cumulative_weights(self) -> List[float]:
        if not isinstance(self.configuration, Pymatmc2Configuration):
            raise TypeError('the configuration attribute has not been set')
        mutation_weights = self.mutation_weights
        cumulative_weights = []
        
        n_mutation_types = len(mutation_weights)
        for k in range(n_mutation_types):
            cum_w = sum(mutation_weights[:k+1]) 
            cumulative_weights.append(cum_w)
        return cumulative_weights

    def determine_mutate_algorithm(self) -> str:
        """
        Returns:
            bool: the mutation type
        """
        probability = np.random.random()

        for i, p in enumerate(self.cumulative_weights):
            if probability < p:
                mutate_type = self.mutation_types[i]
                break
        self.mutate_type = mutate_type
        return mutate_type


    def accept_or_reject(self, 
        multicell_initial: MultiCell, 
        multicell_candidate: MultiCell, 
        temperature: float,
        pressure: float,
        mutate_type = None
    ) -> Tuple[bool, MultiCell]:
        if mutate_type is None:
            if self.mutate_type is None:
                raise TypeError('cannot determine what kind of mutate_type to use')
        else:
            if mutate_type not in self.multicell_mutators:
                raise ValueError('{} is not a valid mutate_type'.format(mutate_type))
            self.mutate_type = mutate_type

        mutator = self.multicell_mutators[self.mutate_type]()
        mutator.configuration = self.configuration
        is_accept, multicell = mutator.accept_or_reject(
            multicell_initial = multicell_initial,
            multicell_candidate = multicell_candidate,
            temperature = temperature,
            pressure = pressure
        )
        return self.mutate_type, multicell

    def mutate_multicell(self, multicell: MultiCell) -> MultiCell:
        if self.mutate_type is None:
            self.mutate_type = self.determine_mutate_algorithm()
        mutator = self.multicell_mutators[self.mutate_type]()
        self.multicell_initial = multicell
        self.multicell_candidate = mutator.mutate_multicell(multicell = multicell)

        return self.mutate_type, self.multicell_candidate
