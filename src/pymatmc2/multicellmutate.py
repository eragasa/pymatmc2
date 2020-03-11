from abc import ABC
import numpy as np
from pymatmc2 import MultiCell


class MultiCellMutateAlgorithm(ABC):
    def __init__(self, multicell: MultiCell):
        self.multicell = multicell

class InterphaseSwap(MultiCellMutateAlgorithm):
    mutate_type = 'interphase_swap'
    pass

class IntraphaseSwap(MultiCellMutateAlgorithm):
    mutate_type = 'intraphase_swap'

class IntraphaseFlip(MultiCellMutateAlgorithm):
    mutate_type = 'intraphase_flip'

class MultiCellMutateAlgorithmFactory(ABC):
    factories = {
        'interphase_swap': InterphaseSwap,
        'intraphase_swap': IntraphaseSwap,
        'interphase_flip': IntraphaseFlip
    }

    def __init__(self):
        self.mutation_types = None
        self.mutation_weights = None
        self.cumulative_weights = None

    def configure(self, configuration):
        self.mutation_types = []
        self.mutation_weights = []
        for k, v in configuration['mutate_weights'].items():
            self.mutation_types.append(k)
            self.mutation_weights.append(v)
        
        self.cumulative_weights = []
        for k in range(len(self.mutation_weights)):
            self.cumulative_weights.append(
                sum(self.mutation_weights[:k+1])
            )

    def determine_mutate_algorithm(self):
        probability = np.random.random()

        for i, p in enumerate(self.cumulative_weights):
            if probability < p:
                mutate_type = self.mutation_types[i]
        return mutate_type

        
