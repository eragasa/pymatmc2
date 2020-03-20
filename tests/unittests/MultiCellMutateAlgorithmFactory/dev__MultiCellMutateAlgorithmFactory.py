import os
import numpy as np
from pymatmc2 import Pymatmc2Configuration
from pymatmc2.multicellmutate import MultiCellMutateAlgorithmFactory

configuration_path = 'pymatmc2.config'
configuration = Pymatmc2Configuration()
configuration.read(configuration_path)

def dev__initialize():
    mcmutate = MultiCellMutateAlgorithmFactory()
    assert mcmutate.mutation_weights is None
    assert mcmutate.mutation_types is None
    assert mcmutate.cumulative_weights is None

def dev__configure():
    print(80*'-')
    print('MultiCellMutateAlgorithmFactory.configure')
    print(80*'-')

    mcmutate = MultiCellMutateAlgorithmFactory()
    mcmutate.configure(configuration=configuration)

    print('mutation_types:', mcmutate.mutation_types)
    print('mutation_weights:', mcmutate.mutation_weights)
    print('cumulative_weights:', mcmutate.cumulative_weights)

def dev__determine_mutate_algorithm():
    print(80*'-')
    print('MultiCellMutateAlgorithmFactory.determine_mutate_algorithm')
    print(80*'-')
    
    mcmutate = MultiCellMutateAlgorithmFactory()
    mcmutate.configure(configuration=configuration)
    mutate_type = mcmutate.determine_mutate_algorithm()
    
    print('mutate_type:', mutate_type)

if __name__ == "__main__":
    dev__initialize()
    dev__configure()
    dev__determine_mutate_algorithm()
