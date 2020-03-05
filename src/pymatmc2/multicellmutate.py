from pymatmc2 import MultiCell

from abc import ABC

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