# coding: utf-8
# Copyright (c) Eugene J. Ragasa
# Distributed under the terms of the MIT License

""" CellClip class

the module implements the CellFlip class
"""

__author__ = "Eugene J. Ragasa"
__email__ = "ragasa.2@osu.edu"
__copyright__ = "Copyright 2020, Eugene J. Ragasa"
__maintainer__ = "Eugene J. Ragasa"
__date__ = "2020/04/20"

from pymatmc2.mutator import BaseMultiCellMutator

class IntraphaseSwapMutator(BaseMultiCellMutator):
    mutate_type = 'intraphase_swap'
    
    def accept_or_reject(self):
        pass

    def mutate_multicell(self):
        pass
    