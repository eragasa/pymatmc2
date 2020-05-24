# coding: utf-8
# Copyright (c) Eugene J. Ragasa
# Distributed under the terms of the MIT License

""" the module implements the MultiCellMutator base class """

__author__ = "Eugene J. Ragasa"
__email__ = "ragasa.2@osu.edu"
__copyright__ = "Copyright 2020, Eugene J. Ragasa"
__maintainer__ = "Eugene J. Ragasa"
__date__ = "2020/04/20"

from typing import Tuple
from abc import ABC
from abc import abstractmethod
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell


class BaseMultiCellMutator(ABC):
    mutate_type = 'base_mc_mutator'

    """ Implements an abstract base class for MC composition mutation type """
    def __init__(self):
        self._configuration = None
        self._multicell_initial = None
        self._multicell_candidate = None
        self._multicell_final = None
        self._is_accept = None

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
        if not isinstance(mc, MultiCell):
            msg = 'multicell_initial must be type pymatmc2.MultiCell'
            raise TypeError(msg)
        self._multicell_initial = mc

    @property
    def multicell_candidate(self) -> MultiCell:
        return self._multicell_candidate

    @multicell_candidate.setter
    def multicell_candidate(self, mc: MultiCell):
        if not isinstance(mc, MultiCell):
            msg = 'multicell_candidate must be type pymatmc2.MultiCell'
            raise TypeError(msg)
        self._multicell_candidate = mc

    @property
    def multicell_final(self) -> MultiCell:
        return self._multicell_final

    @multicell_final.setter
    def multicell_final(self, mc: MultiCell):
        if not isinstance(mc, MultiCell):
            msg = 'multicell_final must be type pymatmc2.MultiCell'
            raise TypeError(msg)

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
        temperature: float,
        pressure: float
    ) -> Tuple[bool, MultiCell]:
        raise 
