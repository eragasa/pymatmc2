from abc import ABC, abstractmethod

from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell


class MultiCellMutateAlgorithm(ABC):
    def __init__(self):
        self.multicell = None
        self._configuration = None
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
