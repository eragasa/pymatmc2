import pytest
import os
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell
from pymatmc2.mutator import MultiCellMutatorFactory

configuration_path = 'pymatmc2.config'
mc_initial_path = os.path.join('simulations', '800K_0GPa', '00000')

if __name__ == "__main__":
    configuration = Pymatmc2Configuration()
    configuration.read(path=configuration_path)

    mc_initial = MultiCell()
    mc_initial.configuration = configuration
    mc_initial.read(path=mc_initial_path)

    o = MultiCellMutatorFactory()
    o.configuration = configuration
    o.determine_mutate_algorithm()
    mutate_type, mc_candidate = o.mutate_multicell(multicell=mc_initial)
