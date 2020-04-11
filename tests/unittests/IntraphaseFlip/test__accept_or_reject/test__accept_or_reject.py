import pytest
import os
import shutil
from numpy import linalg
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell
from pymatmc2.multicellmutate import IntraphaseFlip

if __name__ == "__main__":
    configuration_path = os.path.join('resources','pymatmc2.config')
    configuration = Pymatmc2Configuration()
    configuration.read(path=configuration_path)

    multicell_initial_path = os.path.join('simulations','00000')
    multicell_initial = MultiCell()
    multicell_initial.configuration = configuration
    multicell_initial.read(path=multicell_initial_path)

    multicell_candidate_path = os.path.join('simulations','00001')
    multicell_candidate = MultiCell()
    multicell_candidate.configuration = configuration
    multicell_candidate.read(path=multicell_candidate_path)
    
    mutator = IntraphaseFlip()
    mutator.accept_or_reject(
        multicell_initial = multicell_initial,
        multicell_candidate = multicell_candidate,
        temperature = configuration.temperature
    )

