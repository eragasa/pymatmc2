import pytest
import os

from pymatmc2 import Pymatmc2Configuration
from pymatmc2.mutator import IntraphaseFlipMutator

resources_path = os.path.join('resources', 'test__mutate__min1_max1')
configuration_path = os.path.join(resources_path, 'pymatmc2.config')
configuration = Pymatmc2Configuration()
configuration.read(configuration_path)

mutator = IntraphaseFlipMutator()
mutator.configuration = configuration