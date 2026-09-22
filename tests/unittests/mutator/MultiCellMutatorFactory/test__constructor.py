import pytest
import os
import shutil
from pymatmc2 import MultiCell
from pymatmc2 import Pymatmc2Configuration
from pymatmc2.mutator import IntraphaseSwapMutator
from pymatmc2.mutator import IntraphaseFlipMutator
from pymatmc2.mutator import MultiCellMutatorFactory

src_resource_path = os.path.join('resources', 'test__constructor')
src_configuration_path = os.path.join(src_resource_path, 'pymatmc2.config')

dst_configuration_path = 'pymatmc2.config'

def setup():
    shutil.copy(
        src=src_configuration_path,
        dst=dst_configuration_path
    )

def test__default_constructor():
    mutator = MultiCellMutatorFactory()

def test__property__configuration__setter():
    configuration_path = dst_configuration_path
    configuration = Pymatmc2Configuration()
    configuration.read(path=configuration_path)

    mutator = MultiCellMutatorFactory()
    mutator.configuration = configuration

def test__property__configuration__getter():
    configuration_path = dst_configuration_path
    configuration = Pymatmc2Configuration()
    configuration.read(path=configuration_path)

    mutator = MultiCellMutatorFactory()
    mutator.configuration = configuration

    assert isinstance(mutator.configuration, Pymatmc2Configuration)

def test__property__multicell_algorithms__getter():
    expected_mutator_dict = {
        'intraphase_swap': IntraphaseSwapMutator,
        'intraphase_flip': IntraphaseFlipMutator
    }
    mutator = MultiCellMutatorFactory()
    mutator_dict = mutator.multicell_mutators

    for k, v in expected_mutator_dict.items():
        assert k in mutator_dict
        assert mutator_dict[k] == v

def test__property__multicell_algorithms__get_intraphase_swap():
    mutator_factory = MultiCellMutatorFactory()
    mutator = mutator_factory.multicell_mutators['intraphase_swap']()

    assert isinstance(mutator, IntraphaseSwapMutator)

def test__property__multicell_algorithms__get_intraphase_flip():
    mutator_factory = MultiCellMutatorFactory()
    mutator = mutator_factory.multicell_mutators['intraphase_flip']()

    assert isinstance(mutator, IntraphaseFlipMutator)
