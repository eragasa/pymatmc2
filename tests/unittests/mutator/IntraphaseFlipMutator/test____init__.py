import pytest
import os

from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell
from pymatmc2.mutator import BaseMultiCellMutator
from pymatmc2.mutator import IntraphaseFlipMutator

pymatmc2_resource_AuPt_path = os.path.join('resources', 'test____init__')
pymatmc2_configuration_AuPt_path = os.path.join(
    pymatmc2_configuration_AuPt_path,
    'pymatmc2.config'
)


def test__static_property():
    expected_mutator_type = 'intraphase_flip'
    assert IntraphaseFlipMutator.mutator_type == 'intraphase_flip'

def test__implements_base_class():
    for abstractmethod in BaseMultiCellMutator.__abstractmethods__:
        assert abstractmethod not in IntraphaseFlipMutator.__abstractmethods__

def test__constructor():
    configuration_path = 'pymatmc2.config'
    configuration = Pymatmc2Configuration()
    configuration.read(path=configuration_path)

    mutator = IntraphaseFlipMutator()
    assert issubclass(IntraphaseFlipMutator, MultiCellMutator)

def is_overridden_function(func):
    obj = func.__self__
    print_m = getattr(super(type(obj), obj), func.__name__)
    return func.__func__ != print_m.__func__

def dev__implements_base_class():
    for abstractmethod in BaseMultiCellMutator.__abstractmethods__:
        abstract_method_name = str(abstractmethod)
        print(abstract_method_name)
    
    for method in IntraphaseFlipMutator.__abstractmethods__:
        print(method)

if __name__ == "__main__":
    dev__implements_base_class()
    
