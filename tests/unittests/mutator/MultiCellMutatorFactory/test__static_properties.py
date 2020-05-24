import pytest
from pymatmc2.mutator import BaseMultiCellMutator
from pymatmc2.mutator import MultiCellMutatorFactory

expected_mutator_types = ['intraphase_flip', 'intraphase_swap']

def test__factories__expected_factories():
    assert isinstance(MultiCellMutatorFactory.factories, dict)
    assert len(expected_mutator_types) == len(MultiCellMutatorFactory.factories)
    for mutator_type in expected_mutator_types:
        assert mutator_type in MultiCellMutatorFactory.factories
        assert issubclass(MultiCellMutatorFactory.factories[mutator_type], BaseMultiCellMutator)

def dev__factories():
    print(MultiCellMutatorFactory.factories)
    for mutator_type in expected_mutator_types:
        implements_base_class = issubclass(MultiCellMutatorFactory.factories[mutator_type], BaseMultiCellMutator)
        print(
            'mutator_type:{}'.format(mutator_type), 
            MultiCellMutatorFactory.factories[mutator_type].__name__,
            'implements_base_class:{}'.format(implements_base_class)
        )

if __name__ == "__main__":
    dev__factories()