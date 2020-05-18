import pytest
import os
import shutil
from mexm.structure import SimulationCell
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell

def test__constructor__default():
    o_mc = MultiCell()

def test__constructor__from_configuration():
    src_resource_path = os.path.join('HfZrTaNb', 'resources')
    src_config_path = os.path.join(src_resource_path, 'pymatmc2.config')

    configuration = Pymatmc2Configuration()
    configuration.read(path=src_config_path)

    o_mc = MultiCell.initialize_from_configuration(configuration=configuration)
    for k, v in o_mc.cells.items():
        assert isinstance(k ,str)
        assert issubclass(type(v), SimulationCell)
def dev__constructor__from_configuration():
    src_resource_path = os.path.join('HfZrTaNb', 'resources')
    src_config_path = os.path.join(src_resource_path, 'pymatmc2.config')

    configuration = Pymatmc2Configuration()
    configuration.read(path=src_config_path)

    o_mc = MultiCell.initialize_from_configuration(configuration=configuration)
    for k, v in o_mc.cells.items():
        print(k, type(v))

def test__property__symbols__getters():
    src_resource_path = os.path.join('HfZrTaNb', 'resources')
    src_config_path = os.path.join(src_resource_path, 'pymatmc2.config')

    configuration = Pymatmc2Configuration()
    configuration.read(path=src_config_path)

    o_mc = MultiCell.initialize_from_configuration(configuration=configuration)
    assert isinstance(o_mc.symbols, list)

def dev__property__symbols__getters():
    src_resource_path = os.path.join('HfZrTaNb', 'resources')
    src_config_path = os.path.join(src_resource_path, 'pymatmc2.config')

    configuration = Pymatmc2Configuration()
    configuration.read(path=src_config_path)

    o_mc = MultiCell.initialize_from_configuration(configuration=configuration)
    assert isinstance(o_mc.symbols, list)
    print(o_mc.symbols)

def test__property__cell_names():
    expected_cell_names = ['bcc1', 'bcc2', 'bcc3', 'bcc4']
    src_resource_path = os.path.join('HfZrTaNb', 'resources')
    src_config_path = os.path.join(src_resource_path, 'pymatmc2.config')

    configuration = Pymatmc2Configuration()
    configuration.read(path=src_config_path)

    o_mc = MultiCell.initialize_from_configuration(configuration=configuration)
    assert isinstance(o_mc.cell_names, list)
    assert o_mc.cell_names == expected_cell_names

def dev__property__cell_names():
    print(80*'-')
    print('property__cell_names')
    print(80*'-')
    src_resource_path = os.path.join('HfZrTaNb', 'resources')
    src_config_path = os.path.join(src_resource_path, 'pymatmc2.config')

    configuration = Pymatmc2Configuration()
    configuration.read(path=src_config_path)

    o_mc = MultiCell.initialize_from_configuration(configuration=configuration)
    assert isinstance(o_mc.cell_names, list)
    print(o_mc.cell_names)


def dev__singular_matrix_detection__HfZrTaNb():
    src_resource_path = os.path.join('HfZrTaNb', 'resources')
    src_config_path = os.path.join(src_resource_path, 'pymatmc2.config')

    configuration = Pymatmc2Configuration()
    configuration.read(path=src_config_path)

    o_mc = MultiCell.initialize_from_configuration(configuration=configuration)
    
    from numpy import linalg
    rank = linalg.matrix_rank(o_mc.cell_concentration_matrix)
    n_rows, n_cols = o_mc.cell_concentration_matrix.shape
    print('rank:{}'.format(rank))
    print('nrows:{}, ncols:{}'.format(n_rows, n_cols))

if __name__ == "__main__":
    #dev__constructor__from_configuration()
    #dev__property__symbols__getters()
    #dev__property__cell_names()
    dev__singular_matrix_detection__HfZrTaNb()
