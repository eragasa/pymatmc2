import pytest

import os
import shutil

from numpy import linalg

# abstract classes
from mexm.structure import SimulationCell
from mexm.simulation import AtomicSimulation

# concrete classes
from mexm.io.vasp import Poscar
from mexm.simulation import VaspSimulation

from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell

expected_cell_names = ['bcc1', 'bcc2', 'bcc3', 'bcc4']
expected_symbols = ['Hf', 'Zr', 'Ta', 'Nb']
expected_structure_concrete_class = Poscar
expected_simulation_concrete_class = VaspSimulation

@pytest.fixture
def src_resource_path():
    return os.path.join('HfZrTaNb', 'resources')

@pytest.fixture
def src_config_path(src_resource_path):
    return os.path.join(src_resource_path, 'pymatmc2.config')

@pytest.fixture
def configuration(src_config_path) -> Pymatmc2Configuration:
    configuration_ = Pymatmc2Configuration()
    configuration_.read(path=src_config_path)
    return configuration_

@pytest.fixture
def obj_multicell(configuration) -> MultiCell:
    o_mc = MultiCell.initialize_from_configuration(configuration=configuration)
    return o_mc

def test__constructor__from_configuration(obj_multicell):
    assert isinstance(obj_multicell.configuration, Pymatmc2Configuration)
    assert len(obj_multicell.cells) == len(expected_cell_names)
    for k in expected_cell_names:
        assert k in obj_multicell.cells

def test__property_cells__is_simulationcell_subclass(obj_multicell):
    for cell_name, obj_cell in obj_multicell.cells.items():
        assert issubclass(type(obj_cell), SimulationCell)

def test__property_cells__is_expected_concrete_class(obj_multicell):
    for cell_name, obj_cell in obj_multicell.cells.items():
        assert isinstance(obj_cell, expected_structure_concrete_class)

def test__property_simulations__is_simulation_subclass(obj_multicell):
    """ check to see if concrete class is a subclass of mexm.simulation.Simulations """
    for cell_name, obj_simulation in obj_multicell.simulations.items():
        assert issubclass(type(obj_simulation), AtomicSimulation)

def test__property_simulations__is_expected_concrete_class(obj_multicell):
    """ check to see if concrete class is an instance of mexm.simulation.VaspSimulation """
    for cell_name, obj_simulation in obj_multicell.simulations.items():
        assert isinstance(obj_simulation, VaspSimulation)  

def test__property__symbols__getters(obj_multicell):
    assert isinstance(obj_multicell.symbols, list)
    assert len(obj_multicell.symbols) == len(expected_symbols)
    for symbol in expected_symbols:
        assert symbol in obj_multicell.symbols

def test__property__cell_names__getters(obj_multicell):
    assert isinstance(obj_multicell.cell_names, list)
    assert obj_multicell.cell_names == expected_cell_names

def dev__property__cell_names():
    src_resource_path_ = os.path.join('HfZrTaNb', 'resources')
    src_config_path_ = os.path.join(src_resource_path_, 'pymatmc2.config')
    
    configuration_ = Pymatmc2Configuration()
    configuration_.read(path=src_config_path_)

def dev__constructor__from_configuration():
    src_resource_path_ = os.path.join('HfZrTaNb', 'resources')
    src_config_path_ = os.path.join(src_resource_path_, 'pymatmc2.config')

    configuration_ = Pymatmc2Configuration()
    configuration_.read(path=src_config_path_)

    obj_multicell_ = MultiCell.initialize_from_configuration(configuration=configuration_)

    print(80*'-')
    print('MultiCell.initialize_from_configuration()')
    print('src_config_path:{}'.format(src_resource_path_))
    for k, v in obj_multicell_.cells.items():
        print('type(Multicell.cells[{}]={}'.format(k, type(v)))
    for k, v in obj_multicell_.simulations.items():
        print(k, v)

if __name__ == "__main__":
    dev__constructor__from_configuration()
