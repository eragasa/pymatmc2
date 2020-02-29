import pytest
from collections import OrderedDict
from mexm.io.vasp import Poscar
from pymatmc2.multicell import MultiCell

def test__MultiCell__init():
    o = MultiCell()

def test__MultiCell__read_poscar_files():
    o = MultiCell(n_cells = 2)
    o.cells[0].read('POSCAR_1')
    o.cells[1].read('POSCAR_2')

    assert isinstance(o.cells[0], Poscar)
    assert isinstance(o.cells[1], Poscar)

def test__MultiCell__symbols():
    """ testing symbols property """
    o = MultiCell(n_cells = 2)
    o.cells[0].read('POSCAR_1')
    o.cells[1].read('POSCAR_2')

    assert o.symbols == list(set(o.cells[0].symbols +  o.cells[1].symbols))

def test__MultiCell__get_number_of_atoms():
    o = MultiCell(n_cells=2)

    o.cells[0].read('POSCAR_1')
    o.cells[1].read('POSCAR_2')

    for symbol in o.symbols:
        assert o.get_number_of_atoms(symbol) == o.cells[0].get_number_of_atoms(symbol) + o.cells[1].get_number_of_atoms(symbol)
    
    assert o.get_number_of_atoms() == o.cells[0].get_number_of_atoms() + o.cells[1].get_number_of_atoms()
