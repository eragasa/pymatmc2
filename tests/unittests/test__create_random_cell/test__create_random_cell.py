import pytest

import numpy as np

from mexm.elements import ELEMENTS
from mexm.structure import SimulationCell
from mexm.io.vasp import Poscar

from tempsrc import create_random_cell

def test__create_random_cell__fcc():
    cell = create_random_cell(
        cell_type='fcc',
        composition={
            'Ag':0.50,
            'Pt':0.50
        },
        supercell=[3,3,3]
    )

    ELEMENTS['Ag'].atmrad
    ELEMENTS['Pt'].atmrad

    assert isinstance(cell, SimulationCell)

def dev__create_random_cell__fcc():
    kwargs = {
        'cell_type':'fcc',
        'composition':{
            'Ag':0.50,
            'Pt':0.50
        },
        'supercell':[3,3,3]
    }
    cell = create_random_cell(**kwargs)
    print('Ag:',cell.get_number_of_atoms('Ag'))
    print('Pt:',cell.get_number_of_atoms('Pt'))

    poscar = Poscar.initialize_from_object(obj=cell)
    poscar.write('POSCAR')

def dev__create_random_cell__bcc():
    kwargs = {
        'cell_type':'bcc',
        'composition':{
            'Ag':0.50,
            'Pt':0.50
        },
        'supercell':[3,3,3]
    }

    cell = create_random_cell(**kwargs)
    print('Ag:',cell.get_number_of_atoms('Ag'))
    print('Pt:',cell.get_number_of_atoms('Pt'))

    poscar = Poscar.initialize_from_object(obj=cell)
    poscar.write('POSCAR')

if __name__ == "__main__":
    dev__create_random_cell__fcc()
    # dev__create_random_cell__bcc()
