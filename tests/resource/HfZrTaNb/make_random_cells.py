import pytest
import os
from mexm.io.vasp import Poscar
from pymatmc2 import Pymatmc2Configuration
from pymatmc2.multicell import get_bcc_cell
from pymatmc2.multicell import create_random_cell

def dev__make_random_cell__bcc__HfZrTaNb():
    configuration_path = os.path.join('resources', 'pymatmc2.config')
    cell_type = 'bcc'
    sc = [3,3,3]

    configuration = Pymatmc2Configuration()
    configuration.read(path=configuration_path)
    total_concentration = configuration.total_concentration
    
    cell = create_random_cell(cell_type=cell_type, composition=total_concentration, supercell=sc)

def dev__make_random_cell__bcc__HfZrTaNb__to_poscar():
    configuration_path = os.path.join('resources', 'pymatmc2.config')
    cell_type = 'bcc'
    sc = [3,3,3]

    configuration = Pymatmc2Configuration()
    configuration.read(path=configuration_path)
    total_concentration = configuration.total_concentration
    
    for phase_name in configuration.cell_names:
        cell = create_random_cell(cell_type=cell_type, composition=total_concentration, supercell=sc)
        o_poscar = Poscar.initialize_from_object(cell)
        poscar_path = os.path.join('resources', phase_name, 'POSCAR')
        print(poscar_path)
        o_poscar.write(poscar_path)

if __name__ == "__main__":
    dev__make_random_cell__bcc__HfZrTaNb__to_poscar()