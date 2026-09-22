import pytest

import os
from pymatmc2 import Pymatmc2Configuration

@pytest.fixture
def pymatmc2_obj():
    obj = Pymatmc2Configuration()
    return obj

def test_path__relative(pymatmc2_obj):
    # setup of local path to test
    path = 'local'
    assert not os.path.isabs(path)

    # testing this code
    pymatmc2_obj.path = path

    # expected behavior
    assert os.path.isabs(pymatmc2_obj.path)
    assert pymatmc2_obj.path == os.path.abspath(path)

def test_path__absolute(pymatmc2_obj):
    # setup of absolute path to test
    path = os.path.abspath('local')
    assert os.path.isabs(path)

    # testing this code
    pymatmc2_obj.path = path
    # expected behavior
    assert os.path.isabs(pymatmc2_obj.path)
    assert pymatmc2_obj.path == path

if __name__ == "__file__":
    configuration = {
        'calculator':'vasp',
        'calculation_type':'vasp_min_all',
        'cells':[
            {
                'incar':'INCAR',
                'kpoints':'KPOINTS',
                'poscar':'POSCAR_1'
            },
            {
                'incar':'INCAR',
                'kpoints':'KPOINTS',
                'poscar':'POSCAR_2'
            }
        ],
        'temperatures': [400.0],
        'pressures': [0.0],
    }

