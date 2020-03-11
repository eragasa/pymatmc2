import pytest
import os
from distutils import dir_util
from collections import OrderedDict
from mexm.io.vasp import Poscar
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell

@pytest.fixture
def resource_path(request):
    test_path = os.path.join(
        os.path.dirname(request.module.__file__),
        '..', # MultiCell dir
        '..', # unittests
        '..', # tests
        'resource', 'AgPt', 'input'
    )
    return test_path

@pytest.fixture
def configuration(resource_path, testingdir):
    config_path = os.path.join(resource_path, 'pymatmc2.config')
    obj = Pymatmc2Configuration()
    obj.read(path=config_path)
    for _, v in obj.simulation_cells.items():
        v['poscar'] = os.path.join(str(testingdir), v['poscar'])
    return obj

@pytest.fixture
def testingdir(tmpdir, resource_path):
    tmpdir.mkdir("resource")
    assert os.path.isdir(
        os.path.join(
            str(tmpdir), 
            'resource'
        )
    )
    dir_util.copy_tree(
        resource_path, 
        os.path.join(str(tmpdir),'resource')
    )
    print(os.listdir(os.path.join(str(tmpdir),'resource')))
    return tmpdir

def test___init__():
    obj = MultiCell()
    assert isinstance(obj, MultiCell)

def test_initialize_from_pymatmc2_configuration(configuration, testingdir):
    obj = MultiCell.initialize_from_pymatmc2_configuration(
        configuration=configuration
    )
    assert isinstance(obj, MultiCell)

    for k in configuration.simulation_cells.keys():
        assert isinstance(obj.cells[k], Poscar)
