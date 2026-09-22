import pytest
import os
from distutils import dir_util
from collections import OrderedDict
from mexm.io.vasp import Poscar
from mexm.simulation import VaspSimulation
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell

@pytest.fixture
def resource_path(request):
    file_path = request.module.__file__
    test_path = os.path.join(
        os.path.dirname(file_path),
        '..', # MultiCell dir
        '..', # unittests
        '..', # tests
        'resource', 'AgPt'
    )
    return test_path

def get_resource_path():
    file_path = os.path.abspath(__file__)
    resource_path = os.path.join(
        os.path.dirname(file_path),
        '..', # MultiCell dir
        '..', # unittests
        '..', # tests
        'resource', 'AgPt' 
    )
    return os.path.abspath(resource_path)

def get_configuration():
    config_path = os.path.join(get_resource_path(), 'input', 'pymatmc2.config')
    obj = Pymatmc2Configuration()
    obj.read(path=config_path)

    #monkeypatch the paths
    for _, v in obj.simulation_cells.items():
        for fn in  ['incar', 'poscar', 'potcar', 'kpoints']:
            v[fn] = os.path.join(get_resource_path(), v[fn])
    return obj

@pytest.fixture
def configuration(resource_path, testingdir):
    config_path = os.path.join(str(resource_path), 'pymatmc2.config')
    obj = Pymatmc2Configuration()
    obj.read(path=config_path)

    #monkeypatch the paths
    for _, v in obj.simulation_cells.items():
        for fn in  ['incar', 'poscar', 'potcar', 'kpoints']:
            v['fn'] = os.path.join(str(resource_path), v[fn.upper])
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

def dev_initialize_from_pymatmc2_configuration():
    print(80*'-')
    print('initialize_from_pymatmc2_configuration')
    print(80*'-')

    o_configuration = get_configuration()

    o_multicell = MultiCell.initialize_from_pymatmc2_configuration(
        configuration=o_configuration
    )
    print('cell_names:', o_multicell.cell_names)
    assert o_multicell.cell_names == ['fcc1', 'fcc2']
    for k in ['fcc1', 'fcc2']:
        assert isinstance(o_multicell.simulations[k], VaspSimulation)

    print('molar_fraction_total:')
    for k, v in o_multicell.molar_fraction_total.items():
        print("    {}: {}".format(k, v))


def dev_configure__vasp():
    print(80*'-')
    print('configure__vasp')
    print(80*'-')

    o_configuration = get_configuration()
    o_multicell = MultiCell()
    o_multicell.configure(configuration=o_configuration)

    print('cell_names:', o_multicell.cell_names)
    assert o_multicell.cell_names == ['fcc1', 'fcc2']
    for k in ['fcc1', 'fcc2']:
        assert isinstance(o_multicell.simulations[k], VaspSimulation)

    print('molar_fraction_total:')
    for k, v in o_multicell.molar_fraction_total.items():
        print("  {}: {}".format(k, v))

    print('cell_fraction_total:')
    for cell_name, obj_cell in o_multicell.cell_molar_fraction.items():
        print('  {}:'.format(cell_name))
        for symbol, percentage in obj_cell.items():
            print('    {}: {}'.format(symbol, percentage))

    print('phase_molar_fraction:')
    print(o_multicell.phase_molar_fraction)

if __name__ == "__main__":
    obj = MultiCell()
    resource_path = get_resource_path()
    print('resource_path:{}'.format(
        os.path.abspath(resource_path)
    ))

    dev_initialize_from_pymatmc2_configuration()
    dev_configure__vasp()