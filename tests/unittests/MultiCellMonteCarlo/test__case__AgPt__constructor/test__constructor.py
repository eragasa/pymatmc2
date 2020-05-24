import pytest
import os
import shutil
from time import sleep
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCellMonteCarlo

@pytest.fixture
def configuration_path():
    import platform
    if platform.system() == 'Linux':
        return os.path.join('resources', 'pymatmc2.config.linux')
    elif platform.system() == 'Mac':
        return os.path.join('resources', 'pymatmc2.config.mac')
    elif platform.system() == 'Windows':
        return os.path.join('resources', 'pymatmc2.config.windows')
    else:
        raise ValueError('platform.system() == {}, is not supported'.format(platform.system()))

@pytest.fixture
def kwargs_mc2(configuration_path):


    kwargs = {
        'configuration_path': configuration_path,
        'results_path':'results',
        'logfile_path':'pymatmc2.log',
        'simulations_path':'simulations',
        'is_restart':False
    }
    return kwargs

@pytest.fixture
def configuration(kwargs_mc2) -> Pymatmc2Configuration:
    configuration = Pymatmc2Configuration()
    configuration.read(path=kwargs_mc2['configuration_path'])
    return configuration

@pytest.fixture
def filesystem_fixture(scope='module'):
    if os.path.isdir('simulations'):
        shutil.rmtree('simulations')
    if os.path.isdir('results'):
        shutil.rmtree('results')

    yield filesystem_fixture

    shutil.rmtree('simulations')
    shutil.rmtree('results')

@pytest.fixture
def obj_mc2(kwargs_mc2):
    o_mc2 = MultiCellMonteCarlo(**kwargs_mc2)
    return o_mc2

def test__simulations_path_created(configuration, obj_mc2):
    simulations_path = 'simulations'
    assert os.path.isdir(simulations_path)

def test__simulation_phase_point_path_created(configuration, obj_mc2):
    simulations_phase_point_path = os.path.join('simulations', configuration.phase_point_string)
    assert os.path.isdir(simulations_phase_point_path)

def test__constructor__default(configuration, kwargs_mc2):

    o_mc2 = MultiCellMonteCarlo(**kwargs_mc2)
    print(configuration.phase_point_string)
    simulations_path = 'simulations'
    simulations_phase_point_path = os.path.join('simulations', '400K_0GPa')

    assert os.path.isdir(simulations_path)
    assert os.path.isdir(simulations_phase_point_path)

    # cleanup
    shutil.rmtree('simulations')
    shutil.rmtree('results')
