import pytest

import os
from pymatmc2 import Pymatmc2Configuration

@pytest.fixture
def pymatmc2_obj():
    obj = Pymatmc2Configuration()
    return obj

@pytest.fixture
def configuration():
    configuration = {}
    configuration['calculator'] = {
        'calculator':'vasp',
        'simulation_type':'vasp_min_all'
    }
    configuration['job_submission'] = {
        'hpc_type':'torque',
        'hpc_config_path':'osc.pitzer'
    }
    configuration['atomic_configuration'] = {
        'molar_fraction_total':{'Au':24, 'Pt':8},
        'simulation_cells':{
            'fcc1':{
                'incar':'INCAR',
                'kpoints':'KPOINTS',
                'poscar':'POSCAR_1'
            },
            'fcc2':{
                'incar':'INCAR',
                'kpoints':'KPOINTS',
                'poscar':'POSCAR_2'
            }
        }
    }
    configuration['environment_variables'] = {
        'temperatures':[400.0],
        'pressures':[0.0]
    }

    return configuration

def test_write(
    pymatmc2_obj: Pymatmc2Configuration, 
    tmpdir, 
    configuration: dict
):
    # set up the test
    configuration_path = os.path.join(str(tmpdir), 'pymatmc2.config')
    pymatmc2_obj.configure_from_dict(configuration)
    assert isinstance(pymatmc2_obj.configuration, dict)

    # code being tested
    pymatmc2_obj.write(path=configuration_path)

    # tests
    assert os.path.isfile(configuration_path)
    pymatmc2_obj.read(path=configuration_path)

if __name__ == "__main__":
    configuration = {}
    configuration['calculator'] = {
        'calculator':'vasp',
        'simulation_type':'vasp_min_all'
    }
    configuration['job_submission'] = {
        'hpc_type':'torque',
        'hpc_config_path':'osc.pitzer'
    }
    configuration['atomic_configuration'] = {
        'molar_fraction_total':{'Au':24, 'Pt':8},
        'simulation_cells':{
            'fcc1':{
                'incar':os.path.join('resources','INCAR_1'),
                'kpoints':os.path.join('resources','KPOINTS_1'),
                'poscar':os.path.join('resources','POSCAR_1'),
                'potcar':os.path.join('resources','POTCAR')
            },
            'fcc2':{
                'incar':os.path.join('resources','INCAR_1'),
                'kpoints':os.path.join('resources','KPOINTS_1'),
                'poscar':os.path.join('resources','POSCAR_1'),
                'potcar':os.path.join('resources','POTCAR')
            }
        }
    }
    configuration['environment_variables'] = {
        'temperatures':[400.0],
        'pressures':[0.0]
    }
    path = "pymatmc2.config"
    obj = Pymatmc2Configuration()
    obj.configure_from_dict(configuration)
    obj.write(path=path)
    assert os.path.isfile(path)

