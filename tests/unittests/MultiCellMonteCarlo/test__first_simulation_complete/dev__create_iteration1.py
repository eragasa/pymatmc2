import os
import shutil
import subprocess
from pymatmc2 import MultiCellMonteCarlo

def create_configuration_file():
    make_configuration_path = os.path.join(
        '..', # ../MultiCellMonteCarlo
        '..', # ../unittests,
        '..', # ../tests
        'resource', # 
        'AgPt',
        'input',
        'make_configuration.py'
    )

    subprocess.Popen(
        'python {}'.format(make_configuration_path)
    )

def copy_simulation_results():
    simulations_path = os.path.join(
        '..', # ../MultiCellMonteCarlo
        '..', # ../unittests,
        '..', # ../tests
        'resource', # 
        'AgPt',
        'simulations'
    )

    if os.path.isdir('simulations'):
        shutil.rmtree('simulations')
    shutil.copytree(
        src=simulations_path,
        dst='simulations'
    )

def setup():
    create_configuration_file()
    copy_simulation_results()

if __name__ == "__main__":
    setup()
    assert os.path.isdir(
        os.path.join('simulations','400K_0GPa','00000')
    )
    assert os.path.isdir(
        os.path.join('simulations','400K_0GPa','00000','fcc1')
    )
    assert os.path.isdir(
        os.path.join('simulations','400K_0GPa','00000','fcc2')
    )
    kwargs_mc2 = {
        'configuration_path':'pymatmc2.config',
        'results_path':'pymatmc2.results',
        'logfile_path':'pymatmc2.log',
        'simulations_path':'simulations',
        'is_restart':True
    }
    o_mc2 = MultiCellMonteCarlo(**kwargs_mc2)
    o_mc2.run()

