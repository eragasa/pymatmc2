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

def cleanup():
    if os.path.isfile('pymatmc2.log'):
        os.remove('pymatmc2.log')
    if os.path.isfile('pymatmc2.config'):
        os.remove('pymatmc2.config')
    if os.path.isdir('simulations'):
        shutil.rmtree('simulations')


if __name__ == "__main__":
    create_configuration_file()
    kwargs_mc2 = {
        'configuration_path':'pymatmc2.config',
        'results_path':'pymatmc2.results',
        'logfile_path':'pymatmc2.log',
        'simulations_path':'simulations',
        'is_restart':False
    }
    o_mc2 = MultiCellMonteCarlo(**kwargs_mc2)
    o_mc2.create_iteration0_simulations()
    o_mc2.create_submission_scripts(i_iteration=0)
    assert isinstance(o_mc2, MultiCellMonteCarlo)
    assert os.path.isdir('simulations')
