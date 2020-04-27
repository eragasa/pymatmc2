import os
import shutil
import subprocess
from time import sleep
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCellMonteCarlo

def get_src_input_directory():
    src_input_directory = os.path.join('..', '..', '..', 'resource', 'AgPt', 'input')
    return src_input_directory

def copy_resource_directory():
    src_input_directory = get_src_input_directory()
    dst_input_directory = "input"

    if os.path.isdir(dst_input_directory):
        shutil.rmtree(dst_input_directory)
    shutil.copytree(src=src_input_directory, dst=dst_input_directory)

def remove_resource_directory():
    shutil.rmtree('input')

def remove_pymatmc2_configuration():
    configuration_path = 'pymatmc2.config'
    os.remove('pymatmc2.config')

def create_pymatmc2_configuration():
    make_configuration_script_path = os.path.join(
        'input',
        'make_configuration.py'
    )
    
    if not os.path.isfile(make_configuration_script_path):
        sleep(10)

    cmd = 'python {}'.format(make_configuration_script_path)
    subprocess.Popen(cmd)

    assert os.path.isfile('pymatmc2.config')

def setup():
    print('copying the resource directory...')
    copy_resource_directory()
    print('creating the configuration file...')
    create_pymatmc2_configuration()

def cleanup():
    # remove_resource_directory()
    # remove_pymatmc2_configuration()

    shutil.rmtree('results')
    shutil.rmtree('simulations')

        
def dev____init____initial():
    kwargs_mc2 = {
        'configuration_path':os.path.join('input','pymatmc2.config'),
        'results_path':'results',
        'logfile_path':'pymatmc2.log',
        'simulations_path':'simulations',
        'is_restart':False
    }

    o_mc2 = MultiCellMonteCarlo(**kwargs_mc2)

    assert os.path.isdir('simulations')
    assert os.path.isdir(os.path.join('simulations', '400K_0GPa'))

    cleanup()

if __name__ == "__main__":
    dev____init____initial()