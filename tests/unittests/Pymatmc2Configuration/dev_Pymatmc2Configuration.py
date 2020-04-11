import os
import subprocess
from pymatmc2 import Pymatmc2Configuration

def get_config_path():
    config_path = os.path.join('resources','pymatmc2.config')
    return config_path

def create_pymatmc2_config_file():
    original_path = os.path.abspath(os.getcwd())
    os.chdir('resources')
    subprocess.Popen('python make_configuration.py')
    os.chdir(original_path)

def get_obj_configuration():
    config_path = get_config_path()
    configuration = Pymatmc2Configuration()
    configuration.read(path=config_path)
    return configuration

if __name__ == "__main__":
    configuration = get_obj_configuration()