from collections import OrderedDict
import os
from mexm.io.vasp import Poscar
from pymatmc2 import MultiCell
from pymatmc2 import Pymatmc2Configuration

def get_configuration(path: str) -> Pymatmc2Configuration:
    obj = Pymatmc2Configuration()
    obj.read(path=path)
    return obj

def dev_init():
    obj = MultiCell()

def test_initialize_from_pymatmc2_configuration():
    test_path = os.path.dirname(os.path.abspath(__file__))
    print(test_path)
    config_path = os.path.join(test_path, 'resources', 'pymatmc2.in')
    configuration = get_configuration(path=config_path)
    MultiCell.initialize_from_pymatmc2_configuation(
        configuration=configuration
    )

if __name__ == "__main__":
    pass