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

def dev_initialize_from_pymatmc2_configuration():
    test_path = os.path.dirname(os.path.abspath(__file__))
    print(test_path)
    config_path = os.path.join(test_path, 'resource', 'pymatmc2.config')
    configuration = get_configuration(path=config_path)
    obj = MultiCell.initialize_from_pymatmc2_configuration(
        configuration=configuration
    )
    return obj

def dev_property_symbols():
    multicell = dev_initialize_from_pymatmc2_configuration()
    assert isinstance(multicell.symbols, list)
    print(multicell.symbols)

def dev_property_cell_names():
    multicell = dev_initialize_from_pymatmc2_configuration()
    assert isinstance(multicell.cell_names, list)
    print(multicell.cell_names)

def dev_property_cell_molar_fraction():
    multicell = dev_initialize_from_pymatmc2_configuration()
    assert isinstance(multicell.cell_molar_fraction, list)
    print(multicell.cell_molar_fraction)

def dev_property_total_molar_fraction():
    multicell = dev_initialize_from_pymatmc2_configuration()
    assert isinstance(multicell.total_molar_fraction, list)
    print(multicell.total_molar_fraction)

def dev_property_phase_molar_fraction():
    multicell = dev_initialize_from_pymatmc2_configuration()
    assert isinstance(multicell.phase_molar_fraction, list)
    print(multicell.phase_molar_fraction)
if __name__ == "__main__":
    pass