import pytest

from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import Pymatmc2Results

pymatmc2_config_path = 'pymatmc2.config'
pymatmc2_data_path = 'pymatmc2.data'

def test__default_constructor():
    configuration = Pymatmc2Configuration()
    configuration.read(path=pymatmc2_config_path)

    results = Pymatmc2Results()
    results.configuration = configuration

if __name__ == "__main__":
    pass
