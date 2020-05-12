import pytest
import os
import pandas as pd
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import Pymatmc2Data

resources_path = os.path.join('resources', 'test____init__')
pymatmc2_config_path = os.path.join(resources_path, 'pymatmc2.config')
pymatmc2_data_path = os.path.join(resources_path, 'pymatmc2.results')

def test__default_constructor():
    data = Pymatmc2Data()

def test__property__configuration__set():
    configuration = Pymatmc2Configuration()
    configuration.read(path=pymatmc2_config_path)

    data = Pymatmc2Data()
    data.configuration = configuration

    assert isinstance(data.configuration, Pymatmc2Configuration)

def test__read():
    configuration = Pymatmc2Configuration()
    configuration.read(path=pymatmc2_config_path)

    data = Pymatmc2Data()
    data.configuration = configuration
    data.read(path=pymatmc2_data_path)

    assert isinstance(data.df, pd.DataFrame)
