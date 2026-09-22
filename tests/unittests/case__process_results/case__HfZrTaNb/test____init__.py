import pytest
import os
from pymatmc2 import Pymatmc2Configuration
from plots import Pymatmc2Plots

def get_configuration_path() -> str:
    return 'pymatmc2.config'

def get_configuration() -> Pymatmc2Configuration:
    configuration_path = 'pymatmc2.config'
    o_configuration = Pymatmc2Configuration()
    o_configuration.read(path=configuration_path)
    return o_configuration

def get_data_path() -> str:
    return 'pymatmc2.data'

@pytest.fixture
def configuration_path() -> str:
    return get_configuration_path()

@pytest.fixture
def data_path() -> str:
    return get_data_path()

@pytest.fixture
def configuration() -> Pymatmc2Configuration:
    return get_configuration()

def test____init____no_args():
    o_plot = Pymatmc2Plots()
    assert isinstance(o_plot, Pymatmc2Plots)
    assert isinstance(o_plot.configuration, type(None))

def test____init____arg_configuration_as_obj(configuration):
    assert isinstance(configuration, Pymatmc2Configuration)

    o_plot = Pymatmc2Plots(configuration=configuration)
    assert isinstance(o_plot.configuration, Pymatmc2Configuration)

def test____init____arg_configuration_as_str(configuration_path):
    assert isinstance(configuration_path, str)
    assert os.path.isfile(configuration_path)

    plotter = Pymatmc2Plots(configuration=configuration_path)
    assert isinstance(plotter.configuration, Pymatmc2Configuration)

def test__property__configuration__w_obj(configuration):
    assert isinstance(configuration, Pymatmc2Configuration)

    plotter = Pymatmc2Plots()
    plotter.configuration = configuration
    assert isinstance(plotter.configuration, Pymatmc2Configuration)

def test__property__configuration__w_none():
    plotter = Pymatmc2Plots()
    plotter.configuration = None
    assert isinstance(plotter.configuration, type(None))

def test__property__configuration__w_str_valid_path(configuration_path: str):
    assert isinstance(configuration_path, str)
    assert os.path.isfile(configuration_path)

    plotter = Pymatmc2Plots()
    plotter.configuration = configuration_path
    assert isinstance(plotter.configuration, Pymatmc2Configuration)

def test__property__data__w_path(configuration, data_path):
    plotter = Pymatmc2Plots()
    plotter.configuration = configuration
    plot.data = data_path