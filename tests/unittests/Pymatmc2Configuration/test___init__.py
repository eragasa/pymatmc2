import pytest
from pymatmc2 import Pymatmc2Configuration

def test_default_constructor():
    obj = Pymatmc2Configuration()
    assert obj.path is None
    assert obj.config_dict is None
