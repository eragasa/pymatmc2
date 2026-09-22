# coding: utf-8
# Copyright (c) Eugene J. Ragasa
# Distributed under the terms of the MIT License

""" testing Pymatmc2Configuration constructors """ 

__author__ = "Eugene J. Ragasa"
__email__ = "ragasa.2@osu.edu"
__copyright__ = "Copyright 2020, Eugene J. Ragasa"
__maintainer__ = "Eugene J. Ragasa"
__date__ = "2020/02/22"

import pytest
from pymatmc2 import Pymatmc2Configuration

def test_default_constructor():
    obj = Pymatmc2Configuration()
    assert obj.path is None
    assert obj.configuration is None

if __name__ == "__main__":
    obj = Pymatmc2Configuration()
    assert obj.path is None
    assert obj.configuration is None
