# coding: utf-8
# Copyright (c) Eugene J. Ragasa
# Distributed under the terms of the MIT License

"""
This module implements the marshalling and unmarshalling of results to the filesystem
"""
__all__ = ["Pymatmc2Results"]
__author__ = "Eugene J. Ragasa"
__email__ = "ragasa.2@osu.edu"
__copyright__ = "Copyright 2020, Eugene J. Ragasa"
__maintainer__ = "Eugene J. Ragasa"
__date__ = "2020/02/22"

import os
import copy
import yaml
from mexm.io.filesystem import OrderedDictYAMLLoader

class Pymatmc2Results():
    """

    Attributes:
        results (dict): 
    """

    def __init__(self):
        self._path = None
        self.results = None

    @property
    def path(self):
        """(str): absolute path of the results file"""
        
        return self._path

    @path.setter
    def path(self, path):

        # convert to absolute
        if isinstance(path, str):
            if os.path.isabs(path):
                self._path = path
            else:
                self._path = os.path.abspath(path)
        else:
            msg = "path must be a string"
            raise TypeError(msg)


    def read(self, path: str):
        """ read the results file

        Arguments:
            path (str): path to the results file
        """
        self.path
        raise NotImplementedError

    def write(self):
        raise NotImplementedError