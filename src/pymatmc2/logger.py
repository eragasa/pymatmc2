# coding: utf-8
# Copyright (c) Eugene J. Ragasa
# Distributed under the terms of the MIT License

""" the log file 

This module implements a logging facility
"""

__all__ = ["Pymatmc2Log"]
__author__ = "Eugene J. Ragasa"
__email__ = "ragasa.2@osu.edu"
__copyright__ = "Copyright 2020, Eugene J. Ragasa"
__maintainer__ = "Eugene J. Ragasa"
__date__ = "2020/02/22"

from datetime import datetime

class Pymatmc2Log:

    def __init__(self, path: str):
        self.path = path

    def log(self, message: str):
        logstr = '[{timestamp}] {message}\n'.format(
            timestamp=str(datetime.now()),
            message=message
        ) 
        with open(self.path, 'a') as f:
            f.write(logstr)
