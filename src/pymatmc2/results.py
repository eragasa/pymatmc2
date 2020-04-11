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
import tarfile

from mexm.io.filesystem import OrderedDictYAMLLoader
from pymatmc2 import Pymatmc2Configuration

class Pymatmc2Results():
    """

    Attributes:
        results (dict): 
    """

    def __init__(self,
        configuration: Pymatmc2Configuration,
        results_path='pymatmc2.results', 
        tarball_path='pymatmc2.tar'
    ):
        assert isinstance(results_path, str)
        assert isinstance(tarball_path, str)

        self.configuration = configuration
        self.results_path = results_path
        self.tarball_path = tarball_path

    def archive_multicell(self, i_iteration: int, path: str, mc_type: str):
        if os.path.exists(self.tarball_path):
            tf = tarfile.open(self.tarball_path, mode='a')
        else:
            tf = tarfile.open(self.tarball_path, mode='w')

        if self.configuration.calculator_type == 'vasp':
            archived_files = [
                'POSCAR', 'INCAR', 'KPOINTS', 'OUTCAR', 'OSZICAR', 'CONTCAR'
            ]

            for cell_name in self.configuration.cell_names:
                for f in archived_files:
                    src_path = os.path.join(path, f)
                    dst_path = os.path.join(
                        '{:05}'.format(i_iteration),
                        cell_name,
                        mc_type,
                        f
                    )        
                    tf.addfile(
                        tarfile.TarInfo(dst_path),
                        src_path
                    )
        tf.close()

    def get_multicell(self, i_iteration: int):
        pass

    def archive_initial_multicell(self, i_iteration: int, path: str):
        kwargs = {
            'i_iteration':i_iteration,
            'path':path,
            'mc_type':'initial'
        }
        self.archive_multicell(**kwargs)
    
    def archive_candidate_multicell(self, i_iteration: int, path: str):
        kwargs = {
            'i_iteration':i_iteration,
            'path':path,
            'mc_type':'candidate'
        }
        self.archive_multicell(**kwargs)

    def archive_final_multicell(self, i_iteration: int, path: str):
        kwargs = {
            'i_iteration':i_iteration,
            'path':path,
            'mc_type':'final'
        }
        self.archive_multicell(**kwargs)