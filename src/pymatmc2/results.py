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
import ase

from typing import List

from mexm.io.filesystem import OrderedDictYAMLLoader
from pymatmc2 import Pymatmc2Configuration







class Pymatmc2Results():
    """

    Attributes:
        results (dict): 
    """

    def __init__(self,
        configuration: Pymatmc2Configuration,
    ):
        self._configuration = None
        self.configuration = configuration

        self.tar_path = self.configuration.results_tar_path
        self.data_path = self.configuration.results_file_path


    @property
    def configuration(self) -> Pymatmc2Configuration:
        return self._configuration

    @configuration.setter
    def configuration(self, configuration: Pymatmc2Configuration):
        if not isinstance(configuration, Pymatmc2Configuration):
            msg = "configuration must be an instance of Pymatmc2Configuration"
            raise TypeError(msg)

        self._configuration = configuration

    @property
    def column_names(self) -> List[str]:
        """ :List(str)  - a list of the column names"""
        cell_names = self.configuration.cell_names
        symbols = self.configuration.symbols

        column_names = ['iteration', 'mutate_type']
        
        # initial concentrations
        for cell_name in cell_names:
            for symbol in symbols:
                name = 'initial.{}.{}'.format(cell_name, symbol)
                column_names += [name]

        # candidate concentrations
        for cell_name in cell_names:
            for symbol in symbols:
                name = 'candidate.{}.{}'.format(cell_name, symbol)
                column_names += [name]

        # final concentrations
        for cell_name in cell_names:
            for symbol in symbols:
                name = 'final.{}.{}'.format(cell_name, symbol)
                column_names += [name]

        # initial phase molar fractions
        for cell_name in cell_names:
            name = 'initial.{}.f'.format(cell_name)
            column_names += [name]

        #candidate phase molar fractions
        for cell_name in cell_names:
            name = 'candidate.{}.f'.format(cell_name)
            column_names += [name]

        #final phase molar fractions
        for cell_name in cell_names:
            name = 'final.{}.f'.format(cell_name)
            column_names += [name]

        # initial energies
        for cell_name in cell_names:
            name = 'initial.{}.E'.format(cell_name)
            column_names += [name]
        name = 'initial.total.E'
        column_names += [name]

        # candidate_energies
        for cell_name in cell_names:
            name = 'candidate.{}.E'.format(cell_name)
            column_names += [name]
        name = 'candidate.total.E'
        column_names += [name]

        # final energies
        for cell_name in cell_names:
            name = 'final.{}.E'.format(cell_name)
            column_names += [name]
        name = 'final.total.E'

        return column_names

    def get_iteration_path(self, i_iteration: int) -> str:
        iteration_path = os.path.join(
            self.configuration.results_path,
            self.configuration.phase_point_string,
            self.configuration.get_iteration_string(i=i_iteration)
        )
        return iteration_path

    def get_mc_path(self, i_iteration: int, mc_type: str):
        
        mc_types = ['initial', 'candidate', 'final']
        if mc_type not in mc_types:
            msg = "Invalid mc_type: {}".format(mc_type)
            raise ValueError(msg)

        mc_path = os.path.join(
            self.get_iteration_path(i_iteration),
            mc_type)
        
        return mc_path

    def get_mutate_type(self, i_iteration: int) -> str:
        iteration_path = self.get_iteration_path(i_iteration)
        mutate_type_path = os.path.join(iteration_path, 'mutate_type')
        with open(mutate_type_path) as f:
            mutate_type = f.read()
        return mutate_type

    def get_iteration_info(configuration: Pymatmc2Configuration, i_iteration: int):

        mc_types = ['initial', 'candidate', 'final']
        mc_paths = {}
        for mc_type in mc_types:
            mc_path[mc_type] = self.get_mc_path(i_iteration, mc_type)

        mc = {}
        for k,v in mc_paths.items():
            mc[k] = MultiCell()
            mc[k].configuration = configuration
            mc[k].read(path=mc_paths[k])

    def get_is_accepted_path(self, i_iteration: int) -> str:    
        iteration_path = self.get_iteration_path(i_iteration)
        is_accepted_path = os.path.join(iteration_path, 'is_accepted')
        return is_accepted_path
    
    def write_header_line(self):
        str_header = ",".join(self.column_names) + "\n"

        with open(self.data_path, 'w') as f:
            f.write(str_header)


    def archive_multicell(self, i_iteration: int, path: str, mc_type: str):
        if os.path.exists(self.tarball_path):
            tf = tarfile.open(self.tarball_path, mode='a')
        else:
            tf = tarfile.open(self.tarball_path, mode='w')

        if self.configuration.calculator_type == 'vasp':
            archived_files = [
                'POSCAR', 'INCAR', 'KPOINTS', 'OUTCAR', 'OSZICAR', 'CONTCAR', 'vasprun.xml'
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
