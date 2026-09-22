""" module to manage jobs to slurm """
import os, yaml
from copy import copy, deepcopy
from collections import OrderedDict

from mexm.io.slurm.configuration import SlurmConfiguration
from mexm.io.slurm.submissionscript import SlurmSubmissionScript

class VaspSubmissionScript(SlurmSubmissionScript):
    """ class for creating SLURM submission scripts for VASP

    These classes do require that a SlurmConfiguration file have been already
    created, and provided as the path argument.
    """

    def __init__(self, path=None):
        SlurmSubmissionScript.__init__(self, path)
        self.modules = self.configuration.vasp_configuration['modules']
        self.run_string = self.configuration.vasp_configuration['run_string']

class LammpsSubmissionScript(SlurmSubmissionScript):
    """ class for creating SLURM submission scripts for LAMMPS

    These classes do require that a SlurmConfiguration file have been already
    created, and provided as the path argument.
    """

    def __init__(self, path=None):
        SlurmSubmissionScript.__init__(self, path)
        self.modules = self.configuration.lammps_configuration['modules']
        self.run_string = self.configuration.lammps_configuration['run_string']

class PhontSubmissionScript(SlurmSubmissionScript):
    """ class for creating SLURM submission scripts for PHONTS

    These classes do require that a SlurmConfiguration file have been already
    created, and provided as the path argument.
    """
    def __init__(self, path=None):
        SlurmSubmissionScript.__init__(self,path)
        self.modules = slef.configuration.phonts_configuration['modules']
        self.run_string = self.configuration.phonts_configuration['run_string']
