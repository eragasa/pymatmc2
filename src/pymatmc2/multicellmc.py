# coding: utf-8
# Copyright (c) Eugene J. Ragasa
# Distributed under the terms of the MIT License

""" multicell montecarlo class

This module implements a logging facility
"""

__author__ = "Eugene J. Ragasa"
__email__ = "ragasa.2@osu.edu"
__copyright__ = "Copyright 2020, Eugene J. Ragasa"
__maintainer__ = "Eugene J. Ragasa"
__date__ = "2020/02/22"

import os
import shutil
from copy import deepcopy

from mexm.io.vasp import Poscar
from mexm.simulation import VaspSimulation
from mexm.job import JobSubmissionManagerFactory

from pymatmc2 import Pymatmc2Log
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import Pymatmc2Results
from pymatmc2.utils import clear_folders
from pymatmc2.utils import save_log
from pymatmc2.utils import get_ratio
from pymatmc2.utils import get_structure
from pymatmc2.utils import run_job
from pymatmc2.utils import stop_check
from pymatmc2.utils import prepare

# the results file needs to maintain state of the system


class MultiCellMonteCarlo():
    """ Implementation of the (MC)^2

    This class implements the algorithm as described in [1].
    
    Attributes:
        is_restart (bool)
        configuration (Pymatmc2Configuration)
        results (Pymatmc2ResultsFile)
        logfile (Pymatmc2LogFile)

    Notes:
    [1] C. Niu, Y. Rao, Wolfgang Windl, and Maryam Ghazisaeidi. "Multi-cell Monte Carlo method for phase prediction"
    """

    def __init__(
        self,
        configuration_path = 'pymatmc2.config',
        results_path = 'pymatmc2.results',
        logfile_path = 'pymatmc2.log',
        simulations_path = None,
        is_restart = True
    ):
        """

        Arguments:
            configuration_path (str): default: pymatmc2.config
            results_path (str): the path to the default: pypmatmc2.results
            logfile_path (str): default: pymatmc2.log
            simulations_path (str): default: `simulations`
            is_restart (bool): If set to true, then the code will assume that
                it as previously run before.  If set to False, then we will
                will delete previous results.
        """

        self.is_restart = is_restart
        self.simulations_path = simulations_path

        self.start_logging(
            path = logfile_path, 
            is_restart = is_restart
        )
    
        if self.is_restart:
            self.log('restarting pymatmc2')
        else:
            self.log('starting pymatmc2 from scratch.')

        self.start_configuration(
            path = configuration_path
        )
        self.create_simulation_directory(
            path = self.simulations_path,
            is_restart = self.is_restart
        )
        self.start_results(
            path = results_path,
            is_restart = self.is_restart
        )

        assert isinstance(self.configuration, Pymatmc2Configuration)
        assert isinstance(self.results, Pymatmc2Results)
        assert isinstance(self.logfile, Pymatmc2Log)
        assert os.path.isdir(self.simulations_path)
    
 
    def start_logging(self, path: str, is_restart: bool):
        """

        Arguments:
            path (str): path to the logging facility
        """
        if not self.is_restart:
            if os.path.isfile(path):
                os.remove(path)
        self.logfile = Pymatmc2Log(
            path=os.path.abspath(path)
        )


    def start_configuration(
        self,
        path: str
    ):
        """ read the configuration file

        Arguments:
            path (str): path to the configuration file
        """
        if isinstance(path, str):
            try:
                self.configuration = Pymatmc2Configuration()
                self.configuration.read(path=path)
                self.log('loaded configuration file')

            except FileNotFoundError:
                msg = 'cannot find configuration_path: {}'.format(
                    path
                )
                raise
        else:
            msg = "path must be a string"
            self.log(msg)
            raise TypeError(msg)

    def create_simulation_directory(
        self,
        path: str,
        is_restart: bool 
    ):        
        # delete simulation directory on new simulation
        if is_restart:
            pass
        else:
            # remove the existing simulation directory if it exists
            if os.path.isdir(path):
                self.log('removing existing simulation directory')
                shutil.rmtree(path)
            
            # create new simulation directory
            os.mkdir(path)

    def start_results(self, path: str, is_restart: bool):
        """ start the results file

        Arguments:
            path (str): path to the results file
            is_restart (bool): If set to true, it will ignore the 
                file there
        """

        if is_restart:
            self.results = Pymatmc2Results()
            self.results.read()
        else:
            if os.path.isfile(path):
                self.log('removing existing results file')
                os.remove(path)
            else:
                self.results = Pymatmc2Results()

    def log(self, message):
        self.logfile.log(message = message)
    
    def run(self):
        if self.is_restart:
            i_iteration = self.determine_current_iteration()
            print(i_iteration)
        else:
            i_iteration = 0
            self.run_first_iteration()

    def run_first_iteration(self):
        self.log('starting iteration 0')
        i_iteration = 0

        self.create_simulations(i_iteration=i_iteration)
        self.create_submission_scripts(i_iteration=i_iteration)
        self.submit_jobs(i_iteration=i_iteration)

    def create_simulations(self, i_iteration: int):
        # for vasp simulations
        simulations = {}
        simulation_names = []
        for k, v in self.configuration.simulation_cells.items():
            simulation_name = '{cellname}_{iteration:03}_T{temperature}_P{pressure}'.format(
                    cellname = k,
                    iteration = i_iteration,
                    temperature = int(self.configuration.temperature),
                    pressure = int(self.configuration.pressure)
                )
            
            simulation_names.append(simulation_name)
            simulations[k] = VaspSimulation()
            simulations[k].incar.read(v['incar'])
            simulations[k].poscar.read(v['poscar'])
            simulations[k].kpoints.read(v['kpoints'])
            simulations[k].potcar.read(v['potcar'])

            # abstract Simulation.write()
            simulation_path = os.path.join(self.simulations_path, simulation_name)
            print(simulation_path)
            os.mkdir(simulation_path)
            simulations[k].write(simulation_path=simulation_path)

    def create_submission_scripts(self, i_iteration: int):
        simulation_names = []
        for k, v in self.configuration.simulation_cells.items():
            simulation_name = '{cellname}_{iteration:03}_T{temperature}_P{pressure}'.format(
                    cellname = k,
                    iteration = i_iteration,
                    temperature = int(self.configuration.temperature),
                    pressure = int(self.configuration.pressure)
                )
            simulation_names.append(simulation_name)

            configuration_ = self.configuration.configuration
            hpc_type = configuration_['hpc_manager']['type']
            script_kwargs = deepcopy(configuration_['hpc_manager']['configuration'])
            script_kwargs['jobname'] = simulation_name
            script_path = os.path.join(
                self.simulations_path,
                simulation_name,
                'runjob.sh'
            )
            JobSubmissionManagerFactory.write_submission_script(
                hpc_type = hpc_type,
                script_kwargs = script_kwargs,
                script_path = script_path
            )

    def submit_jobs(self, i_iteration: int):
        simulation_name_fmt = '{cellname}_{iteration:03}_T{temperature}_P{pressure}'
        for k, v in self.configuration.simulation_cells.items():
            simulation_name = simulation_name_fmt.format(
                cellname = k,
                iteration = i_iteration,
                temperature = int(self.configuration.temperature),
                pressure = int(self.configuration.pressure)
            )
            simulation_path = os.path.join(
                self.simulations_path,
                simulation_name
            )

        JobSubmissionManagerFactory.submit_job(simulation_path=simulation_path)
    def determine_current_iteration(self) -> int:
        raise NotImplementedError

    def start_next_iteration(self):
        raise NotImplementedError

if __name__ == "__main__":
    mc2 = MultiCellMonteCarlo(is_restart=False)
    exit()