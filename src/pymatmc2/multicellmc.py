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

import os
import shutil
from mexm.io.vasp import Poscar
from mexm.simulation import VaspSimulation

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
        self.simulations_path = 'simulations'

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
        else:
            i_iteration = 0
            self.run_first_iteration()

    def run_first_iteration(self):
        self.log('starting iteration 0')
        i_iteration = 0
        n_cells = self.configuration.n_cells

        simulations = {}
        for k, v in self.configuration.simulation_cells.items():
            simulations[k] = VaspSimulation()
            simulations[k].incar.read(v['incar'])
            simulations[k].poscar.read(v['poscar'])
            simulations[k].kpoints.read(v['kpoints'])
            simulations[k].potcar.read(v['potcar'])
            simulation_path = os.path.join(
                self.simulations_path,
                '{}_{}_T{}_P{}'.format(i_iteration,k,)
            )
            simulations[k].write(path=simulation_path)

    def start_next_iterationse(self):
        pass

if __name__ == "__main__":
    mc2 = MultiCellMonteCarlo(is_restart=False)
    exit()

    k_B = BOLTZMANN
    C = COULOMB
    T = configuration.temperature
    N_cell = configuration.n_cells
    # prepare()

    utils.clear_lock_file()
    while True:
        stop_check()
        folders = get_structure('flip_alt')
        if len(folders) == 0:
            r1 = Results()
            r1.step += 1
            r1.add_results(REJECTED_FILE)
            continue
        for folder in folders:
            run_job(folder)
        
        # why am initializing results here
        r1 = Results()
        r2 = Results()
        r2.read_next()
        save_log('{:>5d} {}\n'.format(r2.step, time.strftime("%Y-%m-%d %H:%M")))

        
        if energy_new < energy_old:
            self.accept_new_configuration()
        else:
            self.calculate_rejection_probability(toten_1, toten_2)
            
        # what is this probability??
        probability = np.exp((r1.total_energy - r2.total_energy) * N_cell *
                             C / k_B / T)
        r2.probability = np.minimum(probability, 1.0)
        if np.random.rand() < r2.probability:
            r2.add_results(RESULTS_FILE)
            r2.tar_file()
            if r2.step >= 99999:
                save_log("Maximum step 99999 reached.")
                exit()
        else:
            r2.add_results(REJECTED_FILE)
