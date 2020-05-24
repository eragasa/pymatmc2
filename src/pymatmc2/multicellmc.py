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
import random
from typing import Dict, Tuple, List
from copy import deepcopy
from collections import OrderedDict
from numpy import linalg

from mexm.io.vasp import Poscar
from mexm.simulation import VaspSimulation
from mexm.job import JobSubmissionManagerFactory

from pymatmc2 import MultiCell
from pymatmc2 import Pymatmc2Log
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import Pymatmc2Results
# from pymatmc2 import Pymatmc2Tarball
from pymatmc2 import Pymatmc2Data

from pymatmc2.utils import clear_folders
from pymatmc2.utils import save_log
from pymatmc2.utils import get_ratio
from pymatmc2.utils import get_structure
from pymatmc2.utils import run_job
from pymatmc2.utils import stop_check
from pymatmc2.utils import prepare

# the results file needs to maintain state of the system
from pymatmc2.mutator import MultiCellMutatorFactory

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
        results_path = 'results',
        logfile_path = 'pymatmc2.log',
        simulations_path = 'simulations',
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
        self.results_path = results_path

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

        self.create_simulations_directory(
            path = self.simulations_path,
            is_restart = self.is_restart
        )

        self.create_results_directory(
            path = self.results_path,
            is_restart = self.is_restart
        )

    @property
    def phase_space_name(self) -> str:
        name = self.get_phase_space_name(
            temperature = self.configuration.temperature,
            pressure = self.configuration.pressure
        )

        return name

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
                msg = 'cannot find configuration_path: {}'
                msg = msg.format(path)
                raise
        else:
            msg = "path must be a string"
            self.log(msg)
            raise TypeError(msg)

    def create_simulations_directory(
        self,
        path: str,
        is_restart: bool 
    ):
        self.simulations_path = path
        # delete simulation directory on new simulation
        if is_restart:
            pass
        else:
            # remove the existing simulations directory if it exists
            if not os.path.isdir(path):
                msg = 'create_simulations directory: {}'
                msg = msg.format(self.simulations_path)
                self.log(msg)
                os.mkdir(path)
            
            simulation_phasespace_path = os.path.join(
                self.simulations_path,
                self.phase_space_name
            )

            if os.path.isdir(simulation_phasespace_path):
                msg = 'removing existing phasespace directory: {}'
                msg = msg.format(simulation_phasespace_path)
                self.log(msg)
                shutil.rmtree(simulation_phasespace_path)
            msg = 'create phasespace directory: {}'
            msg = msg.format(simulation_phasespace_path)
            os.mkdir(simulation_phasespace_path)

    def create_results_directory(self, 
        path: str, 
        is_restart: bool
    ):
        if is_restart:
            pass
        else:
            if not os.path.isdir(path):
                msg = 'create results_directory:{}'
                msg = msg.format(path)
                self.log(msg)
                os.mkdir(path)
 
            results_phasespace_path = os.path.join(
                self.results_path,
                self.phase_space_name
            )

            if os.path.isdir(results_phasespace_path):
                msg = 'removing existing results directory: {}'
                msg = msg.format(results_phasespace_path)
                self.log(msg)
                shutil.rmtree(results_phasespace_path)
            os.mkdir(results_phasespace_path)

    def log(self, message: str):
        self.logfile.log(message = message)
    
    def run(self) -> bool:
        """

        Returns:
            bool: True if the max iteration condition has been met
        """
        is_max_iterations = False
        if not self.is_restart:
            self.run_iteration(i_iteration=0)
        else:
            i_iteration, status = self.determine_current_iteration()

            if status == 'running':
                self.log('iteration {} is still runnning'.format(i_iteration))
                exit()
            else:
                if i_iteration > self.configuration.max_iterations:
                    self.log('maximum iterations reached')
                    is_max_iterations = True
                
                elif i_iteration >= 0:
                    next_iteration = i_iteration + 1
                    self.process_iteration(i_iteration=i_iteration)
                    self.run_iteration(i_iteration=next_iteration)

                else:
                    msg = "how are we at iteration {}".format(i_iteration)
                    raise ValueError(msg)

        return is_max_iterations
    
    def run_iteration(self, i_iteration: int):
        msg = 'starting iteration {}'.format(i_iteration)
        self.log(msg)
        msg = 'creating simulations...'
        self.log(msg)
        self.create_simulations(i_iteration=i_iteration)
        msg = 'submitting jobs...'
        self.log(msg)
        self.create_submission_scripts(i_iteration=i_iteration)
        self.submit_jobs(i_iteration=i_iteration)
    
    def process_iteration(self, i_iteration: int):
        """

        Arguments:
            i_iteration (int): the iteration of the simulations to process
        """

        # the initial configuration for the current iteration is the
        # final configuration of the previous one
        src_mc_initial_path = os.path.join(
            self.results_path,
            self.phase_space_name,
            self.get_iteration_string(i_iteration-1),
            'final'
        )

        # the candidate configuration are the simulation which have just
        # been completed
        src_mc_candidate_path = os.path.join(
            self.simulations_path,
            self.phase_space_name,
            self.get_iteration_string(i_iteration)
        )

        src_mutate_type = os.path.join(
            self.results_path,
            self.phase_space_name,
            self.get_iteration_string(i_iteration),
            'mutate_type'
        )

        dst_path = os.path.join(
            self.results_path,
            self.phase_space_name,
            self.get_iteration_string(i_iteration)
        )

        dst_mc_initial_path = os.path.join(
            self.results_path,
            self.phase_space_name,
            self.get_iteration_string(i_iteration),
            'initial'
        )
        dst_mc_candidate_path = os.path.join(
            self.results_path,
            self.phase_space_name,
            self.get_iteration_string(i_iteration),
            'candidate'
        )
        dst_mc_final_path = os.path.join(
            self.results_path,
            self.phase_space_name,
            self.get_iteration_string(i_iteration),
            'final'
        )
        dst_is_accepted = os.path.join(
            self.results_path,
            self.phase_space_name,
            self.get_iteration_string(i_iteration),
            'is_accepted'
        )

        if i_iteration == 0:
            # no initial_mc
            mc_candidate = MultiCell()
            mc_candidate.configuration = self.configuration
            mc_candidate.read(path=src_mc_candidate_path)

            #no accept or reject
            mc_final = mc_candidate

            if not os.path.isdir(dst_path):
                os.mkdir(dst_path)
            mc_final.archive(dst_path = dst_mc_final_path)

        else:

            # read the initial cell
            mc_initial = MultiCell()
            mc_initial.configuration = self.configuration
            mc_initial.read(path=src_mc_initial_path)

            # read the candidate cell
            mc_candidate = MultiCell()
            mc_candidate.configuration = self.configuration
            mc_candidate.read(path=src_mc_candidate_path)

            # read the mutation type
            with open(src_mutate_type, 'r') as f:
                mutate_type = f.read()
          
            temperature = self.configuration.temperature
            pressure = self.configuration.pressure
            mutator = MultiCellMutateAlgorithmFactory.factories[mutate_type]()
            is_accept, mc_final = mutator.accept_or_reject(
                multicell_initial = mc_initial,
                multicell_candidate = mc_candidate,
                temperature = temperature,
                pressure = pressure
            )

            mc_initial.archive(dst_path=dst_mc_initial_path)
            mc_candidate.archive(dst_path=dst_mc_candidate_path)
            mc_final.archive(dst_path=dst_mc_final_path)
            with open(dst_is_accepted, 'w') as f:
                f.write(str(is_accept))


    def create_simulations(self, i_iteration:int):
        if i_iteration == 0:
            multicell = MultiCell.initialize_from_configuration(
                self.configuration
            )
                
            dst_path = os.path.join(
                self.simulations_path, 
                self.phase_space_name, 
                self.get_iteration_string(i_iteration)
            )

            multicell.write(path = dst_path)
            msg = 'creating simulations in {}'.format(dst_path)
            self.log(msg)
            #try:
            #    multicell.phase_molar_fraction
            #    multicell.write(path = dst_path)
            #except linalg.LinAlgError:
            #    msg = 'concentration matrix is not full rank'
            #    self.log(msg)
            #    msg = 'using intraphase_flip to fix concentration matrix'
            #    self.log(msg)

            #    mutator = IntraphaseFlip()
            #    mutator.configuration = self.configuration
            #    multicell_new = mutator.mutate_multicell(multicell = multicell)

            #    msg = 'found valid multicell'
            #    self.log(msg)
            #    multicell_new.write(path = dst_path)

            #    results_path = os.path.join(
            #        self.results_path,
            #        self.phase_space_name,
            #        self.get_iteration_string(i_iteration)
            #    )
            #    os.mkdir(results_path)

            #    mutate_type_path = os.path.join(results_path, 'mutate_type')
            #    with open(mutate_type_path,'w') as f:
            #        f.write(IntraphaseFlip.mutate_type)

        else:
            src_path = os.path.join(
                self.results_path,
                self.phase_space_name,
                self.get_iteration_string(i_iteration-1),
                'final'
            )

            dst_path = os.path.join(
                self.simulations_path,
                self.phase_space_name,
                self.get_iteration_string(i_iteration)
            )

            # read initial multicell
            mc_initial = MultiCell()
            mc_initial.configuration = self.configuration
            mc_initial.read(path=src_path)

            # create candidate multicell         
            mutator = MultiCellMutateAlgorithmFactory()
            mutator.configure(configuration=self.configuration)
            mutator.configuration = self.configuration
            assert isinstance(mutator.configuration, Pymatmc2Configuration)
            mutate_type, mc_candidate = mutator.mutate_cells(mc_initial)
            mc_candidate.write(path=dst_path)

            results_path = os.path.join(
                self.results_path,
                self.phase_space_name,
                self.get_iteration_string(i_iteration)
            )
            os.mkdir(results_path)

            mutate_type_path = os.path.join(results_path, 'mutate_type')
            with open(mutate_type_path,'w') as f:
                f.write(mutate_type)

    
    def create_submission_scripts(self, i_iteration: int):

        for k in self.configuration.simulation_cells:
            simulation_name = '{:05}_{}'.format(i_iteration, k)

            hpc_type = self.configuration.hpc_manager['type']
            script_kwargs = deepcopy(
                self.configuration.hpc_manager['configuration']
            )
            script_kwargs['jobname'] = simulation_name

            script_path = os.path.join(
                self.simulations_path,
                self.phase_space_name,
                self.get_iteration_string(i_iteration),
                k,
                'runjob.sh'
            )
            JobSubmissionManagerFactory.write_submission_script(
                hpc_type = hpc_type,
                script_kwargs = script_kwargs,
                script_path = script_path
            )

    def submit_jobs(self, i_iteration: int):
        
        hpc_type = self.configuration.hpc_manager['type']
        for k in self.configuration.simulation_cells:

            simulation_path = os.path.join(
                self.simulations_path,
                self.phase_space_name,
                self.get_iteration_string(i_iteration),
                k
            )

            JobSubmissionManagerFactory.submit_job(
                hpc_type=hpc_type,
                simulation_path=simulation_path,
                submission_script_path='runjob.sh'
            )


    def determine_current_iteration(self) -> Tuple[int, str]:
        """ return the currrent iteration and the status of the current iteration

        Returns:
            Tuple[int, str]:
                int: the current iteration
                str: status of the current iteration
        """
        i_iteration = 0
        status = None

        # determining the max iteration
        iteration_paths = os.listdir(
            path=os.path.join(
                self.simulations_path,
                self.phase_space_name
            )
        )
        i_iteration = max([int(k) for k in iteration_paths])

        # determing is the jobs are complete
        jobs_completed_array = []
        for cell_name in self.configuration.cell_names:
            job_complete_path = os.path.join(
                self.simulations_path,
                self.phase_space_name,
                self.get_iteration_string(i_iteration),
                cell_name,
                'jobComplete'
            )
            jobs_completed_array.append(os.path.isfile(job_complete_path))
        if all(jobs_completed_array):
            status = 'complete'
        else:
            status= 'running'

        return i_iteration, status

    def get_phase_space_name(self, 
        temperature: float, 
        pressure: float
    ) -> str:

        T = int(temperature)
        P = int(pressure)

        fmt = '{}K_{}GPa'
        return fmt.format(T, P)

    def get_iteration_string(self, i: int) -> str:
        """
        Arguments:
            i (int): number of the iteration of interest
        Returns:
            (str): formatted iteration string
        """
        fmt = '{:05}'
        return fmt.format(i)

if __name__ == "__main__":
    mc2 = MultiCellMonteCarlo(is_restart=False)
    exit()
