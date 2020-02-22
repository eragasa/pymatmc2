import os
import subprocess
import shutil
import argparse
import tarfile
import time
import numpy as np
import pdb
import glob
from scipy.optimize import nnls

BOLTZMANN   = 1.38064852E-23
COULOMB     = 1.60217662E-19
INPUT_FILE  = 'mc2.in'
RESULTS_FILE = 'Results'
REJECTED_FILE = 'Rejected'
TAR_FILE = 'Files.tar'
LOG_FILE = 'Logs'
STOP_FILE = 'stop'
flip_flag = 0   # global dynamic variable for flip determination

# MU_AUAU = -0.5446
# MU_AUPT = -0.7586
# MU_PTPT = -1.0102

from vasp import Vasp
from results import Results

from utils import clear_folders
from utils import save_log
from utils import get_ratio
from utils import get_structure
from utils import run_job
from utils import stop_check
from utils import prepare

class LogFile:

    def __init__(self, path):
        self.path = path

    def log(self, message): 
        print(message)

class Pymatmc2Configuration():
    """ configuration file for pymatmc2

    'Attributes:'
    'temperature (float): MC2 temperature'
    'flip_step (int): remove because not necessary'
    'iniital_c (float): initial concentrations'
    'n_cell (int): number of cells'
    """
    def __init__(self):
        self.mpi_vasp = None
        self.temperature = None
        self.flip_step = None
        self.initial_c = None
        pass

    def read(self, path):
        self.path = path

        try:
            with open(self.path, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise

        for line in lines:
            line = line.split('#', 1)[0].strip().split('=')

            if line[0].strip() == 'MPI_VASP':
                self.mpi_vasp = line[1].strip()
            elif line[0].strip() == 'TEMPERATURE':
                self.temperature = float(line[1].strip())
            elif line[0].strip() == 'FLIP_STEP':
                self.flip_step = int(line[1].strip())
            elif line[0].strip() == 'INITIAL_C':
                INITIAL_C = list(line[1].strip().split(' '))
                self.initial_c = [float(i) for i in INITIAL_C]
            elif line[0].strip() == 'N_CELL':
                self.n_cells = int(line[1].strip())

    def __repr__(self):
        return {
            'temperature':self.temperature,
            'concentration_0':self.initial_c,
            'n_cells':self.n_cells
        }

    def __str__(self):
        s = (
            'Pymatmc2('
            'temperature={temperature}, '
            'concentration_0={concentration}, '
            '
        )
if __name__ == "__main__":
    pymatmc2_in_path = 'mc2.in'
    pymatmc2_log_path = 'mc2.log'

    #initialize logfile
    logfile = LogFile(path=pymatmc2_log_path)
    
    configuration = Pymatmc2Configuration()
    configuration.read(path=pymatmc2_in_path)

    try:
        with open(pymatmc2_in_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:

    if not os.path.exists(INPUT_FILE):
        save_log("No mc2.in file found!")
        exit()
    with open(INPUT_FILE) as f:
        lines = f.readlines()
    if len(lines) == 0:
        save_log("The file mc2.in is empty!")
        exit()



    prepare()
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
        r1, r2 = Results(), Results()
        r2.read_next()
        save_log('{:>5d} {}\n'.format(r2.step, time.strftime("%Y-%m-%d %H:%M")))
        probability = np.exp((r1.total_energy - r2.total_energy) * N_CELL *
                             COULOMB / BOLTZMANN / TEMPERATURE)
        r2.probability = np.minimum(probability, 1.0)
        if np.random.rand() < r2.probability:
            r2.add_results(RESULTS_FILE)
            r2.tar_file()
            if r2.step >= 99999:
                save_log("Maximum step 99999 reached.")
                exit()
        else:
            r2.add_results(REJECTED_FILE)
