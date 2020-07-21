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

from utils import Poscar
from utils import save_log

class Results(object):
    """Dealing with MC^2 results.
    """

    def __init__(self):
        """Read current results."""
        p = Poscar('OSCAR')
        #init_structs = glob.glob('POSCAR_*')
        self.total_cells = p.total_elements
        self.step = 0
        self.probability = 1.0
        self.compositions = np.zeros((N_CELL, self.total_cells))
        #self.potentials = np.zeros((self.total_cells, self.total_cells))
        self.cell_energy = np.zeros(N_CELL)
        #self.diff_energy = np.zeros(self.total_cells)
        self.ratio = np.zeros(N_CELL)
        self.total_energy = 0.0
        with open(RESULTS_FILE) as f:
            lines = f.readlines()
        if len(lines) == 1:
            self.read_next()
        else:
            self.read_results()

    def read_results(self):
        """Read results from RESULTS_FILE."""
        with open(RESULTS_FILE) as f:
            last = f.readlines()[-1].split()
        self.step = int(last[0])
        for i in range(N_CELL):
            for j in range(self.total_cells):
                k = i * (self.total_cells + 1) + j + 1
                self.compositions[i, j] = float(last[k])
            k = (i + 1) * (self.total_cells + 1)
            self.cell_energy[i] = float(last[k])
        self.ratio = np.array([
            float(x) / 100 for x in last[-N_CELL - 2:-2]
        ])
        self.total_energy, self.probability = [float(x) for x in last[-2:]]

    def read_next(self):
        """Read results from TMP directories."""
        self.step += 1
        for i in range(N_CELL):
            p = Poscar('TMP-{}/POSCAR'.format(i + 1))
            self.compositions[i, :] = p.element_counts
            #self.potentials[i, :], e_fit = p.chemical_potential(3.1, 12)
            file_name = 'TMP-{}/OSZICAR'.format(i + 1)
            if os.path.exists(file_name):
                with open(file_name) as f:
                    last = f.readlines()[-1]
                self.cell_energy[i] = float(last.split()[4])
            #self.diff_energy[i] = e_fit - self.cell_energy[i]
        #p = Poscar('POSCAR')
        #arr = np.array(inital_c)
        arr = np.array([initial_c])
        arr = np.concatenate((self.compositions.T, arr.T), axis=1)
        print(self.compositions)
        print(self.cell_energy)
        print(self.ratio)
        self.ratio = get_ratio(arr)
        self.total_energy = np.sum(self.cell_energy * self.ratio)

    def add_results(self, file_name):
        """Add current results to given file."""
        with open(file_name, 'a+') as f:
            f.write('{:>5d}'.format(self.step))
            for i in range(N_CELL):
                for j in range(self.total_cells):
                    f.write('{:3.0f}'.format(self.compositions[i, j]))
                # for j in range(self.total_cells):
                #     f.write('{:6.2f}'.format(self.potentials[i, j]))
                f.write('{:10.3f}'.format(self.cell_energy[i]))
                # f.write('{:7.3f}'.format(self.diff_energy[i]))
                # f.write('{:10.3f}'.format(self.star_energy[i]))
            for i in range(N_CELL):
                f.write('{:7.2f}'.format(self.ratio[i] * 100))
            f.write('{:10.3f}{:5.2f}\n'.format(self.total_energy, 
                self.probability))

    def tar_file(self):
        """Add files to TAR_FILE"""
        with tarfile.open(TAR_FILE, 'a') as tar:
            for i in range(N_CELL):
                folder = 'TMP-{}'.format(i + 1)
                wanted = ['POSCAR', 'CONTCAR', 'OSZICAR', 'vasp.out']
                for j in wanted:
                    src = '{}/{}'.format(folder, j)
                    if os.path.exists(src):
                        dst = '{}_{}_{}'.format(j, self.step, i + 1)
                        tar.add(src, arcname=dst, recursive=False)
                if not os.path.exists('TMP'):
                    os.mkdir('TMP')
                shutil.copy('{}/CONTCAR'.format(folder),
                            'TMP/POSCAR_{}'.format(i + 1))


def clear_folders():
    """Clear all TMP folders."""
    p = Poscar('OSCAR')
    for i in range(N_CELL):
        folder = 'TMP-{}'.format(i + 1)
        if os.path.exists(folder):
            for j in os.listdir(folder):
                os.remove('{}/{}'.format(folder, j))





def get_ratio(arr):
    """Calculate cell ratios from (n+1)*n input."""
    n = N_CELL
    a = arr[:, :n]
    b = arr[:, n]
    x = nnls(a, b)[0]
    c = np.dot(a, x)
    if all(-1E-8 <= i <= 1+1E-8 for i in x) and all(abs(i) < 1E-6 for i in b-c):
        return abs(x / np.sum(x))
    else:
        return np.zeros(n)


def get_structure(mode):
    """Get structures for next step."""
    global flip_flag
    clear_folders()
    if mode == 'flip_all':
        r = Results()
        folders = []
        for i in range(N_CELL):
            p = Poscar('TMP/POSCAR_{}'.format(i + 1))
            if p.flip():
                folders.append('TMP-{}'.format(i + 1))
            p.output('TMP-{}'.format(i + 1), i + 1)
        flip_flag = 0
        if len(folders) != 0:
            r.read_next()
            if np.sum(r.ratio) == 0:
                folders = []
        return folders
    elif mode == 'flip_one':
        r = Results()
        folders = []
        random_value = np.random.rand()
        sum_value = .0
        for i in range(N_CELL):
            poscar_path = os.path.join('TMP', 'POSCAR_{}'.format(i+1))
            p = Poscar('TMP/POSCAR_{}'.format(i + 1))
            if sum_value < random_value <= sum_value + r.ratio[i] and r.ratio[i] > 0.01:
                if p.flip():
                    folders = ['TMP-{}'.format(i + 1)]
            p.output('TMP-{}'.format(i + 1), i + 1)
            sum_value += r.ratio[i]
        flip_flag += 1
        if len(folders) != 0:
            r.read_next()
            if np.sum(r.ratio) == 0:
                folders = []
        return folders
    elif mode == 'flip_alt':
        if flip_flag == FLIP_STEP:
            return get_structure('flip_all')
        else:
            return get_structure('flip_one')


def run_job(folder):
    """Run VASP under given folder."""
    if not os.path.exists(folder):
        save_log('Cannot find folder {} for VASP.'.format(folder))
        exit()

    #create a context directory
    cwd = os.getcwd()
    os.chdir('{}/{}'.format(cwd, folder))
    os.system('{} > vasp.out'.format(MPI_VASP))

    #return to the context directory
    os.chdir(cwd)
    if not os.path.exists('{}/OSZICAR'.format(folder)):
        save_log('No OSZICAR under {} after running VASP.'.format(folder))
        exit()
    with open('{}/OSZICAR'.format(folder)) as f:
        last = f.readlines()[-1]
    if last.split()[1] != 'F=':
        save_log('VASP under {} is not finished.'.format(folder))
        exit()


def stop_check():
    """Check for manual stop."""
    if os.path.exists(STOP_FILE):
        os.remove(STOP_FILE)
        save_log("Found stop file. Stopping...")
        exit()


def prepare():
    """Make sure all files exist."""
    wanted = ['INCAR', 'KPOINTS', 'POTCAR']
    for i in wanted:
        if not os.path.exists(i):
            save_log('{} not found.'.format(i))
            exit()
    with open(RESULTS_FILE, 'a+') as f:
        lines = f.readlines()
    if len(lines) > 1:
        return
    unwanted = [STOP_FILE, TAR_FILE, LOG_FILE]
    for i in unwanted:
        if os.path.exists(i):
            os.remove(i)
    clear_folders()
    line = '{:>{}}'.format('Step', 5)
    o = Poscar('OSCAR')
    for i in range(N_CELL):
        for j in range(o.total_elements):
            line += '{:>3}'.format(o.element_names[j])
        line += '{:>10}'.format('E_DFT')
    for i in range(N_CELL):
        line += '{:>7}'.format('x_{}%'.format(i + 1))
    line += '{:>10}{:>5}\n'.format('E({:.0f}K)'.format(TEMPERATURE), 'P')
    with open(RESULTS_FILE, 'w') as f:
        f.write(line)
    with open(REJECTED_FILE, 'w') as f:
        f.write(line)
    for i in range(N_CELL):
        p = Poscar('POSCAR_{}'.format(i+1))
        folder = 'TMP-{}'.format(i + 1)
        p.output(folder, i+1)
        run_job(folder)
    r = Results()
    r.add_results(RESULTS_FILE)
    r.tar_file()

#------- this is the beginning of the script
if not os.path.exists(INPUT_FILE):
    save_log("No mc2.in file found!")
    exit()
with open(INPUT_FILE) as f:
    lines = f.readlines()
if len(lines) == 0:
    save_log("The file mc2.in is empty!")
    exit()
for line in lines:
    line = line.split('#', 1)[0].strip().split('=')
    if line[0].strip() == 'MPI_VASP':
        MPI_VASP = line[1].strip()
    elif line[0].strip() == 'TEMPERATURE':
        TEMPERATURE = float(line[1].strip())
    elif line[0].strip() == 'FLIP_STEP':
        FLIP_STEP = int(line[1].strip())
    elif line[0].strip() == 'INITIAL_C':
        INITIAL_C = list(line[1].strip().split(' '))
        initial_c = [float(i) for i in INITIAL_C]
    elif line[0].strip() == 'N_CELL':
        N_CELL = int(line[1].strip())



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
