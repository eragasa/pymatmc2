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



class Poscar(object):
    """Dealing with POSCAR from VASP."""

    def __init__(self, file_name):
        """"Read POSCAR file."""
        with open('OSCAR') as f:
            lines = f.readlines()
        element_names = lines[5].split()
        with open(file_name) as f:
            lines = f.readlines()
        self.header = lines[0]
        self.lattice_constant = float(lines[1].strip())
        self.lattice_vectors = np.zeros((3, 3))
        for i in range(3):
            self.lattice_vectors[i, :] = \
                [float(x) for x in lines[2 + i].split()]
        self.element_names = element_names
        self.element_counts = np.zeros(len(element_names), dtype=int)
        this_element_names = lines[5].split()
        this_element_counts = lines[6].split()
        matched = []
        for i in range(len(this_element_names)):
            if this_element_names[i] not in element_names:
                save_log("Current POSCAR has elements not found in ./POSCAR")
                exit()
            for j in range(len(element_names)):
                if this_element_names[i] == element_names[j]:
                    self.element_counts[j] = int(this_element_counts[i])
                    matched.append(j)
                    break
        self.total_elements = len(this_element_names)
        self.total_atoms = sum(self.element_counts)
        self.coordinate_mode = lines[7]
        self.atom_coordinates = np.zeros((self.total_atoms, 3))
        for i in range(self.total_atoms):
            self.atom_coordinates[i, :] = \
                [float(x) for x in lines[8 + i].split()[0:4]]
        self.atom_types_index = np.zeros(self.total_atoms, dtype=int)
        count = 0
        for i in range(len(this_element_names)):
            for j in range(int(this_element_counts[i])):
                self.atom_types_index[count] = matched[i]
                count += 1
        sort = np.argsort(self.atom_types_index)
        self.atom_coordinates = self.atom_coordinates[sort, :]
        self.atom_types_index = self.atom_types_index[sort]

    def output(self, path, cell):
        """Write POSCAR and POTCAR to given path."""
        if not os.path.exists(path):
            os.mkdir(path)
        matched = [x for x, item in enumerate(self.element_counts) if item > 0]
        with open(path + '/POSCAR', 'w') as f:
            f.write(self.header)
            f.write(str(self.lattice_constant) + '\n')
            for i in range(3):
                f.write('{:14.8f}{:14.8f}{:14.8f}\n'.format(
                    *self.lattice_vectors[i, :]))
            f.write(' '.join([str(self.element_names[x])
                              for x in matched]) + '\n')
            f.write(' '.join([str(self.element_counts[x])
                              for x in matched]) + '\n')
            f.write(self.coordinate_mode)
            for i in range(self.total_atoms):
                f.write('{:14.8f}{:14.8f}{:14.8f}\n'.format(
                    *self.atom_coordinates[i, :]))
        with open('POTCAR') as f:
            lines = f.readlines()
        line_nums = [0]
        line_nums += [x + 1 for x, item in enumerate(lines)
                      if 'End of Dataset' in item]
        with open(path + '/POTCAR', 'w') as f:
            for i in matched:
                for j in range(line_nums[i], line_nums[i + 1]):
                    f.write(lines[j])
        shutil.copy('KPOINTS_{}'.format(cell), '{}/KPOINTS'.format(path))
        shutil.copy('INCAR_{}'.format(cell), '{}/INCAR'.format(path))
        shutil.copy('{}/POSCAR'.format(path), '{}/CONTCAR'.format(path))

    def flip(self):
        """"Flip a random atom in this POSCAR."""
        p = Poscar('OSCAR')
        (r1, r2) = np.random.randint(32, size=2)
        elem1 = self.atom_types_index[r1]
        elem2 = np.random.randint(2)
        if elem1 == elem2:
            return False
        else:
            self.element_counts[elem1] -= 1
            self.element_counts[elem2] += 1
            self.atom_types_index[r1] = elem2
            sort = np.argsort(self.atom_types_index)
            self.atom_coordinates = self.atom_coordinates[sort, :]
            self.atom_types_index = self.atom_types_index[sort]
            return True

    def neighbor(self, r1, r2):
        """Calculate number of pairs between r1/r2."""
        if not (self.coordinate_mode[0] == 'd' or
                self.coordinate_mode[0] == 'D'):
            save_log("POSCAR must be in direct mode.")
            exit()
        total_bonds_per_atom = np.zeros(self.total_atoms, dtype=int)
        bond_lengths = np.zeros((self.total_atoms, 50), dtype=float)
        atomic_index = np.zeros((self.total_atoms, 50), dtype=int)
        element_index = np.zeros((self.total_atoms, 50), dtype=int)
        number_bonds = np.zeros((self.total_atoms, self.total_elements),
                                dtype=int)
        total_bonds = np.zeros((len(self.element_names), self.total_elements),
                               dtype=int)
        for i in range(self.total_atoms):
            xyz1 = np.dot(self.atom_coordinates[i, :], self.lattice_vectors)
            for j in range(self.total_atoms):
                for k in range(-2, 3):
                    for m in range(-2, 3):
                        for n in range(-2, 3):
                            xyz2 = np.dot(self.atom_coordinates[j, :] +
                                          [k, m, n], self.lattice_vectors)
                            dist = np.sum(np.power(xyz1 - xyz2, 2))
                            dist = np.power(dist, .5) * self.lattice_constant
                            if r1 < dist < r2:
                                bond_lengths[i, total_bonds_per_atom[i]] = dist
                                atomic_index[i, total_bonds_per_atom[i]] = j
                                element_index[i, total_bonds_per_atom[i]] = \
                                    self.atom_types_index[j]
                                total_bonds_per_atom[i] += 1
        for i in range(self.total_atoms):
            for j in range(total_bonds_per_atom[i]):
                number_bonds[i, element_index[i, j]] += 1
        count = 0
        for i in range(len(self.element_names)):
            for j in range(len(self.element_names)):
                for k in range(self.element_counts[i]):
                    total_bonds[i, j] += number_bonds[k + count, j]
            count += self.element_counts[i]
        total_bonds /= 2
        return total_bonds, total_bonds_per_atom

def make_vasp_simulation(
        path,
        src_poscar=None,
        src_incar=None,
        src_potcar=None,
        src_kpoints=None):

    dst_poscar_path = "{}/{}".format(path, 'POSCAR')
    dst_incar_path = "{}/{}".format(path, 'INCAR')
    dst_potcar_path = "{}/{}".format(path, 'POTCAR')
    dst_kpoints_path = "{}/{}".format(path, 'KPOINTS')
    
    if isinstance(src_poscar, str):
        shutil.copy(src_poscar, dst_poscar_path)
    elif isinstance(src_poscar, Poscar):
        src_poscar.write(path=dst_poscar_path)
    else:
        raise TypeError('src_poscar must either a Poscar object or a path')

    if instance(src_incar, str):
        shutil.copy(src_incar, dst_incar_path)
    elif isinstance(src_poscar, Incar):
        Incar.write(path=dst_incar_path)
    else:
        raise TypeError('src_incar must either an Incar object or a path')

    if isinstance(src_potcar, str):
        shutil.copy(src_potcar, dst_potcar_path)
    elif isinstance(src_potcar, Potcar):
        src_potcar.write(path=dst_potcar_path)
    else:
        raise TypeError('src_potcar must either be a Potcar object or a path')

    if isinstance(src_kpoints, str):
        shutil.copy(src_kpoints, dst_kpoints-path)
    elif isinstance(src_kpoints, Kpoints):
        src_kpoints.write(path=dst_kpoints_path)
    else:
        raise TypeError('src_kpoints must either be a Kpoints object or a path')

class Results(object):
    """Dealing with MC^2 results."""



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
        self.ratio = np.array([float(x) / 100
                               for x in last[-N_CELL - 2:-2]])
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




def save_log(line):
    """Save a line to the LOG_FILE."""
    with open(LOG_FILE, 'a+') as f:
        f.write(line + '\n')


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
    """Get structures for next step.
    
    
    
    """
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
    cwd = os.getcwd()
    os.chdir('{}/{}'.format(cwd, folder))
    os.system('{} > vasp.out'.format(MPI_VASP))
    os.chdir(cwd)
    if not os.path.exists('{}/OSZICAR'.format(folder)):
        save_log('No OSZICAR under {} after running VASP.'.format(folder))
        exit()
    with open('{}/OSZICAR'.format(folder)) as f:
        last = f.readlines()[-1]
    if last.split()[1] != 'F=':
        save_log('VASP under {} is not finished.'.format(folder))
        exit()





class MultiCellMonteCarloConfiguration:

    def __init__(self):
        self.vasp_mpi_bin = None
        self.temperature = None
        self.flip_step = None
        self.cells = []

    def add_cell(self, path):
        self.cells.append(path)

class MultiCellMonteCarlo:
    calculator_types = ['vasp', 'lammps']
    simulation_types = ['npt', 'nvt', 'static']

    def __init__(self):
        """
        Properties:
                
        """

        # calculator paths
        self.vasp_std_bin_ = None
        self.lammps_serial_bin_ = None
        self.lammps_mpi_bin_ = None

        self.temperature_ = None
        self.flip_step_ = None
        self.n_cells = None
        self.initial_concentration = None
        self.calculator_type_ = 'vasp'
        self.simulation_type_ = 'static'

        self.results_path = 'results.out'

    @property
    def calculator_type(self):
        return self.calculator_type_

    @calculator_type.setter
    def calculator_type(self, calculator_type):
        calculator_type_ = calculator_type.lower()
        if calculator_type_ not in MultiCellMonteCarlo.calculator_types:
            msg = 'unsupported calculator_type: {}'.format(calculator_type)
            raise ValueError(msg)
        self.calculator_type_ = calculator_type

    @property
    def vasp_std_bin(self):
        return self.vasp_std_bin_

    @vasp_std_bin.setter
    def vasp_std_bin(self, vasp_std_bin):
        if not os.path.isfile(vasp_std_bin):
            msg = 'the path to vasp_std_bin is invalid: {}'.format(vasp_std_bin)
            raise FileNotFoundError(msg)
        self.vasp_std_bin_ = vasp_std_bin

    @property
    def lammps_serial_bin(self):
        return self.lammps_serial_bin_

    @lammps_serial_bin.setter
    def lammps_serial_bin(self, lammps_serial_bin):
        if not os.path.isfile(lammps_serial_bin):
            msg = 'the path to lammps_serial_bin is invalid: {}'.format(lammps_serial_bin)
            raise FileNotFoundError(msg)
        self.lammps_serial_bin_ = lammps_serial_bin
    @property
    def simulation_type(self):
        return self.simulation_type_

    @simulation_type.setter
    def simulation_type(self, simulation_type):
        simulation_type_ = simulation_type.lower()
        if simulation_type_ not in MultiCellMonteCarlo.simulation_types:
            msg = 'unsupported simulation_type: {}'.format(simulation_type_)
            raise ValueError(msg)
        self.simulation_type_ = simulation_type_

    @property
    def temperature(self):
        return self.temperature_

    @temperature.setter
    def temperature(self, temperature):
        temperature_ = float(temperature)
        if temperature_ < 0:
            msg = 'temperature must be greater than 0'
            raise ValueError(msg)
        self.temperature_ = temperature_


    def run(self, config_path='mc2.in'):
        """
        Arguments:
            config_path (str): path to the multicell montecarlo path
        """
        self.read_configuration(path=config_path)
        self.prepare_simulation()

    def read_configuration(self, path='mc2.in'):
        with open(path, 'r') as f:
            lines = f.readlines()
        
        for line in lines:
            line = line.split('#', 1)[0].strip().split('=')
            config_tag = line[0].strip().lower()
            config_value = line[1].strip()

            if config_tag == 'calculator_type':
                self.calculator_type = config_value
            elif config_tag == 'simulation_type':
                self.simulation_type = config_value
            elif config_tag == 'vasp_std_bin':
                self.vasp_std_bin = config_value
            elif config_tag == 'temperature':
                self.temperature = float(config_value)
            elif config_tag == 'flip_step':
                self.flip_step = int(config_value)
            elif config_tag == 'initial_c':
                self.initial_c = [float(k) for k in config_value.split()]
            elif config_tag == 'n_cells':
                self.n_cells = int(config_value)
            else:
                raise ValueError('unknown tag: {}'.format(config_tag))

        # ensure that we can find the vasp binaries required
        if self.calculator_type == 'vasp':
            if self.simulation_type == 'npt':
                raise NotImplementedError()
            elif self.simulation_type == 'nvt':
                raise NotImplementedError()
            elif self.simulation_type == 'static':
                if self.vasp_std_bin is None:
                    self.vasp_std_bin = os.environ['VASP_STD_BIN']
            else:
                raise ValueError('{} is an invalid option for simulation_type')

        # ensure that we can find the lammps binaries are required
        if self.calculator_type == 'lammps':
            if self.simulation_type in ['npt', 'nvt']:
                if self.lammps_mpi_bin is None:
                    self.vasp_mpi_bin  = os.environ['LAMMPS_MPI_BIN']
            elif self.simulation_type in ['static']:
                if self.lammps_serial_bin is None:
                    self.vasp_mpi_bin = os.environ['LAMMPS_SERIAL_BIN']
            else:
                raise ValueError('{} is an invalid option for simulation_type')

    def log(self, msg):
        save_log(msg)

    def remove_old_output_files(self):
        output_file_paths = [
                self.stop_path,
                self.tar_path,
                self.log_path]
        for path in output_file_paths:
            if os.path.isfile(path):
                os.remove(path)

    def remove_temp_folders():
        """Clear all TMP folders."""
        p = Poscar('OSCAR')
        for i in range(self.n_cells):
            path = 'tmp-{}'.format(i)
            if os.path.exists(path):
                shutil.rmtree(path)
    
    def prepare_simulations(self):
        """Make sure all files exist."""
        
        required_files = ['INCAR', 'KPOINTS', 'POTCAR']
        for file_path in required_files:
            if not os.path.isfile(file_path):
                msg = 'cannot find the file, {}'.format(file_path)
                self.log(msg)
                raise ValueError(msg)
       
        # i am not certain what the purpose of this code is
        with open(self.results_path, 'a+') as f:
            lines = f.readlines()
        if len(lines) > 1:
            return

        self.remove_old_output_files()
        # clear_folders()
        line = '{:>{}}'.format('Step', 5)
        o = Poscar('OSCAR')
        for i in range(N_CELL):
            for j in range(o.total_elements):
                line += '{:>3}'.format(o.element_names[j])
            line += '{:>10}'.format('E_DFT')
        for i in range(N_CELL):
            line += '{:>7}'.format('x_{}%'.format(i + 1))
        line += '{:>10}{:>5}\n'.format('E({:.0f}K)'.format(TEMPERATURE), 'P')
        print(line)
        exit()
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


    def check_for_stop_file():
        """Check for manual stop."""
        if os.path.exists(STOP_FILE):
            os.remove(STOP_FILE)
            save_log("Found stop file. Stopping...")
            exit()
if __name__ == "__main__":

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
