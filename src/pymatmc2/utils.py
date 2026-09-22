import os
import tarfile
import scipy.optimize 
import numpy as np
from pymatmc2 import Pymatmc2Log

#### THIS IS THE STUFF I HAVE WRITTEN
def validate_vasp_simulation_input_files_exist(path):
    assert os.path.isdir(path)
    vasp_required_files = ['INCAR', 'KPOINTS', 'POTCAR']
    has_all_files = all([
        os.path.isfile(os.path.join(path, k)) for k in vasp_required_files
    ])
    return has_all_files

#TODO: arr needs to be broken down in arguments a and b 
# once i figure out what it is
#TODO: this algorithm uses a NNLS algoirthm, it should be able to overdetermine
# the system and get statistical information out it!
def get_ratio(
    arr,
    concentration_matrix,
    concentration_total,
    n_cell
    ):
    """Calculate cell ratios from (n+1)*n input.

    f^1 c_1^1 + f^2 c_1^2 + f^3 c_1^3 = c_1
    f^1 c_2^2 + f^2 c_2^2 + f^3 c_2^3 = c_2
    f^1 c_3^2 + f^2 c_3^2 + f^3 c_3^3 = c_2

    C f = c
    C f - c = 0
    Ax = b
    then using an NNLS to solve the  problem.


    Arguments:
        arr (???): i'm assuming this is a numpy array???
        concentration_matrix = a matrix of (m x n) where m is the number of cells,
            and n is the number of chemical species.
        n_cell: I should be able to get this from the dimension of the array

    Returns:
        np.ndarray: ratio for something
    
    Notes:
        Uses the a non-negative least squares solver as described in:
        Lawson C., Hanson R.J., (1987) Solving Least Squares Problems, SIAM
    
    """
    n = n_cell
    concentration_matrix = arr[:, :n]   # is a matrix
    concentration_total = arr[:, n]    # 

    # uses non-negative least squares solver to solve for x
    # Solve argmin_x || Ax - b ||_2 for x>=0.
    cell_molar_fraction = scipy.optimize.nnls(
        concentration_matrix, 
        concentration_total
    )[0]

    # sanity check
    np.array_equal(
        np.dot(concentration_matrix, cell_molar_fraction),
        concentration_total
    )

    # these are establish bounds for determining if 
    ratio =  abs(cell_molar_fraction / np.sum(cell_molar_fraction))

    # TODO: remove if not necessary
    if False:
        if all(-1E-8 <= i <= 1+1E-8 for i in cell_molar_fraction) and all(abs(i) < 1E-6 for i in b-c):
            ratio =  abs(cell_molar_fraction / np.sum(cell_molar_fraction))
        else:
            ratio = np.zeros(n)

    return ratio

#### THIS IS THE STUFF I HAVE NOT WRITTEN


def archive_vasp_simulation(src_path, dst_path):
    assert os.path.isdir(src_path)
    assert not os.path.exists(dst_path)

    vasp_files_wanted = [
        'INCAR',
        'POSCAR',
        'CONTCAR',
    ]
    with tarfile.open(dst_path, 'a') as tarball:
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

def save_log(line):
    """Save a line to the LOG_FILE."""
    with open(LOG_FILE, 'a+') as f:
        f.write(line + '\n')



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

def clear_lock_file():
    lock_file_path = 'pymatmc2_lock'
    if os.path.exists(lock_file_path):
        os.remove(lock_file_path)


def stop_check():
    """Check for manual stop."""
    if os.path.exists(STOP_FILE):
        os.remove(STOP_FILE)
        save_log("Found stop file. Stopping...")
        exit()


def prepare(results_path):
    """Make sure all files exist."""
    vasp_required_files = ['INCAR', 'KPOINTS', 'POTCAR']
    for i in vasp_required_files:
        if not os.path.exists(i):
            save_log('{} not found.'.format(i))
            exit()
    with open(results_path, 'r') as f:
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
