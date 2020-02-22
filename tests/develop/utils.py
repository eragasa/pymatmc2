
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
