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


