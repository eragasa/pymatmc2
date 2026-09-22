class Outcar(object):
    """ object to manage IO to VASP through the OUTCAR file

    Many of the attributes in this class are initialized to None at
    instanciation.  They are processed when the OUTCAR file is read.

    Args:
        filename (str): filename of the OUTCAR file

    Attributes:
        total_energy (float): total energy of the structure in eV.
        elastic_tensor (nd.array): the components of the elastic tensor.  The
            components of the array start at :math:`0`.  So the :math:`c_{11}`
            component would be accessed by elastic_tensor[0,0]
        phonon_eig_val (nd.array): the eigenvalues of the phonons determined
            by lattice dynamics.
        phonon_eig_vec (nd.array): the eigenvectors of the phonons determined
            by lattice dynamics.
    """
    def __init__(self, filename="OUTCAR"):
        self.filename = filename
        self.total_energy = None
        self.encut = None
        self.entropy = None
        self.elastic_tensor = None
        self.phonon_eig_val = None # phonon eigenvalues
        self.phonon_eig_vec = None # phonon eigenvectors
        self.magnetic_moments = None

    def read(self, path=None):
        """ read a VASP outcar file and parse for information

        Args:
            filename (str): filename
            encut (float): energy cutoff for this simulation
            total_energy (float): total energy for this simulation
        """

        if path is not None:
            self.path = path
        with open(self.path) as f:
            for line in f:
                # check if free energy line
                if "TOTEN" in line:
                    try:
                        E = line.strip().split('=')[1].strip().split(' ')[0]
                        E = float(E)
                        self.total_energy = E
                    except ValueError:
                        if type(self.total_energy) is not float:
                            pass
                    except IndexError:
                        if type(self.total_energy) is not float:
                            pass
                elif "ENCUT" in line:
                    E = line.strip().split('=')[1].strip().split(' ')[0]
                    E = float(E)
                    self.encut = E
                elif "EENTRO" in line:
                    try:
                        E = line.strip().split('=')[1].strip()
                        E = float(E)
                        self.entropy = E
                    except ValueError as e:
                        if type(self.entropy) is not float:
                            pass
                    except IndexError as e:
                        if type(self.entropy) is not float:
                            pass
                else:
                    pass

    def __string__():
        print("total_energy[eV] = {}".format(self._ecoh_per_structure))

    def get_phonons(self):
        pass

    def get_time_of_calculation(self):
        pass

    def get_number_of_atoms(self):
        pass

    def get_energy(self, fname="OUTCAR"):
        self.ecoh_per_structure = None
        self.ecoh_per_atom = None

    def get_magnetic_moments(self, filename=None):
        if filename is None:
            filename_ = self.filename
        else:
            self.filename = filename
            filename_ = filename

        with open(filename_) as f:
            line = f.readline()
            while line:
                if "magnetization (x)" in line:
                    self.magnetic_moments = []
                    line = f.readline()
                    line = f.readline()
                    line = line.replace('# of ion','atom_id')
                    line = line.replace('tot','total')
                    labels = line.strip().split()
                    while '---' not in line:
                        line = f.readline()
                    line = f.readline()
                    while '---' not in line:
                        values = [float(k) for k in line.strip().split()]
                        self.magnetic_moments.append(values[1:])
                        line = f.readline()
                line = f.readline()
        return self.magnetic_moments
