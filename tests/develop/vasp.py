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

