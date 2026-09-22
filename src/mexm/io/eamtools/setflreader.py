from mexm.io.eamtools import BaseSetflFile

class SetflReader(BaseSetflFile):
    """Reader for setfl formatted EAM files.

    Args:
        path (str): Path to the setfl file.
    """

    def __init__(self, path=None):

        BaseSetflFile.__init__(self, path=path)
        if path is not None:
            self.read(path=self.path)

    def read(self, path=None):
        if path is not None:
            self.path = path
        path_ = self.path

        with open(path_, "r") as f:
            lines = f.readlines()

        self.comments_ = [lines[i].strip() for i in range(3)]

        symbol_info = lines[3].split()
        self.N_symbols = int(symbol_info[0])
        self.symbols = [str(k) for k in symbol_info[1:]]
        self.symbol_pairs = self.get_symbol_pairs(self.symbols)

        # cutoff info
        cutoff_info = lines[4].split()
        self.N_rho = int(cutoff_info[0])
        self.d_rho = float(cutoff_info[1])
        self.N_r = int(cutoff_info[2])
        self.d_r = float(cutoff_info[3])
        self.r_cutoff = float(cutoff_info[4])

        """Reads a setfl file."""
        body = {
            "lattice_info": {},
            "embedding_function": {},
            "density_function": {},
            "pair_function": {}
        }

        N_body_start = 5  # setfl body ALWAYS begins on line 6
        values = []  # separate lines by spaces for multi-colummn setfl format
        for line in lines[N_body_start:]:
            values += [s for s in line.split()]

        N_element_descriptor = 4
        start = 0  # now start is in reference to individual values; not lines
        for s in self.symbols:
            # read in element information
            body['lattice_info'] = {
                'atomic_number':int(values[start]),
                'atomic_weight':float(values[start+1]),
                'lattice_parameter':float(values[start+2]),
                'lattice_type':str(values[start+3])
            }

            # read in embedding information
            idx_embedding_start = start + N_element_descriptor
            idx_embedding_end = idx_embedding_start + self.N_rho
            body["embedding_function"][s] \
                = [float(k) for k in values[idx_embedding_start:idx_embedding_end]]
            assert len(body['embedding_function'][s]) == self.N_rho

            # read in density function values
            idx_density_start = idx_embedding_end + 1
            idx_density_end = idx_density_start + self.N_r
            body["density_function"][s] \
                = [float(k) for k in values[idx_density_start:idx_density_end]]
            assert len(body['density_function'][s]) == self.N_r
            start += self.N_rho + self.N_r + N_element_descriptor

        start -= N_element_descriptor
        for sp in self.symbol_pairs:
            # read in pair potential values
            idx_pair_start = start
            idx_pair_end = start + self.N_r
            values[idx_pair_start:idx_pair_end]
            body["pair_function"]["".join(sp)] \
                = [float(k) for k in values[idx_pair_start:idx_pair_end]]
            assert len(body['density_function'][s]) == self.N_r
            start += self.N_r

        self._body = body

    @property
    def N_rho(self):
        return self.Nrho

    @N_rho.setter
    def N_rho(self, Nrho):
        assert isinstance(Nrho, int)
        assert Nrho > 0

        self.Nrho = Nrho

    @property
    def N_r(self):
        return self.Nr

    @N_r.setter
    def N_r(self, Nr):
        assert isinstance(Nr, int)
        assert Nr > 0

        self.Nr = Nr


    def embedding_function(self, symbol):
        """Tabulated values of the embedding function.

        Args:
            symbol (str): Elemental symbol.

        Returns:
            list of float
        """
        return self._body["embedding_function"][symbol]

    def density_function(self, symbol):
        """Tabulated values of the density function.

        Args:
            symbol (str): Elemental symbol.

        Returns:
            list of float
        """
        return self._body["density_function"][symbol]

    def pair_function(self, symbol_pair):
        """Tabulated values of the pair function.

        Args:
            symbol_pair (str): PAir of elemental symbols
            - formatted as "{}{}".format(symbol1, symbol2)

        Returns:
            list of float
        """
        return self._body["pair_function"][symbol_pair]
