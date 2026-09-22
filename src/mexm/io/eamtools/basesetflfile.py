import numpy as np

class BaseSetflFile():
    def __init__(self, path=None):
        assert any([isinstance(path, str), path is None])

        self.path = path
        self.comments = []
        self.symbols = []
        self.symbol_pairs = []

        self.Nrho = None
        self.drho = None
        self.maxrho = None
        self.rho = []

        self.Nr = None
        self.dr = None
        self.rmax = None
        self.r = {}

        self.pair = {}
        self.embedding = {}
        self.density = {}
        self.lattice_info = {}

    @staticmethod
    def get_interatomic_distance_vector(rmax, Nr):
        assert isinstance(rmax, float)
        assert isinstance(Nr, float)
        assert Nr % 5 == 0

        r = rmax * np.linspace(1, Nr, Nr)/Nr
        dr = rmax / Nr

        return r, dr

    @staticmethod
    def get_electron_density_vector(rho_max, Nrho):
        assert isinstance(rho_max, float)
        assert isintance(Nr, float)
        assert Nr % 5 == 0

        rho = rho_max * np.linspace(1, Nrho, Nrho)/Nrho
        drho = rho_max/ Nrho
        return rho, drho

    @staticmethod
    def get_symbol_pairs(symbols):
        """determine symbol pairs

        given a list of symbols gives a list of symbol pairs in the appropriate order expected within the pypospack package.

        Args:
            symbols(list of str, str): a list of symbols
        Returns:
            (list of list):a list of symbol pairs
        Exception:
            (TypeError) if symbols is neither a list nor a string

        """

        if isinstance(symbols, str):
            symbols_ = [symbols]
        elif isinstance(symbols, list):
            symbols_ = symbols
        else:
            msg = (
                "symbols argument must either be a list of chemical symbols",
                "or a chemical symbol")
            raise TypeError(msg)

        symbol_pairs = []
        for i1,s1 in enumerate(symbols_):
            for i2,s2 in enumerate(symbols_):
                if i1 <= i2:
                    symbol_pairs.append([s1,s2])

        return symbol_pairs
