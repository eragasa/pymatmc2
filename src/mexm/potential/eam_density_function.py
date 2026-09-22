from mexm.potential import Potential

class EamDensityFunction(Potential):
    potential_type = 'eam_density_base'
    is_base_potential = True
    def __init__(self,
            symbols,
            potential_type='eamdens'):

        super().__init__(symbols=symbols,
                         is_charge=False)

        self.density_evaluations = None

    @classmethod
    def get_parameter_names(cls, symbols, hybrid_format=False):
        assert isinstance(hybrid_format, bool)

        parameter_names = []
        parameter_names += cls.get_parameter_names_1body(symbols, hybrid_format)

        return parameter_names

    def _initialize_parameter_names(self,symbols=None):
        if symbols is None:
            symbols = self.symbols
        self.parameter_names = self.get_parameter_names(symbols, hybrid_format=False)

# from pypospack.potentials.eam_dens_exp import ExponentialDensityFunction
