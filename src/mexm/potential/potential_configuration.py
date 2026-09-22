from collections import OrderedDict

from mexm.potential import Potential

class PotentialConfiguration(object):

    def __init__(self):
        self._potential_dict = OrderedDict()
        self._potential_dict['potential_type'] = None
        self._potential_dict['symbols'] = None
        self._potential_dict['parameters'] = None

    @staticmethod
    def initialize_from_potential(potential):
        assert isinstance(potential, Potential)
        potential_configuration = PotentialConfiguration()
        potential_configuration.potential_type = type(potential).potential_type
        potential_configuration.symbols = potential.symbols
        potential_configuration.parameters = potential.parameters

        return potential_configuration

    @staticmethod
    def initialize_from_dict(potential):
        assert isinstance(potential, dict)
        potential_configuration = PotentialConfiguration()
        potential_configuration.potential_type = potential['potential_type']
        potential_configuration.symbols = potential['symbols']
        if 'parameters' in potential:
            potential_configuration.parameters = potential['parameters']

    @property
    def potential_type(self):
       return self._potential_dict['potential_type']

    @potential_type.setter
    def potential_type(self, potential_type):
        self._potential_dict['potential_type'] = potential_type

    @property
    def symbols(self):
        return self._potential_dict['symbols']

    @symbols.setter
    def symbols(self, symbols):
        self._potential_dict['symbols'] = symbols

    @property
    def parameters(self):
        return self._potential_dict['parameters']

    @parameters.setter
    def parameters(self, parameters):
        self._potential_dict['parameters'] = parameters

    def to_dict(self):
        return_dict = OrderedDict()
        for k,v in self._potential_dict.items():
            if v is None:
                pass
            else:
                return_dict[k] = v
        return return_dict
