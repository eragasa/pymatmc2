from collections import OrderedDict
from copy import deepcopy

from mexm.potential import Potential
from mexm.potential import MEXM_HYBRID_GLOBAL_FMT
from mexm.potential import MEXM_HYBRID_1BODY_FMT

MEXM_CHARGE_TYPES = ['no_charge', 'static_charge', 'qeq_charge']
MEXM_CHARGE_SUMMATION_TYPES = ['ewald']

class Qeq(Potential):
    potential_type = 'qeq'
    is_base_potential = True
    parameter_names_global = ['Nevery', 'cutoff', 'tolerance', 'maxiter']
    parameter_names_1body = ['chi', 'eta', 'gamma', 'zeta', 'qcore']

    def __init__(self, symbols):
        assert isinstance(symbols, list)
        self.symbols = symbols

        self.parameter_names = Qeq.get_parameter_names(symbols)
        self.parameters_ = OrderedDict()

    @property
    def parameters(self):
        return self.parameters

    @parameters.setter
    def parameters(self, parameters):
        assert isinstance(parameters, OrderedDict)
        self.parameters_ = deepcopy(parameters)

    @classmethod
    def get_parameter_names(cls, symbols, hybrid_format=True):
        assert isinstance(hybrid_format, bool)

        parameter_names = []
        parameter_names += cls.get_parameter_names_global(hybrid_format)
        parameter_names += cls.get_parameter_names_1body(symbols, hybrid_format)

        return parameter_names

    @classmethod
    def get_parameter_names_global(cls, hybrid_format=True):
        assert isinstance(hybrid_format, bool)

        parameter_names = []
        for parameter_name in cls.parameter_names_global:
            kwargs = {
                'potential_type':cls.potential_type,
                'parameter_name':parameter_name
            }
            parameter_names.append(MEXM_HYBRID_GLOBAL_FMT.format(**kwargs))

        return parameter_names

    @classmethod
    def get_parameter_names_1body(cls, symbols, hybrid_format=True):
        assert isinstance(symbols, list)
        assert isinstance(hybrid_format, bool)

        parameter_names = []
        for symbol in symbols:
            for parameter_name in cls.parameter_names_1body:
                kwargs = {
                    'symbol':symbol,
                    'potential_type':cls.potential_type,
                    'parameter_name':parameter_name
                }
                parameter_names.append(MEXM_HYBRID_1BODY_FMT.format(**kwargs))

        return parameter_names

    def lammps_fix_qeq_to_string(self, lammps_fix_id, lammps_group_id, qfile_path):
        qeq_style = 'qeq'

        qeq_fix_fmt = "fix {lammps_fix_id} {lammps_group_id} {qeq_style} {Nevery} {cutoff} {tolerance} {maxiter} {qfile}"

        float_fmt = '{:10.6e}'
        kwargs = {
            'lammps_fix_id':str(lammps_fix_id),
            'lamamps_group_id':str(lammps_group_id),
            'qeq_style':Qeq.potential_type,
            'Nevery':int(self.parameters['Nevery']),
            'cutoff':float_fmt.format(self.parameters['cutoff']),
            'tolerance':float_fmt.format(self.parameters['tolerance']),
            'maxiter':float_fmt.format(self.parameters['maxiter']),
            'qfile':qfile_path
        }

        return qeq_fix_fmt.format(**kwargs)
