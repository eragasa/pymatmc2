# -*- coding: utf-8 -*-
"""This module contains classes to interact with GULP."""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017,2018,2019"
__license__ = "Simplified BSD License"
__version__ = "1.0"

from collections import OrderedDict
from mexm.exception import BadParameterException

from mexm.elements import ELEMENTS

from mexm.potential import MEXM_GLOBAL_FMT
from mexm.potential import MEXM_1BODY_FMT
from mexm.potential import MEXM_2BODY_FMT
from mexm.potential import MEXM_3BODY_FMT
from mexm.potential import get_symbol_pairs

class Potential(object):
    """base class for potential

    Args:
        symbols(list of str): a list of symbols to use for the potential
        potential_type(str): a description for the type of potential this is
        is_charge(bool): if set to true, this potential is a charge potential, default: False
    """

    potential_type = "potential"
    is_base_potential = True

    def __init__(self,
                 symbols,
                 is_charge=False):

        # define formatting strings

        # self.potential = None
        self.symbols = symbols
        self.is_charge = is_charge

        # these attributes will be initialized by _init_parameter_names
        self.symbol_pairs = None
        self.parameter_names = None
        self._initialize_parameter_names()

        # these attributes will be initialized by _init_parameter_names
        self.parameters = None
        self._initialize_parameters()

    @classmethod
    def get_parameter_names_global(cls, hybrid_format=True):
        parameter_names = []
        for parameter_name in cls.parameter_names_global:
            kwargs ={
                'potential_type':cls.potential_type,
                'parameter_name':parameter_name
            }
            if hybrid_format:
                parameter_names.append(MEXM_HYBRID_GLOBAL_FMT.format(**kwargs))
            else:
                parameter_names.append(MEXM_GLOBAL_FMT.format(**kwargs))

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
                if hybrid_format:
                    parameter_names.append(MEXM_HYBRID_1BODY_FMT.format(**kwargs))
                else:
                    parameter_names.append(MEXM_1BODY_FMT.format(**kwargs))

        return parameter_names

    @classmethod
    def get_parameter_names_2body(cls, symbols, hybrid_format=True):
        assert isinstance(symbols, list)
        assert isinstance(hybrid_format, bool)

        symbol_pairs = get_symbol_pairs(symbols)
        parameter_names = []
        for symbol_pair in symbol_pairs:
            for parameter_name in cls.parameter_names_2body:
                kwargs = {
                    'symbol1':symbol_pair[0],
                    'symbol2':symbol_pair[1],
                    'potential_type':cls.potential_type,
                    'parameter_name':parameter_name
                }
                if hybrid_format:
                    parameter_names.append(MEXM_HYBRID_2BODY_FMT.format(**kwargs))
                else:
                    parameter_names.append(MEXM_2BODY_FMT.format(**kwargs))
        return parameter_names

    def _initialize_parameter_names(self):
        raise NotImplementedError

    def _initialize_1body_parameter_names(self):
        for s in self.symbols:
            for p in type(self).one_body_parameters:
                parameter_name = MEXM_1BODY_FMT.format(s=s,p=p)
                self.parameter_names.append(parameter_name)


    def _initialize_2body_parameter_names(self):
        symbol_pairs = get_symbol_pairs(self.symbols)
        for sp in symbol_pairs:
            for p in type(self).two_body_parameters:
                parameter_name = MEXM_2BODY_FMT.format(
                    symbol1=sp[0],
                    symbol2=sp[1],
                    parameter_name=p
                )
                self.parameter_names.append(parameter_name)

    def _initialize_3body_parameter_names(self):
        raise NotImplementedError

    def _initialize_3body_parameter_names_bond_order(self):
        raise NotImplementedError

    def _initialize_3body_parameter_names_central_atom(self):
        raise NotImplementedError

    def _initialize_parameters(self):
        self.parameters = OrderedDict()
        for p in self.parameter_names:
            self.parameters[p] = None

    def evaluate(self, r, parameters, r_cut=False):
        """evaluate the potential

        Args:
            r(numpy.ndarray): a numpy array of interatomic distances
            parameters(OrderedDict): an dictionary of parameter values and keys
            r_cut(float,optional): the global cutoff for the potential
        """
        raise NotImplementedError

    def write_lammps_potential_file(self, path):
        """writes the lammps_potential file

        This method exists to write the lammps potential file to disk.  This
        method needs to be overriden to be implemented, in classes that
        inherit from this class.
        """

        raise NotImplementedError

    def lammps_potential_section_to_string(self):
        """generates the lammps string for the lammps potential sections

        Returns:
            str: the string value for the LAMMPS potential section that goes into the potential.mod file.
        """

        raise NotImplementedError

    def lammps_potential_section_set_masses_to_string(self, symbols=None):
        if symbols is None:
            symbols_ = self.symbols
        else:
            symbols_ = symbols

        str_out = ''
        for i, s in enumerate(symbols_):
            group_id = i+1
            amu = self._get_mass(s)
            str_out += "mass {} {}\n".format(group_id, amu)

        return str_out

    def lammps_potential_section_set_group_id_to_string(self, symbols=None):
        if symbols is None:
            symbols_ = self.symbols
        else:
            symbols_ = symbols

        str_out = ''
        for i, s in enumerate(symbols_):
            group_id = i+1
            symbol = s
            str_out += "group {} type {}\n".format(
                symbol,
                group_id
            )

        return str_out

    def lammps_potential_section_set_charges_to_string(self, symbols=None):
        if symbols is None:
            symbols_ = self.symbols
        else:
            symbols_ = symbols

        str_out = ""
        for i,s in enumerate(symbols_):
            charge_parameter_name = '{}_chrg'.format(s)
            charge = self.parameters[charge_parameter_name]
            str_out += "set group {} charge {}\n".format(s,charge)

        return str_out

    def write_gulp_potential_section(self):
        """writes the gulp potential for a GULP string"""
        raise NotImplementedError

    def gulp_potential_section_to_string(self):
        """generates the potential section for a GULP string"""
        raise NotImplementedError

    def _get_mass(self, symbol):
        return ELEMENTS[symbol].mass

    def _get_name(self, symbol):
        return ELEMENTS[symbol].name.lower()

if __name__ == "__main__":
    o = Potential(symbols=['Ni'])
