__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyrirght (C) 2017"
__license__ = "Simplified BSD License"
__version__ = 20171102

import copy
import numpy as np
from collections import OrderedDict

from mexm.potential import PairPotential
from mexm.potential import get_symbol_pairs
from mexm.potential import (
        MEXM_GLOBAL_FMT, MEXM_1BODY_FMT, MEXM_2BODY_FMT,
        MEXM_HYBRID_GLOBAL_FMT, MEXM_HYBRID_1BODY_FMT, MEXM_HYBRID_2BODY_FMT)

class BuckinghamPotential(PairPotential):
    potential_type = 'buckingham'
    is_base_potential = False
    is_charge = True

    parameter_names_global = ['cutoff']
    parameter_names_1body = ['chrg']
    parameter_names_2body = ['A', 'rho', 'C', 'cutoff']

    """ Implementation of the Buckingham Potential

    This class provides an interface for the management of parameters
    to and from different molecular dynamics and lattice dynamics programs.

    Args:
        symbols(list): a list of chemicals

    Notes:
        parameter_names_global:
            cutoff (float): in Angstroms
        single_body_parameters:
            chrg (float): Coulombic charge in electron charge units
            cutoff (float): in Angstroms
        parameter_names_2body:
            A (float):
            rho (float):
            C (float):
    """
    def __init__(self,symbols):

        PairPotential.__init__(self,
                symbols=symbols,
                is_charge=self.is_charge)

    #override Potential method
    def _initialize_parameter_names(self):
        self.symbol_pairs = list(get_symbol_pairs(self.symbols))
        self.parameter_names = BuckinghamPotential.get_parameter_names(
                symbols=self.symbols,
                hybrid_format=False
        )

    #override Potential method
    def _initialize_parameters(self):
        self.parameters = OrderedDict()
        for v in self.parameter_names:
            self.parameters[v] = None

    def evaluate(self,r,parameters,r_cut=False):
        raise NotImplementedError

    def lammps_potential_section_to_string(self,parameters=None):
        """get the string for the lammps potential section

        Args:
            parameters(OrderedDict): a dictionary of parameter name keys with the associated values
            rcut(float): the global cutoff
        Returns:
            str: the string of the LAMMPS potential section
        """

        if parameters is not None:
            self.parameters = parameters

        str_out = self.lammps_potential_section_set_masses_to_string() + "\n"
        str_out += self.lammps_potential_section_set_group_id_to_string() +"\n"
        str_out += self.lammps_potential_section_set_charges_to_string() + "\n"

        global_cutoff = self.parameters['cutoff']
        str_out += 'pair_style buck/coul/long {}\n'.format(global_cutoff)

        for pair in get_symbol_pairs(self.symbols):
            s1 = pair[0]
            s2 = pair[1]
            str_A = '{}{}_A'.format(s1, s2)
            str_rho = '{}{}_rho'.format(s1, s2)
            str_C = '{}{}_C'.format(s1, s2)
            str_cutoff = '{}{}_cutoff'.format(s1, s2)

            group_id_s1 = self.symbols.index(s1) + 1
            group_id_s2 = self.symbols.index(s2) + 1
            A = float(self.parameters[str_A])
            rho = float(self.parameters[str_rho])
            C = float(self.parameters[str_C])
            cutoff = float(self.parameters[str_cutoff])

            args = [group_id_s1, group_id_s2, A, rho, C, cutoff]
            str_out += "pair_coeff {} {} {} {} {} {}\n".format(*args)

        return str_out

    # overrides the parents class
    def gulp_potential_section_to_string(self,parameters=None):
        """ get GULP potential to string

        The buckingham potential is a charged potential and so the charges
        associated with the potential are also part of the potential."
        """

        str_out = 'species\n'
        for s in self.symbols:
            chrg=self.parameters['{}_chrg'.format(s)]
            str_out += "{s} core {chrg}\n".format(s=s,chrg=chrg)
        str_out += "\n"
        str_out += 'buck\n'

        for symbols in self.symbol_pairs:
            s1 = symbols[0]
            s2 = symbols[1]
            sp = "{}{}".format(s1,s2)

            # get parameters
            A = parameters['{}_A'.format(sp)]
            rho = parameters['{}_rho'.format(sp)]
            C = parameters['{}_C'.format(sp)]
            cutoff = parameters['{}_cutoff'.format(sp)]
            str_out += "{s1} core {s2} core {A} {rho} {C} {cutoff}\n".format(
                s1=s1,
                s2=s2,
                A=A,
                rho=rho,
                C=C,
                cutoff=cutoff
            )

        return str_out

    # same as parent class
    def phonts_potential_section_to_string(self):
        raise NotImplementedError

    # same as parent class
    def write_lammps_potential_file(self):
        raise NotImplementedError

    # same as parent class
    def write_gulp_potential_section(self):
        raise NotImplementedError
