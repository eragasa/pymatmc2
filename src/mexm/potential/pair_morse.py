__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = 20171102
"""
this module impplements the morse potential
"""
from copy import deepcopy
from collections import OrderedDict
import numpy as np

from mexm.potential import PairPotential
from mexm.potential import get_symbol_pairs
from mexm.potential import (MEXM_GLOBAL_FMT,
                            MEXM_1BODY_FMT,
                            MEXM_2BODY_FMT,
                            MEXM_HYBRID_GLOBAL_FMT,
                            MEXM_HYBRID_1BODY_FMT,
                            MEXM_HYBRID_2BODY_FMT)
from mexm.potential import PairPotential


def function_morse_potential(r, D0, a, r0, cutoff):
    assert isinstance(r, np.ndarray) or isinstance(r, float)
    assert isinstance(D0, float)
    assert isinstance(a, float)
    assert isinstance(r0, float)
    assert isinstance(cutoff, float)

    if cutoff is None:
        return function_morse_potential_no_cutoff(r, D0, a, r0)
    else:
        V = function_morse_potential_no_cutoff(r, D0, a, r0)
        rcut = np.max(r[np.where(r < cutoff)])
        h = r[1] - r[0]
        V_rc = func_morse(rcut, D0, a, r0)
        V_rc_p1 = func_morse(rcut+_h, D0, a, r0)
        V_rc_m1 = func_morse(rcut-_h, D0, a, r0)
        dVdr_at_rc = (V_rc_p1-V_rc)/h

        # apply the cutoff
        V = V - V_rc - dVdr_at_rc * (r-rcut)
        V[np.where(r >= rcut)] = 0.0

def function_morse_potential_no_cutoff(r, D0, a, r0):
    return D0*(np.exp(-2*a*(r-r0))-2*np.exp(-a*(r-r0)))

class MorsePotential(PairPotential):
    potential_type = 'morse'
    is_charge = False

    parameter_names_global = ['cutoff']
    parameter_names_1body = [None]
    parameter_names_2body = ['D0', 'a', 'r0', 'cutoff']

    parameter_names_global_latex = {
        'cutoff':'r_{{c,g}}'
    }
    parameter_names_2body_latex = {
        'D0':'D_{{0,{symbol1}-{symbol2}}}',
        'a':'a_{{{symbol1}-{symbol2}}}',
        'r0':'r_{{0, {symbol1}-{symbol2}}}',
        'cutoff':'r_{{c, {symbol1}-{symbol2}}}'
    }

    """Implementation of the morse potential

    Args:
        symbols(list): a list of chemical symbols.

    Notes:
        parameter_names_global:
            cutoff (float6): global cutoff for Morse interactions, Angstroms
        parameter_names_1body:
            None
        parameter_names_2body:
            D0 (float): eV
            alpha (float): 1/Angs
            r0 (float): Angs
            cutoff (float): Angs
    """

    def __init__(self, symbols):
        super().__init__(symbols=symbols, is_charge=self.is_charge)


    # this method overrides the parents stub
    def _initialize_parameter_names(self):
        self.symbol_pairs = list(get_symbol_pairs(self.symbols))
        self.parameter_names = MorsePotential.get_parameter_names(
                symbols=self.symbols,
                hybrid_format=False
        )

    # this method overrides the parents stub
    def _initialize_parameters(self):
        self.parameters = OrderedDict()
        for v in self.parameter_names:
            self.parameters[v] = None

    # this method overrides the parent stub
    def evaluate(self, r, parameters):
        """evaluate the pair potential

        This method evaluates the MorsePotential
        Args:
            r(np.ndarray): A numpy array of interatomic distances which to evaluate.
            parameters(OrderedDict): A dictionary of parameters on which to evaluate the interatomic potential.
        Notes:

            >>> symbols=['Ni']
            >>> parameters=OrderedDict([
            >>>    ('NiNi_D0',0.001114),('NiNi_a',3.429506),('NiNi_r0',2.6813)
            >>>  ])
            >>> o = MorsePotential(symbols = s)
            >>> o.evaluate(r,testing_set['parameters'])

        """
        # <----------------------------check arguments are correct
        assert isinstance(r, np.ndarray)
        assert isinstance(parameters, dict)

        # <----------------------------copy a local of the parameters
        for k in self.parameters:
            self.parameters[k] = parameters[k]

        # <----------------------------evaluate the parameters now
        values = {}
        for s in self.symbol_pairs:
            # <------------------------extract the paramters for symbol pair
            parameters_ = {
                'D0':self.parameters['{}{}_D0'.format(s[0], s[1])],
                'a':self.parameters['{}{}_a'.format(s[0], s[1])],
                'r0':self.parameters['{}{}_r0'.format(s[0], s[1])],
                'cutoff':self.parameters['{}{}_cutoff'.format(s[0], s[1])]
            }
            # <------------------------embedded morse function
            pair_name = '{}{}'.format(s[0], s[1])
            values[pair_name] = function_morse_potential(r, **parameters_)

        return values

    # same as parent class
    def lammps_potential_section_to_string(self):
        """needs to be overridden"""
        raise NotImplementedError

    # same as parent class
    def gulp_potential_section_to_string(self):
        """needs to be overridden"""
        raise NotImplementedError

    # same as parent class
    def phonts_potential_section_to_string(self):
        """needs to be overridden"""
        raise NotImplementedError

    # same as parent class
    def write_lammps_potential_file(self):
        """needs to be overridden"""
        raise NotImplementedError

    # same as parent class
    def write_gulp_potential_section(self):
        """needs to be overridden"""
        raise NotImplementedError

    def references(self):
        reference_dict = {}
        reference_dict['LammpsMorse'] = "http://lammps.sandia.gov/doc/pair_morse.html"
        return reference_dict

if __name__ == "__main__":
    o = MorsePotential(symbols=['Ni'])
