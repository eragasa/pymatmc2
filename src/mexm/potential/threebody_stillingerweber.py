# -*- coding: utf-8 -*-
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017,2018"
__license__ = "Simplified BSD License"
__version__ = "1.0"

from collections import OrderedDict
from mexm.potential import MEXM_2BODY_FMT
from mexm.potential import MEXM_3BODY_FMT
from mexm.potential import ThreeBodyPotential

class StillingerWeberPotential(ThreeBodyPotential):
    two_body_parameters = ['A','B','p','q','epsilon','sigma', 'a']
    three_body_parameters = ['lambda', 'costheta0', 'tol']
    potential_type = 'stillingerweber'
    is_base_potential = False
    is_charge = False

    """Implementation of the Stillinger-Weber potential
    Args:
        symbols: list of string
    Attributes:
        symbols
        potential_type
        is_charge
    References:
        http://lammps.sandia.gov/doc/pair_sw.html
    """

    def __init__(self,symbols):
        _is_charge = False

        ThreeBodyPotential.__init__(self,
                symbols=symbols,
                is_charge=_is_charge)

        self.lmps_parameter_filename = "lmps_parameter_filename"

    def _initialize_parameter_names(self):
        self.parameter_names = []
        print('initialize_2body_parameter_names')
        self._initialize_2body_parameter_names()
        print(self.parameter_names)
        self._initialize_3body_parameter_names()
        print(self.parameter_names)

    def _initialize_3body_parameter_names(self):
        self.three_body_triplets = []
        for i, s1 in enumerate(self.symbols):
            for j, s2 in enumerate(self.symbols):
                for k, s3 in enumerate(self.symbols):
                    if j <= k:
                        triplet = [s1, s2, s3]
                        self.three_body_triplets.append(triplet)

        for triplet in self.three_body_triplets:
            for p in StillingerWeberPotential.three_body_parameters:
                parameter_name = MEXM_3BODY_FMT.format(
                    symbol1=triplet[0],
                    symbol2=triplet[1],
                    symbol3=triplet[2],
                    parameter_name=p
                )
                self.parameter_names.append(parameter_name)


    def lammps_potential_section_to_string(self,parameters=None):

        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]

        fname_params= self.lmps_parameter_filename

        str_out = ''
        for i,s in enumerate(self.symbols):
            str_out += "mass {} {}\n".format(i+1,self._get_mass(s))
        str_out += "\n"

        for i,s in enumerate(self.symbols):
            str_out += "group {} type {}\n".format(s,i+1)
        str_out += "\n"

        str_out += "pair_style sw\n"
        for s in self.symbols:
            str_out += "pair_coeff * * {} {}\n".format(fname_params,s)
        str_out += "\n"

        #<--------- neighbor lists moved to pypospack.task.lammps.LammpsTask
        #str_out += "neighbor 1.0 bin\n"
        #str_out += "neigh_modify every 1 delay 0 check yes\n"

        return str_out

    def write_lammps_parameter_file(self,dst_dir,dst_filename):
        assert isinstance(dst_dir,str)
        assert isinstance(dst_filename,str)

        _strout = self.lammps_paramfile_file_to_string()
        with open(os.path.join(dst_dir,dst_filename)) as f:
            f.write(_strout)

    def lammps_parameter_file_to_string(self,parameters=None):
        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]

        str_out = ''
        for i, s1 in enumerate(self.symbols):
            for j, s2 in enumerate(self.symbols):
                for k, s3 in enumerate(self.symbols):
                    s = '{}{}{}'.format(s1,s2,s3)
                    str_out += '{} {} {}'.format(s1,s2,s3)
                    str_out += ' ' + str(self.parameters['{}_epsilon'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_sigma'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_a'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_lambda'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_gamma'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_costheta0'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_A'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_B'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_p'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_q'.format(s)])
                    str_out += ' ' + str(self.parameters['{}_tol'.format(s)])
                    str_out += '\n'

        return str_out
