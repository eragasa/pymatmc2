# coding: utf-8
# Copyright (c) Eugene J. Ragasa
# Distributed under the terms of the MIT License

""" collection of INCAR tags for VASP

This module contains base classes for INCAR tags, not every
INCAR tag has been implemented.  In addition, not every INCAR
tag has been implemented for output in the mexm.io.vasp.Incar
class.  If there is an issue with tag support, please contact me
and it can be implemented fairly quickly

Base Classes:
    IncarBaseTag
    IncarBaseTags (deprecated, use IncarBaseTag instead)
    IncarBaseEnumeratedTag
    IncarBaseIntegerTag
    IncarBaseStringTag
    IncarBaseFloatTag
"""

__author__ = "Eugene J. Ragasa"
__email__ = "ragasa.2@osu.edu"
__copyright__ = "Copyright 2020, Eugene J. Ragasa"
__maintainer__ = "Eugene J. Ragasa"
__date__ = "2020/02/22"

from collections import OrderedDict
from mexm.io.vasp.errors import VaspIncarError

class IncarBaseTag():
    """ abstract base class for all IncarTag """

    # static variables
    tag_name = "BASETAG"
    tag_dictionary = OrderedDict()

    @classmethod
    def is_valid_option(cls, option):
        if option in cls.tag_dictionary:
            return True
        else:
            return False

    @classmethod
    def get_comment(cls, option):
        try:
            return cls.tag_dictionary[option]
        except KeyError:
            msg = "{}unknown option, {}".format(
                str(cls),
                option
            )
            raise VaspIncarError(msg)

    @classmethod
    def cast_option(cls, option):
        raise NotImplementedError

class IncarBaseTags(IncarBaseTag): pass

class IncarBaseEnumeratedTag(IncarBaseTag):
    """ abstract base class for INCAR tags with enumerated values """
    tag_name = 'ENUMERATEDTAG'
    tag_dictionary = {}

    @classmethod
    def is_valid_option(cls, option):
        if option in cls.tag_dictionary:
            return True
        else:
            return False

class IncarBaseIntegerTag(IncarBaseTag):
    tag_name = "INTEGERTAG"
    comment = ""

    @classmethod
    def is_valid_option(cls, option: int):
        """
        Arguments:
            option (integer): value of the tag

        Returns
            bool: whether or not the tag is valid
        """
        return isinstance(option, int)

    @classmethod
    def get_comment(cls, option):
        return cls.comment

class IncarBaseStringTag(IncarBaseTag):
    tag_name = "STRINGTAG"
    comment = ""

    @classmethod
    def is_valid_option(cls, option):
        return isinstance(option, str)

    @classmethod
    def get_comment(cls, option):
        raise NotImplementedError

class IncarBaseFloatTag(IncarBaseTag):
    tag_name = 'fake tag'
    comment = 'fake comment'

    @classmethod
    def is_valid_option(cls, option):
        if any([
                isinstance(option,int),
                isinstance(option,float)
                ]):
            return True
        else:
            return False

    @classmethod
    def get_comment(cls, option):
        return cls.comment

class AlgoTag(IncarBaseEnumeratedTag):
    tag_name = 'ALGO'
    tag_dictionary = OrderedDict([
        ('Normal', 'blocked Davidson iteration scheme'),
        ('VeryFast', 'RMM-DIIS'),
        ('Fast', 'blocked Davidson, followed by RMM_DIIS')
    ])

    @classmethod
    def convert(cls, option):
        return str(option)

class AmixTag(IncarBaseFloatTag):
    tag_name = 'AMIX'
    comment = 'linear mixing parameter mixing scheme'

class AmixmagTag(IncarBaseFloatTag):
    tag_name = "AMIX_MAG"
    comment = "linear mixing parameter for the magnetization density"

class BmixTag(IncarBaseFloatTag):
    tag_name = 'BMIX'
    comment = 'cutoff wave vector for Kerker mixing scheme'

class BmixmagTag(IncarBaseFloatTag):
    tag_name = 'BMIX_MAG'
    comment = 'cutoff wave vector for magnetization density'

class LchargTag(IncarBaseEnumeratedTag):
    tag_name = 'LCHARG'
    tag_dictionary = OrderedDict([
        ('.TRUE.', 'write CHGCR, write CHG'),
        ('.FALSE.', 'no CHGCAR, no CHG')
    ])

    @classmethod
    def convert(cls, option):
        if option in [True, 'T', '.TRUE.']:
            return '.TRUE.'
        elif option in [False, 'F', '.FALSE.']:
            return '.FALSE.'
        else:
            return str(option)


class LrealTags(IncarBaseEnumeratedTag):
    tag_name = 'LREAL'
    tag_dictionary = OrderedDict([
        ('.FALSE.', 'projection done in reciprocal space'),
        ('On', 'method of King-Smith, et al. Phys. Rev B 44, 13063 (1991).'),
        ('Auto', 'unpublished method of G. Kresse')
    ])

    @classmethod
    def convert(cls, option):
        return str(option)

class LorbitTags(IncarBaseEnumeratedTag):
    tag_name = 'LORBIT'
    tag_dictionary = OrderedDict([
        (0, 'DOSCAR and PROCAR'),
        (1, 'DOSCAR and lm-decomposted PROCAR'),
        (2, 'DOSCAR, lm_decomposed PROCAR, phase factors'),
        (5, 'DOSCAR, PROCAR'),
        (10, 'DOSCAR, PROCAR'),
        (11, 'DOSCAR, lm-decomposed PROCAR'),
        (12, 'DOSCAR, lm-decomposed PROCAR, phase factors'),
    ])

    @classmethod
    def convert(cls, option) -> int:
        return int(option)

class LplaneTag(IncarBaseEnumeratedTag):
    tag_name = 'LORBIT'
    tag_dictionary = OrderedDict([
        ('.TRUE.', 'use real space projectors, small number of compute nodes'),
        ('.FALSE.', 'no real space projectors, large number of compute nodes')
    ])

    @classmethod
    def convert(cls, option):
        if option in [True, 'T', '.TRUE.']:
            return '.TRUE.'
        elif option in [False, 'F', '.FALSE.']:
            return '.FALSE.'
        else:
            return str(option)

class LwaveTag(IncarBaseEnumeratedTag):
    tag_name = 'LWAVE'
    tag_dictionary = OrderedDict([
        ('.TRUE.', 'write WAVECAR'),
        ('.FALSE.', 'do not write WAVECAR')
    ])

    @classmethod
    def convert(cls, option):
        if option in [True, 'T', '.TRUE.']:
            return '.TRUE.'
        elif option in [False, 'F', '.FALSE.']:
            return '.FALSE.'
        else:
            return str(option)

class LvtotTag(IncarBaseEnumeratedTag):
    tag_name = 'LVTOT'
    tag_dictionary = OrderedDict([
        ('.TRUE.', 'write LOCPOT'),
        ('.FALSE.', 'no LOCPOT')
    ])

    @classmethod
    def convert(cls, option):
        if option in [True, 'T', '.TRUE.']:
            return '.TRUE.'
        elif option in [False, 'F', '.FALSE.']:
            return '.FALSE.'
        else:
            return str(option)

class EdiffTag(IncarBaseFloatTag):
    tag_name = 'EDIFF'
    comment = 'convergence condition for SC-loop in eV'

class EdiffgTag(IncarBaseFloatTag):
    tag_name = 'EDIFFG'
    comment_force_relaxation = 'force convergence requirements in ev A'
    comment_energy_relaxation = 'energy convergence in eV'

    @classmethod
    def is_valid_option(cls, option):
        try:
            option_ = float(option)
        except ValueError:
            return False

        if option_ == 0:
            return False
        else:
            return True
        
    @classmethod
    def get_comment(cls,option):
        try:
            option_ = float(option)
        except ValueError:
            msg = 'value must be castable into a float'
            raise TypeError(msg)

        if option == 0:
            msg = 'cannot select zero, must select a convergence value'
            raise ValueError(msg)
        elif option_ < 0:
            return cls.comment_force_relaxation
        else:
            return cls.comment_energy_relaxation

class EncutTag(IncarBaseFloatTag):
    tag_name = 'ENCUT'
    comment = 'Cut-off energy for plane wave basis set in eV'

    @classmethod
    def is_valid_option(cls, option):
        try:
            option_ = float(option)
        except ValueError:
            return False
        
        if option > 0:
            return True
        else:
            return False

class IbrionTag(IncarBaseEnumeratedTag):
    tag_name = 'IBRION'
    tag_dictionary = OrderedDict([
        (-1, 'no update'),
        (0, 'molecular dynamics'),
        (1, 'ionic relaxation by RMM-DIIS'),
        (2, 'ionic relaxation by CG'),
        (3, 'ionic relaxation by damped MD'),
        (5, 'phonons, by frozen ion, without symmetry'),
        (6, 'phonons, by frozen ion, with symmetry'),
        (7, 'phonons, by perturbation theory, no symmetry'),
        (8, 'phonons, by perturbtion theory, with symmetry')
    ])

    @classmethod
    def convert(cls, option):
        return int(option)

class IniwaveTag(IncarBaseIntegerTag):
    tag_name = 'INIWAV'
    tag_dictionary = OrderedDict([
        (0, 'jellium wave functions'),
        (1, 'random wave functions')
    ])

    @classmethod
    def is_valid_option(cls, option):
        return option in cls.tag_dictionary

    @classmethod
    def convert(cls, option) -> int:
        return int(option)

class IsifTag(IncarBaseEnumeratedTag):
    tag_name = 'ISIF'
    tag_dictionary = OrderedDict([
        (0, 'relaxation, ions=T, cellshape=F, cellvolume=F'),
        (1, 'relaxation, ions=T, cellshape=F, cellvolume=F'),
        (2, 'relaxation, ions=T, cellshape=F, cellvolume=F'),
        (3, 'relaxation, ions=T, cellshape=T, cellvolume=T'),
        (4, 'relaxation, ions=T, cellshape=T, cellvolume=F'),
        (5, 'relaxation, ions=F, cellshape=T, cellvolume=F'),
        (6, 'relaxation, ions=F, cellshape=T, cellvolume=T'),
        (7, 'relaxation, ions=F, cellshape=F, cellvolume=T')
    ])

    @classmethod
    def convert(cls, option):
        return int(option)


class IstartTag(IncarBaseEnumeratedTag):
    tag_name = 'ISTART'
    tag_dictionary = OrderedDict([
        (0,'begin from scratch'),
        (1,'continuation job, constant energy cutoff'),
        (2,'continuation job, constant basis set'),
    ])

    @classmethod
    def is_valid_option(cls, option):
        try:
            option_ = int(option)
        except:
            option_ = option

        if option_ in cls.tag_dictionary:
            return True
        else:
            return False

    @classmethod
    def convert(cls, option):
        return int(option)
    
class IsymTag(IncarBaseEnumeratedTag):

    tag_dictionary = OrderedDict([
        (-1,'symmetry off'),
        (0,'symmetry_on'),
        (1,'symmetry_on'),
        (2,'symmetry_on, efficient symmetrization'),
        (3,'symmetry_on, only forces and stress tensor')
    ])

    @classmethod
    def convert(cls, option) -> int:
        return int(option)

class IspinTag(IncarBaseEnumeratedTag):
    tag_name = 'ISPIN'
    tag_dictionary = OrderedDict([
        (1, 'non-spin polarized calculations'),
        (2, 'spin polarized calculations')
    ])

    @classmethod
    def convert(cls, option):
        return int(option)


class IchargTag(IncarBaseEnumeratedTag):
    tag_name = 'ICHARG'
    tag_dictionary = OrderedDict([
        (0,'Calculate charge density from initial wave functions.'),
        (1,'Read the charge density from file CHGCAR'),
        (2,'Take superposition of atomic charge densities')
    ])

    @classmethod
    def is_valid_option(cls, option):
        try:
            option_ = int(option)
        except:
            option_ = option

        if option_ in cls.tag_dictionary:
            return True
        else:
            return False

    @classmethod
    def convert(cls, option):
        return int(option)

class IsmearTag(IncarBaseEnumeratedTag):
    tag_name = 'ISMEAR'
    tag_dictionary = OrderedDict([
        (-5,'tetrahedron method with Blochl corrections'),
        (-4,'tetrahedron method'),
        (0,'method of Gaussian smearing'),
        (1,'method of Methfessel-Paxton order 1'),
        (2,'method of Methfessel-Paxton order 2')
    ])

    @classmethod
    def convert(cls, option):
        return int(option)


class KparTag(IncarBaseIntegerTag):
    """ INCAR tag KPAR

    To obtain high efficiency on massively parallel systems or modern multi-
    core machines, it is strongly recommended to use all at the same time. 
    Most algorithms work with any data distribution.
    """
    tag_name = 'KPAR',
    tag_default = 1,
    comment = "number of kpoints to paralellize"

    @classmethod
    def is_valid_option(cls, option: int):
        """ determine if the prposed option is an appropriate tag

        Arguments:
            option (int): this method will see if the option can be cast
               into an int
        """
        try:
            option_ = int(option)
        except ValueError:
            return False

        if option_ > 0:
            return True
        else:
            return False

    @classmethod
    def convert(self, option) -> int:
        try:
            option_ = int(option)
        except ValueError:
            msg = "KPAR tag must be castable into an integer"

        if option_ > 0:
            return option_
        else:
            msg = "KPAR tag must be greater than zero"
            return option_
 
class NcoreTag(IncarBaseIntegerTag):
    tag_name = 'NCORE',
    tag_default = 1,
    comment = "compute cores per orbital"

    @classmethod
    def is_valid_option(cls, option:int):
        """ determine if the proposed option in an approptiate tag

        Arguments:
            option (int): this method will if the option can be cast to an int
        """

        try:
            option_ = int(option)
        except ValueError:
            return False

        if option_ > 0:
            return True
        else:
            return False
    
    @classmethod
    def convert(self, option) -> int:
        try:
            option_ = int(option)
        except ValueError:
            msg = "NCORE tag must be castable into an integer"

        if option_ > 0:
            return option_
        else:
            msg = "NCORE tag must be greater than zero"
            return option_

class NelmTag(IncarBaseIntegerTag):
    tag_name = 'NELM'
    comment = 'maximum number of electronic SCF steps'

    @classmethod
    def is_valid_option(cls, option):
        try:
            option_ = float(option)
        except ValueError:
            return False

        if option_ > 0:
            return True
        else:
            return False

class NelmdlTag(IncarBaseIntegerTag):
    tag_name = 'NELMDL'
    comment = "the number of non-selfconsistent steps at the beginning."

    @classmethod
    def is_valid_option(cls, option):
        try:
            option_ = int(option)
            return True
        except ValueError:
            return False


class NelminTag(IncarBaseIntegerTag):
    tag_name = "NELMIN"
    comment = 'minimum number of electronic SCF steps'

    @classmethod
    def is_valid_option(cls, option):
        try:
            option_ = float(option)
        except ValueError:
            return False

        if option_ > 0:
            return True
        else:
            return False

class PrecTag(IncarBaseEnumeratedTag):
    tag_name = 'PREC'
    tag_dictionary = OrderedDict([
        ('Accurate', 'avoid wrap around errors'),
        ('High', 'avoid wrap around errors')
    ])

    @classmethod
    def convert(cls, option):
        return str(option)

class NparTag(IncarBaseIntegerTag):
    tag_name = 'NPAR'
    comment = 'number of bands to be treated in parallel'
    
    @classmethod
    def is_valid_option(cls, option):
        if not isinstance(option, int):
            return False
        else:
            if option <= 0:
                return False
            else:
                return True

class NsimTag(IncarBaseIntegerTag):
    tag_name = 'NSIM'
    comment = 'number of bands to optimize simultaneously'
    
    @classmethod
    def is_valid_option(cls, option):
        if not isinstance(option, int):
            return False
        else:
            if option <= 0:
                return False
            else:
                return True
                
class PotimTag(IncarBaseFloatTag):
    tag_name = 'POTIM'
    comment = 'scaling factor in relaxation'

class NswTag(IncarBaseIntegerTag):
    tag_name = 'NSW'
    comment = 'maximum number of ionic relaxation steps'

class SigmaTag(IncarBaseFloatTag):
    tag_name = 'SIGMA'
    comment = 'width of the smearing in eV.'

class SystemTag(IncarBaseStringTag):
    tag_name = 'SYSTEM'
    comment = ""
    @classmethod
    def is_valid_option(cls, option: str) -> bool:
        if isinstance(option, str):
            return True
        else:
            return False

class SymprecTag(IncarBaseFloatTag):
    tag_name = 'SYMPREC'
    comment = 'determines how accurate positions must be'


class IncarComments():
    tag_dictionary = {
        'ISTART':IstartTag,
        'ISYM':IsymTag,
        'SYMPREC':SymprecTag
    }

    @staticmethod
    def get_comment(tag_name, tag_option):
        return IncarComments.tag_dictionary[tag_name].get_comment(tag_option)

