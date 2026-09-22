# -*- coding: utf-8 -*-
import os, shutil, pathlib
import re
import copy
#import pyflamestk.base as base
#import pypospack.crystal as crystal
import numpy as np

# *****************************************************************************
# ****    ERROR EXCEPTION HANDLING CLASSES                                 ****
# *****************************************************************************
from mexm.io.vasp.incar import VaspIncarError
from mexm.io.vasp.poscar import VaspPoscarError
from mexm.io.vasp.potcar import VaspPotcarError

# *****************************************************************************
# ****    SOME CONVENIENCE FUNCTIONS                                       ****
# *****************************************************************************

def read_poscar_file(filename):
    poscar = Poscar()
    poscar.read(path=filename)
    return poscar

def write_poscar_file(self,obj_poscar,filename='POSCAR'):
    obj_poscar.write(filename)

def read_incar_file(filename):
    incar = Incar()
    incar.read(path=filename)
    return incar

def make_super_cell(obj, scp):
    sc = base.make_super_cell(copy.deepcopy(obj), list(scp))
    return copy.deepcopy(Poscar(sc))

# input files
from mexm.io.vasp.incar import Incar
from mexm.io.vasp.poscar import Poscar
from mexm.io.vasp.kpoints import Kpoints
from mexm.io.vasp.potcar import Potcar
class VaspInputFiles():
    repository = {
        'INCAR':Incar,
        'POSCAR':Poscar,
        'KPOINTS':Kpoints,
        'POTCAR':Potcar,
    }



class Bsefatband():
    default_filename = 'BSEFATBAND'
    desc = "	BSE eigenvalues used for \"fatband\" plots."

class Chg():
    default_filename = 'CHG'
    description = (
        "Contains charge density, lattice vectors, and atomic coodinates "
        "Should be used for visualization"
    )

class Chgcar():
    default_filename = 'CHGCAR'
    description = (
        "Same as CHG but it contains also one-center occupancies. Should "
        "be used to restart VASP from existing charge density."
    )
from mexm.io.vasp.contcar import Contcar

class Doscar():
    default_filename = 'DOSCAR'
    description = (
        "Contain density of stats and integrated density of states"
    )

class Eigenval():
    default_filename = 'EIGENVAL'
    description = (
        "Contains Kohn-Sham eigenvalues for each kpoint after the end of the "
        "calculation"
    )

class Elfcar():
    default_filename = 'ELFCAR'
    description = (
        'Contains electron localization function.'
    )

class Ibzkpt():
    default_filename = 'IBZKPT'
    description = (
        "Contains k-point coordinates and weights."
        
    )

class Locpot():
    default_filename = 'LOCPOT'
    description = (
        "Contains total local potential in eV."
    )
from mexm.io.vasp.oszicar import Oszicar
from mexm.io.vasp.outcar import Outcar

class Parchg():
    default_filename = 'PARCHG'
    description = (
        'Contains partial charge densities.'
    )
class Procar():

    default_filename = 'PROCAR'
    description = 'Contains spd and site-projected wave function character.'

class Report():
    default_filename = 'REPORT'
    description = 'Contains output of various molecular dynamics calculations'

class TmpCar():
    default_filename = 'TMPCAR'
    description = 'Contains wavefunction and ionic positions of previous ionic step.'

class Vasprunxml():
    default_filename = 'vasprun.xml'
    description = 'main output file in xml format'

class BseDiagonalScreenedExchange():
    default_filename = 'WXXXX.tmp'
    description = (
        "Contains diagonal elements of screened exchange in BSE calculations."
    )

class Wavecar():
    default_filename = 'WAVECAR'
    description = (
        "Binary file containing information such as wave function coefficients," "eigenvalues, Fermi weights, etc."
    )
class Waveder():
    default_filename = 'WAVEDER'
    description = (
        'Contains derivative of wave functions with respect to k point.'
    )
class BseScreenedExchange():
    default_filename = 'WFULLxxxx.tmp'
    description = 'Store full screened exchange in BSE calculations.'
class Xdatcar():
    default_filename = 'XDATCAR'
    description = 'Contains ionic configuration for each output step of molecular dynamics simulations.'

# output files
class VaspOutputFiles():
    repository = {
        'BSEFATBAND':Bsefatband,
        'CHG':Chg,
        'CHGCAR':Chgcar,
        'CONTCAR':Contcar,
        'DOSCAR':Doscar,
        'EIGENVAL':Eigenval,
        'ELFCAR':Elfcar,
        'IBZKPT':Ibzkpt,
        'LOCPOT':Locpot,
        'OSZICAR':Oszicar,
        'OUTCAR':Outcar
    }

class VaspSimulation():
    def __init__(self):
        self.input = VaspInputFiles()
        self.output = VaspOutputFiles()
