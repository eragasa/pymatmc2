from collections import OrderedDict
import copy
import numpy as np
from mexm.potential import (Potential,
                            PairPotential,
                            EamDensityFunction,
                            EamEmbeddingFunction)
from mexm.exception import MexmException

class MexmPotentialError(Exception): pass

class EamPotential(Potential):
    potential_type = 'eam'
    is_base_potential = False
    is_charge = False

    MEXM_EAM_PAIR_FORMAT = "pair_{symbol1}{symbol2}_{parameter_name}"
    MEXM_EAM_DENS_FORMAT = "dens_{symbol}{parameter_name}"
    MEXM_EAM_EMBED_FORMAT = "embed_{symbol}_{parameter_name}"
    """embedded energy method potential
    This class is for the modelling of an EAM potential
    Args:
       symbols(list of str):
       func_pairpotential(str):
       func_density(str):
       func_embedding(str):
       filename(str):

    Attributes:
       symbols(list of str): the list of symbols
       obj_pair(OrderedDict of PairPotential)
       obj_density(OrderedDict of EamDensityFunction)
       obj_embedding(OrderedDict of EamEmbeddingFunction)
       N_r(int): number of radial points, the distance between two atoms
       r_max(float): the maximum distance between two atoms
       r_cut(float): the cutoff distance between two atoms
       N_rho(int): the number of points in the electron density evaluation
       rho_max(float): the maximum electron density
    """
    def __init__(self,
            symbols,
            pair_type=None,
            density_type=None,
            embedding_type=None):

        self.pair = {
            'type':pair_type
        }
        self.density = {
            'type':density_type
        }
        self.embedding = {
            'type':embedding_type
        }

        super().__init__(symbols=symbols, is_charge=False)

        self.r = None
        self.N_r = None
        self.r_max = None

        self.rho = None
        self.N_rho = None
        self.rho_max = None


        self.setfl_filename_src = None
        self.setfl_filename_dst = "{symbols}.eam.alloy".format(
                symbols="".join(self.symbols))
        self.setfl = None

        self._initialize_pair_potential(pair_type)
        self._initialize_density_function(density_type)
        self._initialize_embedding_function(embedding_type)

    @property
    def pair_type(self): return self.pair['type']

    @property
    def density_type(self): return self.embedding['type']

    @property
    def embedding_type(self): return self.embedding['type']

    @classmethod
    def get_parameter_names(cls,
                            symbols,
                            pair_type,
                            density_type,
                            embedding_type):
        from mexm.manager import PotentialManager

        parameter_names = ['rhocut', 'rcut']
        potential_types = [k.potential_type for k in PotentialManager.get_potential_types()]
        assert pair_type in potential_types
        assert density_type in potential_types
        assert embedding_type in potential_types

        for potential in PotentialManager.get_potential_types():
            if potential.potential_type == pair_type:
                parameter_names += [
                    'pair_{}'.format(k)
                    for k in potential.get_parameter_names(symbols)
                ]
        for potential in PotentialManager.get_potential_types():
            if potential.potential_type == density_type:
                parameter_names += [
                    'dens_{}'.format(k)
                    for k in potential.get_parameter_names(symbols)
                ]
        for potential in PotentialManager.get_potential_types():
            if potential.potential_type == embedding_type:
                parameter_names += [
                    'embed_{}'.format(k)
                    for k in potential.get_parameter_names(symbols)
                ]

        return parameter_names

    def _initialize_parameter_names(self):
        kwargs = {
            'symbols':self.symbols,
            'pair_type':self.pair['type'],
            'density_type':self.density['type'],
            'embedding_type':self.embedding['type']
        }
        self.parameter_names = self.get_parameter_names(**kwargs)

    def _initialize_pair_potential(self, pair_type):
        from mexm.manager import PotentialManager

        assert isinstance(pair_type, str)
        self.pair['type'] = pair_type
        self.pair['obj'] = PotentialManager.get_potential_by_name(
                pair_type,
                self.symbols
        )

    def _initialize_density_function(self, density_type):
        from mexm.manager import PotentialManager

        assert isinstance(density_type, str)
        self.density['type'] = density_type
        self.density['obj'] = PotentialManager.get_potential_by_name(
                density_type,
                self.symbols
        )

    def _initialize_embedding_function(self, embedding_type):
        from mexm.manager import PotentialManager

        assert isinstance(embedding_type, str)
        self.embedding['type'] = embedding_type
        self.embedding['obj'] = PotentialManager.get_potential_by_name(
                embedding_type,
                self.symbols
        )

    def lammps_potential_section_to_string(self,setfl_dst_filename):
        """provide string for the potential section
        Args:
            setfl_dst_filename(str):the path that LAMMPS should read for the SETFL file
        Returns:
            str: the string of the lammps potential section
        """

        # set masses
        str_out = ''
        for i,s in enumerate(self.symbols):
            str_out += "mass {} {}\n".format(i+1,self._get_mass(s))
        str_out += "\n"

        # set groups
        for i,s in enumerate(self.symbols):
            str_out += "group {} type {}\n".format(s,i+1)
        str_out += "\n"

        str_out += "pair_style eam/alloy\n"
        str_out += "pair_coeff * * {setfl_dst_filename} {str_symbols}\n".format(
                    setfl_dst_filename=setfl_dst_filename,
                    str_symbols=" ".join(self.symbols))

        return str_out


    def determine_r_max(self,a0,latt_type):
        _a0 = a0

        if latt_type == 'fcc':
            _d_1NN = 0.707 * _a0
            _d_2NN = 1.000 * _a0
            _d_3NN = 1.225 * _a0
            _d_4NN = 1.414 * _a0
            _d_5NN = 1.581 * _a0

        # r_max should be between the 3NN and 4NN
        _rcut = 0.5 * (_d_3NN + _d_4NN)
        return _rcut

    def determine_rho_max(self,a0,latt_type):
        _a0 = a0

        #extract parameters for the density function
        _parameters = OrderedDict()
        for k,v in self.parameters.items():
            if k.startswith('d_'):
                _parameter_name = k[2:]
                _parameters[_parameter_name] = v

        if latt_type == 'fcc':
            _d_1NN = 0.707 * _a0
            _d_2NN = 1.000 * _a0
            _d_3NN = 1.225 * _a0
            _d_O = 0.866 * _a0
            _d_T = 0.433 * _a0

            _natoms_1NN = 12
            _natoms_2NN = 6
            _natoms_O = 4
            _natoms_T = 2

        for s in self.symbols:
            _rhomax = _natoms_1NN * self.obj_density.evaluate(_d_1NN,_parameters)[s]
            _rhomax += _natoms_2NN * self.obj_density.evaluate(_d_2NN,_parameters)[s]
            _rhomax += _natoms_O * self.obj_density.evaluate(_d_O,_parameters)[s]
            _rhomax += _natoms_T * self.obj_density.evaluate(_d_T,_parameters)[s]

        _rhomax = float(_rhomax)
        return _rhomax

    def write_setfl_file(self,filename,symbols,
            Nr,rmax,rcut,
            Nrho,rhomax,
            parameters):
        assert type(filename) is str
        assert type(Nr) is int
        assert type(rmax) in [int,float]
        assert type(rcut) in [int,float]
        assert type(Nrho) is int
        assert type(rhomax) in [int,float]

        r = rmax * np.linspace(1,Nr,Nr)/Nr
        rho = rhomax * np.linspace(1,Nrho,Nrho)/Nrho

        self.evaluate(
                r=r,
                rho=rho,
                rcut=rcut,
                parameters=parameters)

        setfl_file = EamSetflFile()
        setfl_file.write(
                filename=filename,
                symbols=symbols,
                r=self.r,
                rho=self.rho,
                rcut=self.r_cut,
                pair=self.pair,
                density=self.density,
                embedding=self.embedding)

    def evaluate(self,r,rho,rcut,parameters):
        assert isinstance(r,np.ndarray)
        assert isinstance(rho,np.ndarray)
        assert type(rcut) in [float,int,type(None)]
        assert type(parameters) in [OrderedDict,dict]

        for p in self.parameters:
            self.parameters[p] = parameters[p]

        self.rcut = rcut
        self.r = copy.deepcopy(r)
        self.rho = copy.deepcopy(rho)
        self.evaluate_pair(r=r,rcut=rcut)
        self.evaluate_density(r=r,rcut=rcut)
        self.evaluate_embedding(rho=rho, rhocut=rcut)

    def _log(self,msg):
        print(msg)



    def evaluate_pair(self,r,parameters=None,rcut=None):
        assert isinstance(r,np.ndarray)

        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]

        # pair potential parameters are prepended with 'p_'
        p_params = [p for p in self.parameters if p.startswith('p_')]

        # subselect the parameters required for the pair potential
        _parameters = OrderedDict()
        for p in p_params:
            _parameter_name = p.partition('_')[2]
            _parameters[_parameter_name] = self.parameters[p]

        self.obj_pair.evaluate(
                r=r,
                parameters=_parameters,
                r_cut=rcut)

        self.pair = copy.deepcopy(self.obj_pair.potential_evaluations)

    def evaluate_density(self,
            r,
            rcut=None,
            parameters=None):
        assert isinstance(r,np.ndarray) or isinstance(r,float)
        assert isinstance(rcut,float) or rcut is None
        assert isinstance(parameters,dict) or parameters is None

        #<--- check arguments of the function
        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]

        #<--- grab the parameters of the density function
        d_params = [p for p in self.parameters if p.startswith('d_')]

        _parameters = OrderedDict()
        for p in d_params:
            _parameter_name = p.partition('_')[2]
            _parameters[_parameter_name] = self.parameters[p]

        _dens_eval= self.obj_density.evaluate(
                r=r,
                parameters=_parameters,
                r_cut=rcut)

        #<--- set the internal attribute
        self.density = copy.deepcopy(_dens_eval)

        return _dens_eval

    def evaluate_embedding(self,
            rho=None,
            parameters=None):
        #<--- check arguments of the function
        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]
        if rho is not None:
            if isinstance(rho,np.ndarray):
                self.rho = np.copy(rho)
            else:
                raise ValueError("r must be a numpy.ndarray")
            self.rho = np.copy(rho)
        #<--- grab the parameters for the embedding function
        e_params = [p for p in self.parameters if p.startswith('e_')]

        _parameters = OrderedDict()
        for p in e_params:
            _parameter_name = p.partition('_')[2]
            _parameters[_parameter_name] = self.parameters[p]

        #<--- for solving embedding function implicitly from the RoseEquationOfState
        if type(self.obj_embedding) == RoseEquationOfStateEmbeddingFunction:
            if parameters is not None:
                self.obj_embedding.evaluate(
                    rho=self.rho,
                    parameters=parameters, # pass in all parameters not just a subset
                    o_pair=self.obj_pair,
                    o_density=self.obj_density)
            else:
                self.obj_embedding.evaluate(
                    rho=self.rho,
                    r=self.r,
                    parameters=self.parameters, # pass in all parameters not just a subset
                    o_pair=self.obj_pair,
                    o_density=self.obj_density)

        else:
            self.obj_embedding.evaluate(
                rho=self.rho,
                parameters=_parameters)

        #<--- set the internal attribute
        self.embedding = copy.deepcopy(self.obj_embedding.embedding_evaluations)
