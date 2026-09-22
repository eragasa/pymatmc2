""" Implementation of AbstractLammpsSimulation

"""
import os, importlib, subprocess, shutil
import inspect

from copy import deepcopy
from collections import OrderedDict

from mexm.io.vasp import Poscar
from mexm.io.lammps import LammpsStructure
from mexm.structure import SimulationCell

from mexm.simulation import Simulation
from mexm.io.eamtools import SetflFile

from mexm.potential import Potential,EamPotential
from mexm.potential import StillingerWeberPotential
from mexm.manager import PotentialManager

from mexm.exception import LammpsSimulationError

class LammpsSimulation(Simulation):
    simulation_type = 'lammps_base'
    is_base_class = True
    results_names = [
        'toten', 'natoms',
        'a11', 'a12', 'a13', 
        'a21', 'a22', 'a23', 
        'a31', 'a32', 'a33',
        'totpress',
        'p11', 'p12', 'p13', 
        'p21', 'p22', 'p23', 
        'p31', 'p32', 'p33'
    ]

    """ Calculates cohesive energy

    This is an abstract data class which defines the attributes and methods
    necessary to interact with the Workflow manager.  A default implementation
    has been created for this class.  This class calculates the energy of
    a simulation cell.  No structural or positions are calculated.

    Args:
        task_name(str): unique id for the task name being defined
        simulation_path(str): the directory where this task will create
            input and output files for LAMMPS.
        structure_path(str):
        simulation_type(str):
        restart(bool)
        fullauto(bool)
    Attributes:
        task_name(str)
        path(str)
        task_type(str)
        is_restart(bool)
        is_fullauto(bool)
        lammps_input_path(str)
        lammps_output_path(str)
        lammps_structure_path(str)
        lammps_eam_filename(str)
        potential(pypospack.potential.Potential): the potential class
        structure_path(str)
        structure(pypospack.io.vasp.Poscar): the structure
        config(:obj:'list' of :obj:'str'): a list of attributes required to
            configure this LAMMPS task.
        config_map(dict):
        potential_map(dict):
        results(dict): results of the simulation
        lammps_bin(str): location of the serial lammps binary
    """

    def __init__(self,
                 name,
                 path,
                 structure_path='POSCAR',
                 bulk_structure_name=None,
                 fullauto=False,
                 use_mpi=False,
                 lammps_bin=None):


        super().__init__(name=name, path=path)

        self._bulk_structure_name = None
        self._is_fullauto = None

        if bulk_structure_name is not None:
            self.bulk_structure_name = bulk_structure_name
        self.is_fullauto = fullauto

        self.configuration =  None
        self.potential = None
        self.structure = None
        self.lammps_script = None

        self.structure_path = structure_path
        self.read_structure_file(path=self.structure_path)
        self.bulk_structure_name = bulk_structure_name

        self.lammps_input_path = os.path.join(self.path, 'lammps.in')
        self.lammps_output_path = os.path.join(self.path, 'lammps.out')
        self.lammps_structure_path = os.path.join(self.path, 'lammps.structure')
        self.lammps_potentialmod_path = os.path.join(self.path, 'potential.mod')
        self.lammps_setfl_path = None

        self.use_mpi = use_mpi
        self.set_lammps_bin(lammps_bin)

        # flowcontrol filename
        self.results = None
        self.results_path = 'mexm.{}.out'.format(name)

        self.process = None
        self.results_processed = None
        # configuration
        self.configuration = OrderedDict()

    def read_structure_file(self, path):
        if path.endswith('.vasp') or path.endswith('POSCAR'):
            self.structure = Poscar()
            self.structure.read(path=path)
        else:
            raise LammpsSimulationError('cannot read structure file: {}'.format(path))

    def set_lammps_bin(self, lammps_bin=None):
        if lammps_bin is not None:
            self.lammps_bin = lammps_bin
        elif self.use_mpi:
            self.lammps_bin = os.environ['LAMMPS_MPI_BIN']
        else:
            try:
                self.lammps_bin = os.environ['LAMMPS_SERIAL_BIN']
            except KeyError:
                self.lammps_bin = None

    @property
    def bulk_structure_name(self):
        return self._bulk_structure_name
    
    @bulk_structure_name.setter
    def bulk_structure_name(self, name):
        if isinstance(name, str):
            self._bulk_structure_name = name
        elif name is None:
            self._bulk_structure_name = name
        else:
            raise TypeError('name must be a str')

    @property
    def is_fullauto(self):
        return self._is_fullauto

    @is_fullauto.setter
    def is_fullauto(self, is_fullauto):
        if isinstance(is_fullauto, bool):
            self._is_fullauto = is_fullauto
        else:
            raise TypeError('is_fullauto must be a string')

    @property
    def parameters(self):
        if self.potential is None:
            return None
        else:
            return self.potential.parameters

    @parameters.setter
    def parameters(self,parameters):
        self.potential.parameters = parameters

    @property
    def symbols(self):
        return self.potential.symbols

    def set_potential_parameters(self,parameters=None):
        """
        Sets the parameters of the potential

        Args:
            parameters(OrderedDict)
        """
        if self.potential is None:
            raise LammpsSimulationError('cannot set potential parameters '
                'because the potential attribute has not been initialized')

        if parameters is None:
            self.parameters = self.configuration['potential']['parameters']
        else:
            assert isinstance(parameters, dict)
            self.parameters = deepcopy(parameters)

    def on_init(self, configuration=None):
        if configuration is not None:
            self.configuration = deepcopy(configuration)

        self.configure_potential(potential=self.configuration['potential'])

        if isinstance(self.potential, EamPotential):
            self.write_eam_potential_file()
        self.set_potential_parameters()

        if self.structure is None:
            if self.structure_path is not None:
                self.read_structure_file(path=self.structure_path)
                assert isinstance(self.structure, SimulationCell)
                
            else:
                msg = (
                    "simulation cannot continue unless either the structure "
                    "can be resolved"
                )
                raise LammpsSimulationError(msg)

        self.update_status()
        if self.is_fullauto:
            self.on_update_status()

    def initialize_potential(self, potential = None):
        """ initialize the potential attribute

        Args:
            potential (dict)
        """

        # process the arguments
        if potential is None:
            # maybe it is in self.configuration
            pass
        elif isinstance(potential, dict):
            # add to the configuration if passed
            self.configuration['potential'] = deepcopy(potential)
        else:
            raise TypeError('potential must be a dict')

        # can't initialize the potential, if it is not in the configuration
        if 'potential' not in self.configuration:
            return
        else:
            if not isinstance(self.configuration['potential'], dict):
                msg = (
                    'LammpsSimulation.configuration[\'potential\') should be '
                    'a dict'
                )
                raise TypeError(msg)

        potential = self.configuration['potential']
        self.configure_potential(potential=potential)

    def on_config(self, configuration=None, results=None):
        if configuration is not None:
            self.configuration = deepcopy(configuration)

        assert isinstance(results, dict)
        self.results = deepcopy(results)

        # initialize the potential
        if 'potential' in self.configuration:
            self.initialize_potential()

        # assign the parameters on the potential
        if isinstance(self.potential, Potential):
            if 'parameters' in self.configuration['potential']:
                parameters_ = self.configuration['potential']['parameters']
                self.potential.parameters = parameters_

        # writing eam potential files
        if type(self.potential) is EamPotential:
            if self.lammps_setfl_path is None:
                self.lammps_setfl_path = "{}.eam.alloy".format(
                        "".join(self.potential.symbols))

            # if setfl_filename_src is set, then we just copy the
            # EAM potential file.
            if all([self.potential.obj_pair is None,
                    self.potential.obj_density is None,
                    self.potential.obj_embedding is None,
                    self.potential.setfl_filename_src is not None]):
                _eam_setfl_filename_src = self.potential.setfl_filename_src
                _eam_setfl_filename_dst = os.path.join(
                        self.path,
                        self.lammps_setfl_path)
                shutil.copyfile(
                        src=_eam_setfl_filename_src,
                        dst=_eam_setfl_filename_dst)
            elif all([self.potential.obj_pair is not None,
                      self.potential.obj_density is not None,
                      self.potential.obj_embedding is not None]):
                pass
            else:
                msg_err = (
                    "EamPotential must be either be parameterized by setting "
                    "pair,density,and embedding formalisms through the "
                    "constructor or a setfl filename must be provided through "
                    "the filename argument\n"
                    "obj_pair:{obj_pair}\n"
                    "obj_density:{obj_pensity}\n"
                    "obj_embedding:{obj_embedding}\n"
                    "setfl_filename:{setfl_filename\n"
                    ).format(
                            obj_pair=str(type(self.potential.obj_pair)),
                            obj_density=str(type(self.potential.obj_density)),
                            obj_embedding=str(type(self.potential.obj_embedding)),
                            setfl_filename=str(self.potential.setfl_filename))
                raise ValueError(msg_err)

        if isinstance(self.potential,EamPotential):
            if self.potential.setfl_filename_src is None:
                if self.lammps_setfl_path is None:
                    self.lammps_setfl_path = '{}.eam.alloy'.format(
                            "".join(self.potential.symbols))
                _eam_setfl_filename_dst = os.path.join(
                        self.path,
                        self.lammps_setfl_path)
                self.write_eam_setfl_file(
                        filename=_eam_setfl_filename_dst)

        self.update_status()
        if self.is_fullauto:
            self.on_update_status()

    def on_ready(self,configuration=None, results=None):
        if configuration is not None:
            self.configuration = deepcopy(configuration)
        self.write_lammps_input_file()
        self.write_structure_file(results=results)
        self.write_potential_file(results=results)
        
        self.run()

        self.update_status()
        if self.is_fullauto:
            self.on_update_status()

    def on_running(self,configuration=None):
        self.update_status()
        if self.is_fullauto:
            self.on_update_status()

    def on_post(self,configuration=None):
        self.results_processed = True
        self.update_status()
        if self.is_fullauto:
            self.on_update_status()

    def on_finished(self,configuration=None):
        # doing nothing here
        self.update_status()
        if self.is_fullauto:
            self.on_update_status()

    def on_error(self,configuration=None):
        raise ValueError()

    def get_conditions_init(self):
        self.conditions_INIT = OrderedDict()
        self.conditions_INIT['simulation_directory_created']\
                = os.path.isdir(self.path)
        return self.conditions_INIT

    def is_potential_initialized(self):
        return isinstance(self.potential, Potential)

    def is_potential_parameters_processed(self):
        if isinstance(self.potential, EamPotential):
            if isinstance(self.potential.setfl_src_path, str):
                return True
            else:
                return all([
                    v is not None for k, v in self.potential.parameters.items()
                ])
        else:
            return all([
                v is not None for k, v in self.potential.parameters.items()
            ])

    def get_conditions_config(self):
        self.conditions_CONFIG = {
            'potential_initialized':self.is_potential_initialized(),
            'parameters_processed':self.is_potential_parameters_processed()
        }
        return self.conditions_CONFIG

    def get_conditions_ready(self):
        self.conditions_READY = OrderedDict()
        return self.conditions_READY

    def get_conditions_running(self):
        self.conditions_RUNNING = OrderedDict()
        self.conditions_RUNNING['process_initialized'] = self.process is not None

    def get_conditions_post(self):
        self.conditions_POST = OrderedDict()

        if self.process is None:
            _process_finished = False
        else:
            _poll_result = self.process.poll()
            if _poll_result is None:
                _process_finished = False
            elif _poll_result == 0:
                _process_finished = True
            elif _poll_result == 1:
                if self.conditions_ERROR is None:
                    self.conditions_ERROR=OrderedDict()

                lammps_out_fn = os.path.join(self.path,'lammps.out')
                with open(lammps_out_fn) as f:
                    lines = f.readlines()
                    last_line = lines[len(lines)-1]

                if "Neighbor list overflow" in last_line:
                    m  = "Neighbor list overflow"
                else:
                    m = "Lammps excited with status {}.  If running an EAM "
                    m += "potential this is most likely caused by an out-of-index "
                    m += "exception because the electron density is too high when "
                    m += "evaluating the embedding function.  The code for modifying "
                    m += "max_rho for the embedding function is in "
                    m += "pypospack.potential.EamPotential"

                    m = m.format(_poll_result)

                self.conditions_ERROR['lmps_bin_err'] = m
                raise LammpsSimulationError(m, parameters=self.potential.parameters)
            else:
                if self.conditions_ERROR is None:
                    self.conditions_ERROR= OrderedDict()
                m = 'Lammps exited with status {}.'.format(_poll_result)
                self.conditions_ERROR['lmps_bin_err'] = m
                raise LammpsSimulationError(m, parameters=self.potential.parameters)

        self.conditions_POST['process_finished'] = _process_finished


    def get_conditions_finished(self):
        self.conditions_FINISHED = OrderedDict()
        self.conditions_FINISHED['is_processed'] = self.results_processed

    def get_conditions_error(self):
        if self.conditions_ERROR is None:
            self.conditions_ERROR = OrderedDict()

    def restart(self):
        raise NotImplementedError

    def run(self):

        _lammps_bin = self.lammps_bin
        cmd_str = '{} -i lammps.in > lammps.out'.format(_lammps_bin)

        # change context directory

        _cwd = os.getcwd()
        os.chdir(self.path)

        # https://stackoverflow.com/questions/4789837/how-to-terminate-a-python-subprocess-launched-with-shell-true/4791612#4791612
        #self.process = subprocess.Popen(
        #        cmd_str,
        #        shell=True,
        #        cwd=self.path,
        #        preexec_fn=os.getpgrp)
        #        #preexec_fn=os.setsid)
        #self.process_info = OrderedDict()

        # NEW CODE
        self.process = subprocess.Popen("exec " + cmd_str, shell=True)
        os.chdir(_cwd)

    def write_eam_setfl_file(self,filename):
        _setfl_dst_filename = filename
        _Nr = self.configuration['potential']['N_r']
        _Nrho = self.configuration['potential']['N_rho']
        _rmax = self.configuration['potential']['r_max']
        _rhomax = self.configuration['potential']['rho_max']
        _rcut = self.configuration['potential']['r_cut']
        _parameters = self.configuration['parameters']

        _a0 = self.configuration['potential']['a0']
        _latt_type = self.configuration['potential']['lattice_type']

        _rcut = self.potential.determine_r_max(a0=_a0,latt_type='fcc')
        _rhomax = self.potential.determine_rho_max(a0=_a0,latt_type='fcc')
        self.configuration['potential']['rho_max'] = _rhomax
        self.configuration['potential']['rcut'] = _rcut

        self.potential.write_setfl_file(
                filename=_setfl_dst_filename,
                symbols=self.potential.symbols,
                Nr=_Nr,
                rmax=_rmax,
                rcut=_rcut,
                Nrho=_Nrho,
                rhomax=_rhomax,
                parameters=_parameters)

    def write_potential_file(self, results=None):
        if self.potential is None:
            raise ValueError

        self.potential.parameters = deepcopy(self.configuration['potential']['parameters'])

        _setfl_dst_filename = None

        # <-------- FOR EAM POTENTIALS
        if isinstance(self.potential,EamPotential):
            _symbols = "".join(self.potential.symbols)
            _filename = "{}.eam.alloy".format(_symbols)
            _setfl_dst_filename = os.path.join(self.path,_filename)
            _str_out = self.potential.lammps_potential_section_to_string(
                setfl_dst_filename=_setfl_dst_filename)

        # <-------- FOR STILLINGER WEBER POTENTIALS
        elif isinstance(self.potential,StillingerWeberPotential):
            # define the filename --- SiO.parameters, Si.parameters
            _symbols_str = "".join(self.potential.symbols)
            _p_fname = "{}.parameters".format(_symbols_str)

            # set the name of the output file
            self.potential.lmps_parameter_filename = _p_fname

            # get the string of potential.mod
            _str_out = self.potential.lammps_potential_section_to_string()

            # write out the potential parameter file
            _str_lmps_params = self.potential.lammps_parameter_file_to_string()

            _p_fname_dst = os.path.join(self.path,_p_fname)
            with open(_p_fname_dst,'w') as f:
                f.write(_str_lmps_params)

        else:
            _str_out = self.potential.lammps_potential_section_to_string()

        _str_out += "\n"

        # <-------- EWALD CHARGE SUMMATION METHOD
        if self.potential.is_charge:
            _str_out += "kspace_style pppm 1.0e-5\n"
            _str_out += "\n"

        # <-------- TREATMENT OF NEAREST NEIGHBORS
        _str_out += "neighbor 1.0 bin\n"
        _str_out += "neigh_modify every 1 delay 0 check yes\n"

        # <-------- WRITE POTENTIAL.MOD TO FILESYSTEM
        _lammps_potentialmod_path = os.path.join(
                self.path,
                self.lammps_potentialmod_path)
        with open(_lammps_potentialmod_path,'w') as f:
            f.write(_str_out)

    def write_lammps_input_file(self,filename='lammps.in'):
        """ writes LAMMPS input file

        Args:
            filename (str): name of the input file for LAMMPS. Default is
                'lammps.in'.
        """
        str_out = self.lammps_input_file_to_string()
        filename = os.path.join(self.path,filename)
        with open(filename,'w') as f:
            f.write(str_out)

    def get_atom_style(self):
        if self.potential.is_charge:
            atom_style = 'charge'
        else:
            atom_style = 'atomic'
        return atom_style

    def modify_structure_file(self, results=None):
        pass

    def write_structure_file(self, lammps_structure_path=None, results=None):
        if lammps_structure_path is not None:
            self.lammmps_structure_path = lammps_structure_path

        kwargs = {
            'path':os.path.join(self.path, self.lammps_structure_path),
            'atom_style':self.get_atomic_style()
        }

        self.lammps_structure = LammpsStructure.initialize_from_mexm(self.structure)
        assert isinstance(self.lammps_structure, LammpsStructure)
        self.modify_structure_file(results=results)
        self.lammps_structure.write(**kwargs)

    def get_atomic_style(self):
        if self.potential.is_charge:
            return 'charge'
        else:
            return 'atomic'

    def modify_structure(self, new_info):
        if new_info[0] == 'a0':
            self.structure.a0 = new_info[1]


    def configure_potential(self,potential=None):
        """
            Args:
                potential(dict): potential is a dictionary which has the
                    necessary keywords to configure the object.

            Notes:
                For buckingham potential,
                    potential['name'] = 'buckingham'
                    potential['symbols'] = ['Mg', 'O']
                For eam potentials,
                    potential['potential_type'] = 'eam'
                    potential['pair_type'] = 'morse'
                    potential['density_type'] = 'eam_exp_dens'
                    potential['embedding_type'] = 'eam_universal_embedding'
        """

        # checking arguments
        if potential is None:
            pass
        elif isinstance(potential, dict):
            pass
        else:
            msg = 'potential argument must either be None or a dictionary'
            raise TypeError(msg)


        if isinstance(potential, dict):
            self.configuration['potential']=deepcopy(potential)

        if potential['name'] == 'eam':
            potential_args = inspect.getfullargspec(PotentialManager.get_potential_by_name)
            potential_dict = {}
            for k in potential_args.args:
                try:
                    potential_dict[k] = self.configuration['potential'][k]
                except KeyError:
                    if k == 'potential_name':
                        potential_dict[k] = self.configuration['potential']['name']

        else:
            potential_args = inspect.getfullargspec(PotentialManager.get_potential_by_name)
            potential_dict = {}
            for k in potential_args.args:
                try:
                    potential_dict[k] = self.configuration['potential'][k]
                except KeyError:
                    if k == 'potential_name':
                        potential_dict[k] = self.configuration['potential']['name']
                    elif k in ['pair_type', 'density_type', 'embedding_type']:
                        'only necessary for EAM'
                    else:
                        raise
        self.potential = PotentialManager.get_potential_by_name(**potential_dict)


    def lammps_input_file_to_string(self):
        """ string for the LAMMPS input file """

        str_out = "".join([\
                self._lammps_input_initialization_section(),
                self._lammps_input_create_atoms(),
                self._lammps_input_define_potential(),
                self._lammps_input_run_minimization(),
                self._lammps_input_out_section()])
        return(str_out)

    # private functions for building lammps input files
    def _lammps_input_initialization_section(self):
        if self.potential.is_charge:
            _atom_style = 'charge'
        else:
            _atom_style = 'atomic'
        str_out = (
            '# ---- initialize simulations\n'
            'clear\n'
            'units metal\n'
            'dimension 3\n'
            'boundary p p p\n'
            'atom_style {atom_style}\n'
            'atom_modify map array\n'
        ).format(
            atom_style=_atom_style
        )
        return str_out

    def _lammps_input_create_atoms(self):
        """ creates create atoms section of the LAMMPS input file """

        _structure_path = os.path.join(
                self.lammps_structure_path)

        # this provides LAMMPS with the structure
        str_out = (
            '# ---- create atoms\n'
            'read_data {structure_path}\n'
        ).format(
                structure_path=_structure_path
        )

        return str_out

    def _lammps_input_define_potential(self):
        str_out = (
            '# ---- define interatomic potential\n'
            'include potential.mod\n')
        return str_out

    def _lammps_input_run_minimization(self):
        str_out = (
            '# ---- define settings\n'
            'compute eng all pe/atom\n'
            'compute eatoms all reduce sum c_eng\n'
            '# ---- run minimization\n'
            'reset_timestep 0\n'
            'fix 1 all box/relax aniso 0.0 vmax 0.001\n'
            'thermo 10\n'
            'thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n'
            # 'thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n'
            'min_style cg\n'
            'minimize 1e-25 1e-25 5000 10000\n'
            )
        return str_out

    def _lammps_input_out_section(self):
        str_out = (
            '# ---- define output variables ----\n'
            'variable natoms equal "count(all)"\n'
            'variable tot_energy equal "c_eatoms"\n'
            'variable a11 equal "xhi-xlo"\n'
            'variable a22 equal "yhi-ylo"\n'
            'variable a33 equal "zhi-zlo"\n'
            'variable tilt_xy equal "xy"\n'
            'variable tilt_xz equal "xz"\n'
            'variable tilt_yz equal "yz"\n'
            'variable tot_press equal "press"\n'
            'variable press_xx equal "pxx"\n'
            'variable press_yy equal "pyy"\n'
            'variable press_zz equal "pzz"\n'
            'variable press_xy equal "pxy"\n'
            'variable press_xz equal "pxz"\n'
            'variable press_yz equal "pyz"\n'
            '\n'
            '# ---- output ----\n'
            'print \"pypospack:output_section:begin\"\n'
            'print \"toten = ${tot_energy}\"\n'
            'print \"natoms = ${natoms}"\n'
            'print \"a11 = ${a11}\"\n'
            'print \"a22 = ${a22}\"\n'
            'print \"a33 = ${a33}\"\n'
            'print \"a12 = ${tilt_xy}\"\n'
            'print \"a13 = ${tilt_xz}\"\n'
            'print \"a23 = ${tilt_yz}\"\n'
            'print \"totpress = ${tot_press}\"\n'
            'print \"p11 = ${press_xx}\"\n'
            'print \"p22 = ${press_yy}\"\n'
            'print \"p33 = ${press_zz}\"\n'
            'print \"p12 = ${press_xy}\"\n'
            'print \"p13 = ${press_xz}\"\n'
            'print \"p23 = ${press_yz}\"\n'
            'print \"pypospack:output_section:done\"\n'
            'print \"pypospack:lammps_sim:done\"\n'
                  )
        return str_out

    def __get_results_from_lammps_outputfile(self):
        path_ = os.path.join(self.path, 'lammps.out')
        with open(path_, 'r') as f:
            lines = f.readlines()

        results_names_ = self.results_names
        results_ = OrderedDict()

        for i, line in enumerate(lines):
            for name in results_names_:
                if line.startswith('{} = '.format(name)):
                    results_[name] = float(line.split('=')[1].strip())

                if line.startswith('ERROR:'):
                    print('name:{}'.format(name))
                    print('line:{}'.format(line.strip))
                    raise NotImplementedError

        # store the results in the obj attribute
        self.results = OrderedDict()
        for k in self.results_names:
            self.results['{}.{}'.format(self.name, k)] = results_[k]
    def write_eam_potential_file(self):
        if self.lammps_setfl_path is None:
            self.lammps_setfl_path = "{}.eam.alloy".format(
                    "".join(self.potential.symbols))
            if self.lammps_setfl_path is None:
                self.lammps_setfl_path = "{}.eam.alloy".format(
                        "".join(self.potential.symbols))

            # if setfl_filename_src is set, then we just copy the
            # EAM potential file.

            if all([self.potential.obj_pair is None,
                    self.potential.obj_density is None,
                    self.potential.obj_embedding is None,
                    isinstance(self.potential.setfl_filename_src,str)]):
                eam_setfl_filename_src = self.potential.setfl_filename_src
                eam_setfl_filename_dst = os.path.join(
                        self.path,
                        self.lammps_setfl_path)
                try:
                    shutil.copyfile(
                            src=eam_setfl_filename_src,
                            dst=eam_setfl_filename_dst)
                except FileNotFoundError as e:
                    msg  = 'o.potential.obj_pair:{}\n'.format(self.potential.obj_pair)
                    msg += 'o.potential.obj_density:{}\n'.format(self.potential.obj_density)
                    msg += 'o.potential.obj_embedding:{}\n'.format(self.potential.obj_embedding)
                    msg += 'eam_setfl_filename_src:{}\n'.format(eam_setfl_filename_src)
                    msg += 'eam_setfl_filename_dst:{}\n'.format(eam_setfl_filename_dst)
                    raise

            elif all([self.potential.obj_pair is not None,
                      self.potential.obj_density is not None,
                      self.potential.obj_embedding is not None]):
                pass
            else:
                msg_err = (
                    "EamPotential must be either be parameterized by setting "
                    "pair,density,and embedding formalisms through the "
                    "constructor or a setfl filename must be provided through "
                    "the filename argument\n"
                    "obj_pair:{obj_pair}\n"
                    "obj_density:{obj_pensity}\n"
                    "obj_embedding:{obj_embedding}\n"
                    "setfl_filename:{setfl_filename\n"
                    ).format(
                            obj_pair=str(type(self.potential.obj_pair)),
                            obj_density=str(type(self.potential.obj_density)),
                            obj_embedding=str(type(self.potential.obj_embedding)),
                            setfl_filename=str(self.potential.setfl_filename))
                raise ValueError(msg_err)
