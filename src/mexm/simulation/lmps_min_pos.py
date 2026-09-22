import os,copy
from collections import OrderedDict
from mexm.io.vasp import Poscar
from mexm.simulation import LammpsSimulation
from mexm.simulation import PositionMinimization

class LammpsPositionMinimization(LammpsSimulation, PositionMinimization):
    """ Class for LAMMPS structural minimization

    This data class defines additional attributes and methods necessary to
    interact with the Workflow manager.

    Attributes:
        bulk_structure_name (str)
        bulk_structure_path (str)
        config_map
    """
    simulation_type = 'lammps_min_pos'
    is_base_class = False
    results_name = [
        'toten', 'natoms',
        # lattice elements
        'a11', 'a12', 'a13', 
        'a21', 'a22', 'a23', 
        'a31', 'a32', 'a33',
        'totpress',
        # pressure tensor
        'p11', 'p12', 'p13', 
        'p21', 'p22', 'p23', 
        'p31', 'p32', 'p33'
    ]

    def __init__(
        self,
        name,
        simulation_path,
        structure_path,
        bulk_structure_name=None
    ):
        """ default constructor

        Args:
            name (str)
            simulation_path (str)
            structure_path (str)
            bulk_structure_name (str)
        """

        self.bulk_structure_name = None
        self.bulk_structure_path = None
        self.bulk_structure_lattice = None


        LammpsSimulation.__init__(self,
                                  name=name,
                                  path=simulation_path,
                                  structure_path=structure_path,
                                  bulk_structure_name=bulk_structure_name)

    def lammps_input_file_to_string(self):
        str_out = "".join([\
                self._lammps_input_initialization_section(),
                self._lammps_input_create_atoms(),
                self._lammps_input_define_potential(),
                self._lammps_input_run_minimization(),
                self._lammps_input_out_section()])
        return(str_out)

    def on_init(self,configuration=None,results=None):
        LammpsSimulation.on_init(self,configuration=configuration)

        if 'bulk_structure' in configuration:
            self.bulk_structure_name = configuration['bulk_structure']
            self.bulk_structure_path = configuration['bulk_structure_filename']
            self.bulk_structure_lattice = OrderedDict()

            _lattice_parameter_variables = [
                    'lmps_min_all.a11',
                    'lmps_min_all.a12',
                    'lmps_min_all.a13',
                    'lmps_min_all.a21',
                    'lmps_min_all.a22',
                    'lmps_min_all.a23',
                    'lmps_min_all.a31',
                    'lmps_min_all.a32',
                    'lmps_min_all.a33']

            self.bulk_lattice_components = OrderedDict()
            for k in _lattice_parameter_variables:
                _k = '{}.{}'.format(self.bulk_structure_name,k)
                self.bulk_lattice_components[_k] = None

    def on_config(self,configuration,results=None):
        LammpsSimulation.on_config(self,configuration=configuration,results=results)

    def get_conditions_ready(self):
        LammpsSimulation.get_conditions_ready(self)
        a0_name = '{}.{}.a0'.format(self.bulk_structure_name, 'lmps_min_all')
        try:
            self.conditions_READY['bulk_structure_results_available'] \
                = a0_name in self.results
        except TypeError as e:
            if self.results is None:
                self.conditions_READY['bulk_structure_results_available'] \
                    = False
            else:
                raise
        return self.conditions_READY

    def on_ready(self,configuration=None,results=None):
        self.results = deepcopy(results)
        LammpsSimulation.on_ready(self,
                                  configuration=configuration,
                                  results=self.results)

    def on_post(self,configuration=None):
        self.__getresults__from_lammps_outputfile()
        LammpsSimulation.on_post(self,configuration=configuration)

    def modify_structure_file(self, results):
        a0_name = '{}.{}.a0'.format(self.bulk_structure_name, 'lmps_min_all')
        self.lammps_structure.a0 = results[a0_name]


    def __getresults__from_lammps_outputfile(self):
        _filename = os.path.join(
                self.path,
                'lammps.out')
        with open(_filename,'r') as f:
            lines = f.readlines()

        _variables = [
                'tot_energy',
                'num_atoms',
                'a11','a12','a13','a22','a23','a33',
                'tot_press',
                'pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz',
                ]
        results_ = OrderedDict()

        for i,line in enumerate(lines):
            for name in _variables:
                if line.startswith('{} = '.format(name)):
                    results_[name] = float(line.split('=')[1].strip())

                if line.startswith('ERROR:'):
                    print('name:{}'.format(name))
                    print('line:{}'.format(line.strip))
                    raise NotImplementedError

        name_ = self.name
        self.results = OrderedDict()
        self.results['{}.{}'.format(name_,'toten')] = results_['tot_energy']
        self.results['{}.{}'.format(name_,'natoms')] = results_['num_atoms']
        # this only works for orthogonal cells
        self.results['{}.{}'.format(name_,'a11')] = results_['a11']
        self.results['{}.{}'.format(name_,'a12')] = results_['a12']
        self.results['{}.{}'.format(name_,'a13')] = results_['a13']
        self.results['{}.{}'.format(name_,'a21')] = 0
        self.results['{}.{}'.format(name_,'a22')] = results_['a22']
        self.results['{}.{}'.format(name_,'a23')] = results_['a23']
        self.results['{}.{}'.format(name_,'a31')] = 0
        self.results['{}.{}'.format(name_,'a32')] = 0
        self.results['{}.{}'.format(name_,'a33')] = results_['a33']
        self.results['{}.{}'.format(name_,'totpress')] = results_['tot_press']
        self.results['{}.{}'.format(name_,'p11')] = results_['pxx']
        self.results['{}.{}'.format(name_,'p12')] = results_['pxy']
        self.results['{}.{}'.format(name_,'p13')] = results_['pxz']
        self.results['{}.{}'.format(name_,'p21')] = results_['pxy']
        self.results['{}.{}'.format(name_,'p22')] = results_['pyy']
        self.results['{}.{}'.format(name_,'p23')] = results_['pyz'] #pyz=pzy
        self.results['{}.{}'.format(name_,'p31')] = results_['pxz'] #pxz=pzx
        self.results['{}.{}'.format(name_,'p32')] = results_['pyz']
        self.results['{}.{}'.format(name_,'p33')] = results_['pzz']

    def _lammps_input_run_minimization(
        self,
        etol = 1e-20,
        ftol = 1e-20,
        maxiter = 1000,
        maxeval = 10000):
        str_out = (
            '# ---- define settings\n'
            'compute eng all pe/atom\n'
            'compute eatoms all reduce sum c_eng\n'
            '\n'
            '# ---- run minimization\n'
            'reset_timestep 0\n'
            'thermo 1\n'
            'thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n'
            'min_style cg\n'
            'minimize 1e-20 1e-20 1000 100000\n'
            '\n'
            )
        return str_out
