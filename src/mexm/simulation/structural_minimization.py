import os,copy
from collections import OrderedDict
from mexm.simulation import LammpsSimulation
from mexm.simulation import StructuralMinimization

class LammpsStructuralMinimization(LammpsSimulation, StructuralMinimization):
    """ Class for LAMMPS structural minimization

    This data class defines additional attributes and methods necessary to 
    interact with the Workflow manager.

    Args:
        task_name(str): unique id for the task name being define
        task_directory(str): the directory where this task will create
            input and output files for LAMMPS
    """

    simulation_type = 'lmps_min_all'
    is_base_class= False
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
    ):
        """ default constructor
        Args:
            task_name(str): unique id for the task name being define
            task_directory(str): the directory where this task will create
                input and output files for LAMMPS

        """
        LammpsSimulation.__init__(
            self,
            name=name,
            simulation_path=simulation_path,
            strucutre_path=structure_path
        )

    def lammps_input_file_to_string(self):
        str_out = "".join([\
                self._lammps_input_initialization_section(),
                self._lammps_input_create_atoms(),
                self._lammps_input_define_potential(),
                self._lammps_input_run_minimization(),
                self._lammps_input_out_section()])
        return(str_out)

    def on_post(self,configuration=None):
        self.__get_results_from_lammps_outputfile()
        LammpsSimulation.on_post(self,configuration=configuration)
    
    def __get_results_from_lammps_outputfile(self):
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

    def _lammps_input_run_minimization(self):
        str_out = (
            '# ---- define settings\n'
            'compute eng all pe/atom\n'
            'compute eatoms all reduce sum c_eng\n'
            '# ---- run minimization\n'            
            'reset_timestep 0\n'
            'fix 1 all box/relax iso 0.0 vmax 0.001\n'
            'thermo 10\n'
            'thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n'
            # 'thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n'
            'min_style cg\n'
            'minimize 1e-25 1e-25 5000 10000\n'
            )
        return str_out

