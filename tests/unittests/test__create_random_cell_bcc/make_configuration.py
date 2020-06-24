import os
from pymatmc2 import Pymatmc2Configuration

configuration = {}
configuration['hpc_manager'] = {}
configuration['hpc_manager']['type'] = 'torque'
configuration['hpc_manager']['configuration'] = {
    'account':'PAA0028',
    'walltime':1,
    'n_nodes':4,
    'ppn':40,
    'jobname':'Si_job_name',
    'errpath':'job.err',
    'stdpath':'job.out',
    'modules':['intel/19.0.5', 'intelmpi/2019.3'],
    'cmd':'mpiexec $VASP_STD_BIN > vasp.out'
}

# configuration of results output
configuration['results'] = {}
configuration['results']['dir'] = 'results'
configuration['results']['tar_path'] = 'pymatmc2.results.tar'
configuration['results']['file_path'] = 'pymatmc2.results'

configuration['calculator'] = {
    'calculator':'vasp',
    'simulation_type':'vasp_min_all'
}

configuration['job_submission'] = {
    'hpc_type':'torque',
    'hpc_config_path':'osc.pitzer'
}

configuration['atomic_configuration'] = {
    'total_concentration':{
       'Hf':1/4,
       'Zr':1/4,
       'Ta':1/4,
       'Nb':1/4
    },
    'simulation_cells':{
        'bcc1':{
            'incar':  os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc1', 'INCAR'),
            'kpoints':os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc1', 'KPOINTS'),
            'poscar': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc1', 'POSCAR'),
            'potcar': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc1', 'POTCAR')
        },
        'bcc2':{
            'incar':  os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc2', 'INCAR'),
            'kpoints':os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc2', 'KPOINTS'),
            'poscar': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc2', 'POSCAR'),
            'potcar': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc2', 'POTCAR')
        },
        'bcc3':{
            'incar':  os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc3', 'INCAR'),
            'kpoints':os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc3', 'KPOINTS'),
            'poscar': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc3', 'POSCAR'),
            'potcar': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc3', 'POTCAR')
        },
        'bcc4':{
            'incar':  os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc4', 'INCAR'),
            'kpoints':os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc4', 'KPOINTS'),
            'poscar': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc4', 'POSCAR'),
            'potcar': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bcc4', 'POTCAR')
        }
    }
}
configuration['max_iterations'] = 1000
configuration['mutation_weights'] = {
    'intraphase_flip':1.0
}
configuration['environment_variables'] = {
    'temperature':800.0,
    'pressure':0.0
}
path = "pymatmc2.config"

obj = Pymatmc2Configuration()
obj.configure_from_dict(configuration)
obj.write(path=path)

assert os.path.isfile(path)