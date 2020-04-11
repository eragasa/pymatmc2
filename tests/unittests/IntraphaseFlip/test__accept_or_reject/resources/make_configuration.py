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

configuration['calculator'] = {
    'calculator':'vasp',
    'simulation_type':'vasp_min_all'
}

configuration['job_submission'] = {
    'hpc_type':'torque',
    'hpc_config_path':'osc.pitzer'
}

configuration['atomic_configuration'] = {
    'molar_fraction_total':{'Au':0.50, 'Pt':0.50},
    'simulation_cells':{
        'fcc1':{
            'incar':  os.path.join('fcc1', os.path.dirname(os.path.abspath(__file__)), 'INCAR'),
            'kpoints':os.path.join('fcc1', os.path.dirname(os.path.abspath(__file__)), 'KPOINTS'),
            'poscar': os.path.join('fcc1', os.path.dirname(os.path.abspath(__file__)), 'POSCAR'),
            'potcar': os.path.join('fcc1', os.path.dirname(os.path.abspath(__file__)), 'POTCAR')
        },
        'fcc2':{
            'incar':  os.path.join('fcc2', os.path.dirname(os.path.abspath(__file__)), 'INCAR'),
            'kpoints':os.path.join('fcc2', os.path.dirname(os.path.abspath(__file__)), 'KPOINTS'),
            'poscar': os.path.join('fcc2', os.path.dirname(os.path.abspath(__file__)), 'POSCAR'),
            'potcar': os.path.join('fcc2', os.path.dirname(os.path.abspath(__file__)), 'POTCAR')
        }
    }
}
configuration['max_iterations'] = 10
configuration['mutation_weights'] = {
    'intraphase_swap':0.5,
    'intraphase_flip':0.5
}
configuration['environment_variables'] = {
    'temperature':400.0,
    'pressure':0.0
}
path = "pymatmc2.config"

obj = Pymatmc2Configuration()
obj.configure_from_dict(configuration)
obj.write(path=path)

assert os.path.isfile(path)
