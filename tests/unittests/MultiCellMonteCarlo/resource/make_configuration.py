import os
from pymatmc2 import Pymatmc2Configuration

configuration = {}
configuration['hpc_manager'] = {}
configuration['hpc_manager']['type'] = 'torque'
configuration['hpc_manager']['configuration'] = {
    'account':'PAA0028',
    'walltime':72,
    'n_nodes':1,
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
    'molar_fraction_total':{'Au':24, 'Pt':8},
    'simulation_cells':{
        'fcc1':{
            'incar':os.path.join('resource','INCAR_1'),
            'kpoints':os.path.join('resource','KPOINTS_1'),
            'poscar':os.path.join('resource','POSCAR_1'),
            'potcar':os.path.join('resource','POTCAR')
        },
        'fcc2':{
            'incar':os.path.join('resource','INCAR_2'),
            'kpoints':os.path.join('resource','KPOINTS_2'),
            'poscar':os.path.join('resource','POSCAR_2'),
            'potcar':os.path.join('resource','POTCAR')
        }
    }
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
