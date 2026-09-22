import os
import platform
from copy import deepcopy
from pymatmc2 import Pymatmc2Configuration

# generic hpc configuration for OSC
configuration_osc_generic = {}
configuration_osc_generic['type'] = 'torque'
configuration_osc_generic['configuration'] = {
    'account':'PAA0028',
    'walltime':1,
    'n_nodes':2,
    'ppn':40,
    'jobname':'Si_job_name',
    'errpath':'job.err',
    'stdpath':'job.out',
    'modules':['intel/19.0.5', 'intelmpi/2019.3'],
    'cmd':'mpiexec $VASP_STD_BIN > vasp.out'
}

# specific hpc cluster configurations, only differs by the number
# of processing units per node
configuration_pitzer = deepcopy(configuration_osc_generic)
configuration_ptizer['configuration']['n_nodes'] = 2
configuration_pitzer['configuration']['ppn'] = 40
configuration_owens = deepcopy(configuration_osc_generic)
configuration_ptizer['configuration']['n_nodes'] = 2
configuration_owens['configuration']['ppn'] = 28
configuration_ruby = deepcopy(configuration_osc_generic)
configuration_ruby['configuration']['n_nodes'] = 2
configuration_ruby['configuraton']['ppn'] = 20

hpc_configurations = {
    'pitzer': configuration_pitzer,
    'owens': configuration_owens,
    'ruby': configuration_pitzer
}

login_nodenames = {}
login_nodenames['ruby'] = [
    'ruby01.osc.edu',
    'ruby02.osc.edu'
]
login_nodenames['owens'] = [
    'owens-login01.hpc.osc.edu', 
    'owens-login02.hpc.osc.edu', 
    'owens-login03.hpc.osc.edu',
    'owens-login04.hpc.osc.edu'
]
login_nodenames['pitzer'] = [
    'pitzer-login01.hpc.osc.edu', 
    'pitzer-login02.hpc.osc.edu', 
    'pitzer-login03.hpc.osc.edu',
    'pitzer-login04.hpc.osc.edu'   
]

def get_hpc_cluster_name():
    hostname = platform.node()
    for k, v in login_nodenames.items():
        if hostname in v:
            return k
    
configuration = {}

# choosing hpc manager
hpc_cluster_name = get_hpc_cluster_name()
configuration['hpc_manager'] = hpc_configurations[hpc_cluster_name]

# configuration of results output
configuration['results'] = {}
configuration['results']['dir'] = 'results'
configuration['results']['tar_path'] = 'pymatmc2.results.tar'
configuration['results']['file_path'] = 'pymatmc2.results'

# run_scheme options: mpi_serial, mpi_crontab, mpi_spawn
configuration['calculator'] = {
    'calculator':'vasp',
    'simulation_type':'vasp_min_all',
    'run_scheme': 'mpi_crontab'
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
