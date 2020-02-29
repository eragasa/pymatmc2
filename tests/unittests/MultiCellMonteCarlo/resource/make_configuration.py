import os
from pymatmc2 import Pymatmc2Configuration

configuration = {}
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
            'incar':os.path.join('resources','INCAR_1'),
            'kpoints':os.path.join('resources','KPOINTS_1'),
            'poscar':os.path.join('resources','POSCAR_1'),
            'potcar':os.path.join('resources','POTCAR')
        },
        'fcc2':{
            'incar':os.path.join('resources','INCAR_2'),
            'kpoints':os.path.join('resources','KPOINTS_2'),
            'poscar':os.path.join('resources','POSCAR_2'),
            'potcar':os.path.join('resources','POTCAR')
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
