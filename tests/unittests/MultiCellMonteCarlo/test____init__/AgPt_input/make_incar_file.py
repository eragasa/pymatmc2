import os
import subprocess
from mexm.io.vasp import Incar
from mexm.io.vasp import Kpoints
from mexm.io.vasp import Poscar

from pymatmc2 import Pymatmc2Configuration
from pymatmc2.multicell import create_random_cell

def make_pymatmc2_configuration_file(configuration_script_path: str):
    cmd = 'python {}'.format(configuration_script_path)
    subprocess.Popen(cmd)

def make_incar_files(src_incar_path, phase_paths):
    o_incar = Incar()
    o_incar.read(path=src_incar_path)

    for phase in phase_paths:
        dst_incar_path = os.path.join(phase,'INCAR')
        print('writing {}'.format(dst_incar_path))
        o_incar.write(path=dst_incar_path)

def make_kpoint_files(src_kpoint_path, phase_paths):
    o_kpoints = Kpoints()
    o_kpoints.read(path=src_kpoint_path)

    for phase in phase_paths:
        dst_kpoints_path = os.path.join(phase, 'KPOINTS')
        print('writing {}'.format(dst_kpoints_path))
        o_kpoints.write(path=dst_kpoints_path)

def make_poscar_files(src_configuration_path, phase_paths):
    o_configuration = Pymatmc2Configuration()
    o_configuration.read(path=src_configuration_path)

    total_concentration = o_configuration.concentration
    print(total_concentration)
    for phase in phase_paths:
        cell_type = 'fcc'
        cell = create_random_cell(
            cell_type=cell_type,
            composition=total_concentration,
            supercell=[2,2,2]
        )

        dst_poscar_path = os.path.join(phase, 'POSCAR')
        o_poscar = Poscar.initialize_from_mexm(cell)
        o_poscar.write(path=dst_poscar_path)
        
if __name__ == "__main__":
    # configuration_script_path = 'make_configuration.py'
    src_incar_path = 'INCAR'
    src_kpoint_path = 'KPOINTS'
    src_configuration_path = 'pymatmc2.config'
    
    phase_paths = ['fcc1', 'fcc2']

    # make_pymatmc2_configuration_file(configuration_script_path=configuration_script_path)
    make_incar_files(src_incar_path=src_incar_path, phase_paths=phase_paths)
    make_kpoint_files(src_kpoint_path=src_kpoint_path, phase_paths=phase_paths)
    make_poscar_files(src_configuration_path=src_configuration_path, phase_paths=phase_paths)