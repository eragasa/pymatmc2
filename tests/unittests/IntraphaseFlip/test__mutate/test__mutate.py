import pytest
import os
import shutil
from numpy import linalg
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell
from pymatmc2.multicellmutate import IntraphaseFlip

if __name__ == "__main__":
    configuration_path = os.path.join('resources','pymatmc2.config')
    configuration = Pymatmc2Configuration()
    configuration.read(path=configuration_path)

    multicell_initial_path = os.path.join('simulations','00000')
    multicell_initial = MultiCell()
    multicell_initial.configuration = configuration
    multicell_initial.read(path=multicell_initial_path)

    mutator = IntraphaseFlip()
    mutator.mutate_multicell(multicell=multicell_initial)

    assert isinstance(mutator.multicell_initial, MultiCell)
    print('initial_cell:')
    print(mutator.multicell_initial.concentration)
    print(mutator.multicell_initial.cell_concentration)
    try:
        print(mutator.multicell_initial.phase_molar_fraction)
    except linalg.LinAlgError as e:
        print('phase_molar_fraction:{}'.format(e))
    assert isinstance(mutator.multicell_candidate, MultiCell)
    print('candidate_cell:')
    print(mutator.multicell_candidate.concentration)
    print(mutator.multicell_candidate.cell_concentration)
    print(mutator.multicell_candidate.phase_molar_fraction)

    print(mutator.multicell_initial.src_path)
    print(mutator.multicell_candidate.src_path)
    next_iteration = 1

    results_path = os.path.join('results')
    results_iteration_path = os.path.join(
        results_path,'{:05}'.format(next_iteration)
    )
    dst_initial_path = os.path.join(
        results_iteration_path, 'initial'
    )
    dst_candidate_path = os.path.join(
        results_iteration_path, 'candidate'
    )


    # this block is just testing code
    if os.path.isdir(results_path):
        shutil.rmtree(results_path)
    os.mkdir(results_path)

    # this block is to archive previous simulations
    os.mkdir(results_iteration_path)
    print('archiving_initial_multicell:{}'.format(dst_initial_path))
    src_initial_path = mutator.multicell_initial.src_path
    print('initial_src_path:{}'.format(src_initial_path))

    os.mkdir(dst_initial_path)

    initial_files_to_archive = [
        'POSCAR', 'INCAR', 'KPOINTS', 'OSZICAR', 'CONTCAR', 'OUTCAR'
    ]
    for cn in mutator.multicell_initial.cell_names:
        os.mkdir(os.path.join(dst_initial_path,cn))
        for fn in initial_files_to_archive:
            shutil.copy(
                src=os.path.join(src_initial_path, cn, fn),
                dst=os.path.join(dst_initial_path, cn, fn)
            )
    
    os.mkdir(dst_candidate_path)
    src_candidate_path = mutator.multicell_initial.src_path
    print('src_candidate_path:{}'.format(src_candidate_path))
    print('archiving_candidate_multicell:{}'.format(dst_candidate_path))
    mutator.multicell_candidate.write(path=dst_candidate_path)
    for cn in mutator.multicell_candidate.cell_names:
        os.remove(os.path.join(dst_candidate_path, cn, 'POTCAR'))

    # write next simulation
    dst_candidate_path = os.path.join(
        os.path.join('simulations','{:05}'.format(next_iteration))
    )
    mutator.multicell_candidate.write(path=dst_candidate_path)
    exit()
 
