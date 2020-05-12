import pytest
import os
import shutil
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell

from pymatmc2.mutator import IntraphaseFlipMutator

src_configuration_path = os.path.join('resources','test__first_iteration', 'pymatmc2.config')
src_results_path = os.path.join('resources', 'test__first_iteration', 'results')
dst_configuration_path = 'pymatmc2.config'
dst_results_path = os.path.join('results')

def get_iteration_path(configuration: Pymatmc2Configuration, i_iteration: int) -> str:
    path = os.path.join(
        configuration.results_path, 
        configuration.phase_point_string, 
        configuration.get_iteration_string(i=i_iteration)
    )
    return path

def setup():
    if not os.path.isdir(dst_results_path):
        shutil.copytree(src=src_results_path, dst=dst_results_path)
        while not os.path.isdir(dst_results_path):
            pass
    if not os.path.isfile(dst_configuration_path):
        shutil.copy(src=src_configuration_path, dst=dst_configuration_path)
        while not os.path.isfile(dst_configuration_path):
            pass

def cleanup():
    shutil.rmtree(dst_results_path)
    while os.path.isdir(dst_results_path):
        pass

    os.remove(dst_configuration_path)

def test__first_iteration():
    configuration_path = dst_configuration_path
    configuration = Pymatmc2Configuration()
    configuration.read(path=configuration_path)

    mutator = IntraphaseFlipMutator()
    mutator.configuration = configuration

    iteration = 1 
    
    mc_types = ['initial', 'candidate', 'final']
    mc = {}
    for mc_type in mc_types:
        mc_path = os.path.join(
            get_iteration_path(mutator.configuration, iteration),
            mc_type
        )
        mc[mc_type] = MultiCell()
        mc[mc_type].configuration = mutator.configuration
        mc[mc_type].read(mc_path)        
    
    mutator.accept_or_reject(
        multicell_initial=mc['initial'], 
        multicell_candidate=mc['candidate'],
        temperature = mutator.configuration.temperature,
        pressure = mutator.configuration.pressure
    )
    cleanup()

def dev__first_iteration():
    configuration = Pymatmc2Configuration()
    configuraiton.read(path=configuration_path)

    mutator = IntraphaseFlipMutator()
    mutator.configuration = configuration

    iteration = 1 

    mc_types = ['initial', 'candidate', 'final']
    mc = {}
    for mc_type in mc_types:
        mc_path = os.path.join(
            get_iteration_path(mutator.configuration, iteration),
            mc_type
        )
        mc[mc_type] = MultiCell()
        mc[mc_type].configuration = mutator.configuration
        mc[mc_type].read(mc_path)        
    
    mutator.accept_or_reject(
        multicell_initial=mc['initial'], 
        multicell_candidate=mc['candidate'],
        temperature = mutator.configuration.temperature,
        pressure = mutator.configuration.pressure
    )

    print(mutator.multicell_initial.cell_concentration)
    print(mutator.multicell_candidate.cell_concentration)
    print(mutator.is_accept)
    print(mutator.p_accept)
    print(mutator.p_accept_rnd)
if __name__ == "__main__":
    dev__first_iteration()

