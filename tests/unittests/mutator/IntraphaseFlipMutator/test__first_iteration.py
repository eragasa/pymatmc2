import os
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell

from pymatmc2.mutator import IntraphaseFlipMutator

pymatmc2_configuration_AuPt_
def get_configuration() -> Pymatmc2Configuration:
    configuration_path = 'pymatmc2.config'
    configuration = Pymatmc2Configuration()
    configuration.read(configuration_path)
    return configuration

def get_mutator() -> IntraphaseFlipMutator:
    configuration = get_configuration()
    mutator = IntraphaseFlipMutator()
    mutator.configuration = configuration
    return mutator

def get_iteration_path(configuration: Pymatmc2Configuration, i_iteration: int) -> str:
    path = os.path.join(
        configuration.results_path, 
        configuration.phase_point_string, 
        configuration.get_iteration_string(i=i_iteration)
    )
    return path

def test__first_iteration():
    mutator = get_mutator()


def dev__first_iteration():
    mutator = get_mutator()

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

