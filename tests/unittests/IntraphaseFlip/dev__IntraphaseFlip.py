import os
from collections import OrderedDict
from mexm.simulation import VaspSimulation
from pymatmc2 import MultiCell
from pymatmc2.multicellmutate import IntraphaseFlip

"""
Currently the testing example uses simulation from an intraphase FLIP, these
cases need to be rewritten for an intraphase SWAP example.  (e.g. mass balance 
isn't preserved in this example)
"""

def get_simulation_path() -> str:

    simulation_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '..', # unittests
        '..', # tests
        'resource',
        'AgPt',
        'simulations'
    )

    return os.path.abspath(simulation_path)

def get_multicell_iteration_0() -> MultiCell:
    multicell = MultiCell()
    multicell.simulations = OrderedDict()
    multicell.simulations['fcc1'] = VaspSimulation()
    multicell.simulations['fcc1'].read(
        os.path.join(get_simulation_path(), 'fcc1_000_T400_P0')
    )
    multicell.simulations['fcc2'] = VaspSimulation()
    multicell.simulations['fcc2'].read(
        os.path.join(get_simulation_path(), 'fcc2_000_T400_P0')
    )
    return multicell

def get_multicell_iteration_1() -> MultiCell:
    multicell = MultiCell()
    multicell.simulations = OrderedDict()
    multicell.simulations['fcc1'] = VaspSimulation()
    multicell.simulations['fcc1'].read(
        os.path.join(get_simulation_path(), 'fcc1_001_T400_P0')
    )
    multicell.simulations['fcc2'] = VaspSimulation()
    multicell.simulations['fcc2'].read(
        os.path.join(get_simulation_path(), 'fcc2_001_T400_P0')
    )
    return multicell

def dev__accept_or_reject():
    print(80*'-')
    print('IntraphaseFlip.accept_or_reject()')
    print(80*'-')
    simulation_path = get_simulation_path()
    multicell0 = get_multicell_iteration_0()
    multicell1 = get_multicell_iteration_1()
    temperature = 400
    print('simulation_path:{}'.format(simulation_path))
    obj = IntraphaseFlip()
    obj.accept_or_reject(
        multicell_initial = multicell0,
        multicell_candidate = multicell1,
        temperature = temperature
    )

def dev__mutate_multicell():
    print(80*'-')
    print('IntraphaseFlip.mutate_multicell()')
    print(80*'-')
    multicell0 = get_multicell_iteration_0()
    obj = IntraphaseFlip()
    multicell1 = obj.mutate_multicell(multicell = multicell0)
    
    assert isinstance(multicell0, MultiCell)
    assert isinstance(multicell1, MultiCell)

if __name__ == "__main__":
    dev__accept_or_reject()
    dev__mutate_multicell()
