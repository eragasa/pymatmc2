import pytest
import os
import shutil
from time import sleep
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCellMonteCarlo

def test__constructor__default__AgPt():
    kwargs_mc2 = {
        'configuration_path':os.path.join('AgPt_input','pymatmc2.config'),
        'results_path':'results',
        'logfile_path':'pymatmc2.log',
        'simulations_path':'simulations',
        'is_restart':False
    }

    o_mc2 = MultiCellMonteCarlo(**kwargs_mc2)

    simulations_path = 'simulations'
    simulations_phase_point_path = os.path.join('simulations', '400K_0GPa')
    assert os.path.isdir(simulations_path)
    assert os.path.isdir(simulations_phase_point_path)

    # cleanup
    shutil.rmtree('simulations')
    shutil.rmtree('results')

def test__constructor__default__HfZrTaNb():
    kwargs_mc2 = {
        'configuration_path':os.path.join('HfZrTaNb_input','pymatmc2.config'),
        'results_path':'results',
        'logfile_path':'pymatmc2.log',
        'simulations_path':'simulations',
        'is_restart':False
    }

    o_mc2 = MultiCellMonteCarlo(**kwargs_mc2)

    simulations_path = 'simulations'
    simulations_phase_point_path = os.path.join('simulations', '800K_0GPa')
    assert os.path.isdir(simulations_path)
    assert os.path.isdir(simulations_phase_point_path)

    #cleanup
    shutil.rmtree('simulations')
    shutil.rmtree('results')

