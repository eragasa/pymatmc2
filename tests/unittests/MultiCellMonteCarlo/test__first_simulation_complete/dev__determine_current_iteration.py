import os
from mexm.simulation import VaspSimulation
from pymatmc2 import MultiCellMonteCarlo

if __name__ == "__main__":
    kwargs_mc2 = {
        'configuration_path':os.path.join('resource','pymatmc2.config'),
        'results_path':'pymatmc2.results',
        'logfile_path':'pymatmc2.log',
        'simulations_path':os.path.join('simulations'),
        'is_restart':True
    }
    o_mc2 = MultiCellMonteCarlo(**kwargs_mc2)
    assert isinstance(o_mc2, MultiCellMonteCarlo)
    i_iteration, status = o_mc2.determine_current_iteration()
    o_mc2.create_next_structures(i_iteration=i_iteration)
    assert os.path.isdir('simulations')