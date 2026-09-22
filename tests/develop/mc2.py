"mc2 script"

from scipy.optimize import nnls
from pymatmc2 import MultiCellMonteCarlo
from pymatmc2.utils import clear_folders
from pymatmc2.utils import save_log
from pymatmc2.utils import get_ratio
from pymatmc2.utils import get_structure
from pymatmc2.utils import run_job
from pymatmc2.utils import stop_check
from pymatmc2.utils import prepare

if __name__ == "__main__":
    mc2 = MultiCellMonteCarlo(is_restart=False)
    mc2.run()
    exit()

    k_B = BOLTZMANN
    C = COULOMB
    T = configuration.temperature
    N_cell = configuration.n_cells
    # prepare()

    utils.clear_lock_file()
    while True:
        stop_check()
        folders = get_structure('flip_alt')
        if len(folders) == 0:
            r1 = Results()
            r1.step += 1
            r1.add_results(REJECTED_FILE)
            continue
        for folder in folders:
            run_job(folder)
        
        # why am initializing results here
        r1 = Results()
        r2 = Results()
        r2.read_next()
        save_log('{:>5d} {}\n'.format(r2.step, time.strftime("%Y-%m-%d %H:%M")))

        # what is this probability??
        probability = np.exp((r1.total_energy - r2.total_energy) * N_cell *
                             C / k_B / T)
        r2.probability = np.minimum(probability, 1.0)
        if np.random.rand() < r2.probability:
            r2.add_results(RESULTS_FILE)
            r2.tar_file()
            if r2.step >= 99999:
                save_log("Maximum step 99999 reached.")
                exit()
        else:
            r2.add_results(REJECTED_FILE)
