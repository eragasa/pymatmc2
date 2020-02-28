from mexm.simulation import VaspSimulation
from mexm import MultiCellMonteCarlo

if __name__ == "__main__":
    o_mc2 = MultiCellMonteCarlo()
    o_mc2.read_configuration(path='mc2.in')
    print('calculator_type:',o_mc2.calculator_type)
    print('vasp_std_bin:',o_mc2.vasp_std_bin)
    print('lammps_serial_bin:',o_mc2.lammps_serial_bin)
    print('lammps_mpi_bin:',o_mc2.lammps_mpi_bin)
    o_mc2.prepare_simulations()


    while True:
        o_mc2.check_for_stop_file()
        folders = get_structure('flip_alt')
        if len(folders) == 0:
            r1 = Results()
            r1.step += 1
            r1.add_results(REJECTED_FILE)
            continue
        for folder in folders:
            run_job(folder)
        r1, r2 = Results(), Results()
        r2.read_next()
        save_log('{:>5d} {}\n'.format(r2.step, time.strftime("%Y-%m-%d %H:%M")))
        probability = np.exp((r1.total_energy - r2.total_energy) * N_CELL *
                             COULOMB / BOLTZMANN / TEMPERATURE)
        r2.probability = np.minimum(probability, 1.0)
        if np.random.rand() < r2.probability:
            r2.add_results(RESULTS_FILE)
            r2.tar_file()
            if r2.step >= 99999:
                save_log("Maximum step 99999 reached.")
                exit()
        else:
            r2.add_results(REJECTED_FILE)
