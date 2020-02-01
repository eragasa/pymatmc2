from collections import OrderedDict
from mexm.io.vasp import Poscar
from pymatmc2.simulationcells import SimulationCells

if __name__ == "__main__":
    sim_cells = SimulationCells(n_cells=2)
    sim_cells.cells[0].read('POSCAR_1')
    sim_cells.cells[1].read('POSCAR_2')
    print(sim_cells.cells[0].symbols)
    print(sim_cells.cells[1].symbols)
    print(sim_cells.symbols)

    for symbol in sim_cells.symbols:
        print(symbol, sim_cells.cells[0].get_number_of_atoms(symbol))
        print(symbol, sim_cells.cells[1].get_number_of_atoms(symbol))
        print(symbol, sim_cells.get_number_of_atoms(symbol))
    print('all', sim_cells.cells[0].get_number_of_atoms())
    print('all', sim_cells.cells[1].get_number_of_atoms())
    print('all', sim_cells.get_number_of_atoms())
