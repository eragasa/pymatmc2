import os
from mexm.io.vasp import Outcar
from mexm.simulation import VaspSimulation
from pymatmc2 import MultiCell

multicell_dict = {
    'fcc1':os.path.join('simulations','fcc1_000_T400_P0'),
    'fcc2':os.path.join('simulations','fcc2_000_T400_P0')
}

simulations = {}
for k, v in multicell_dict.items():
    cell_name = k
    simulation_path = v
    simulations[cell_name] = VaspSimulation()
    simulations[cell_name].read(simulation_path=simulation_path)

for k, v in simulations.items():
    assert isinstance(v.outcar, Outcar)
    toten = v.outcar.total_energy
    n_atoms = v.poscar.n_atoms
    energy_per_atom = toten/n_atoms
    print(energy_per_atom)

