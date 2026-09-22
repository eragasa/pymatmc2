from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

from mexm.crystal import SimulationCell

def convert_structure_mexm_to_phonopy(simulation_cell):
    assert isinstance(simulation_cell, SimulationCell)

    kwargs = {
        symbols = [atom.symbol for atom in simulation_cell.atomic_basis],
        cell= (simulation_cell.H * simulation_cell.a0),
        scaled_positions=[atom.positions for atom in simulation_cell.atomic_basis]
    }
    phonopy_cell = PhonopyAtoms(**kwargs)

def convert_structure_phonopy_to_mexm(simulation_cell):
    assert isinstance(simulation_cell, PhonopyAtoms)

    kwargs = {
        symbols = []
    }
