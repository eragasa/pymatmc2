import math
import numpy as np
from typing import Dict, List
from mexm.elements import ELEMENTS
from mexm.structure import SimulationCell
from mexm.structure import make_super_cell

def get_bcc_cell(r: float) -> SimulationCell:
    #TODO:

    a0 = 4 / math.sqrt(3) * r
    simulation_cell = SimulationCell()
    simulation_cell.a0 = a0
    simulation_cell.add_atom('Ni',[0.0, 0.0, 0.0])
    simulation_cell.add_atom('Ni',[0.5, 0.5, 0.5])
    return simulation_cell

def get_fcc_cell(r: float) -> SimulationCell:
    #TODO:
    
    a0 = math.sqrt(8) * r
    simulation_cell = SimulationCell()
    simulation_cell.a0 = a0
    simulation_cell.add_atom('Ni',[0.0, 0.0, 0.0])
    simulation_cell.add_atom('Ni',[0.5, 0.5, 0.0])
    simulation_cell.add_atom('Ni',[0.5, 0.0, 0.5])
    simulation_cell.add_atom('Ni',[0.0, 0.5, 0.5])
    return simulation_cell

def get_hcp_cell():
    #TODO:
    simulation_cell = SimulationCell()
    return simulation_cell

def create_random_cell(
    cell_type: str,
    composition: Dict[str, float],
    supercell: List[int]
) -> SimulationCell:
    """
    Arguments:
        cell_type (str): the type of base cell structure to use
        composition (List[str, float]): the composition of the alloy
        supercell (List[int]): the supercell 
    """

    assert isinstance(cell_type, str)
    assert isinstance(composition, dict)
    assert isinstance(supercell, list)
    assert all([isinstance(k, int) for k in supercell])

    sum_composition = sum(composition.values())
    composition_ = {
        k:v/sum_composition for k, v in composition.items()
    }
    max_atmrad = max([ELEMENTS[k].atmrad for k in composition_])
    cell_type_to_unit_cell_map = {
        'bcc': get_bcc_cell,
        'fcc': get_fcc_cell,
        'hcp': get_hcp_cell
    } 
    cell = cell_type_to_unit_cell_map[cell_type](r=max_atmrad)
    cell = make_super_cell(structure=cell, sc=supercell)

    idx_atoms_all = [k for k in range(cell.n_atoms)]
    idx_atoms = {}

    atoms_processed = []
    for symbol in composition_:
        idx_atoms[symbol] = []
        
        sum_composition = sum([
            v for k, v in composition_.items() if k not in atoms_processed
        ])
        probability = composition_[symbol]/sum_composition

        n_atoms = int(cell.n_atoms * probability) \
            - sum([len(v) for v in idx_atoms.values()])

        for i_atom in range(n_atoms):
            idx = np.random.choice(idx_atoms_all, 1).tolist()[0]
            idx_atoms[symbol].append(idx)
            idx_atoms_all.remove(idx)
            atoms_processed.append(symbol)
    
    for symbol in idx_atoms:
        for idx_atom in idx_atoms[symbol]:
            cell.atomic_basis[idx_atom].symbol = symbol

    return cell
