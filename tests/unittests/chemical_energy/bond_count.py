import numpy as np
from mexm.structure import SimulationCell
from mexm.structure import make_super_cell

def get_distance(cell, index_1, index_2, is_debug=False, max_distance=1e6):
    x1 = np.array(cell.atomic_basis[index_1].position) * cell.lattice.a0
    x2 = np.array(cell.atomic_basis[index_2].position) * cell.lattice.a0

    d1 = max_distance
    for ia1 in [-1,0,1]:
        for ia2 in [-1,0,1]:
            for ia3 in [-1,0,1]:
                shift = ia1 * supercell.a1 + ia2 * supercell.a2 + ia3 * supercell.a3
                x2_ = x2 + shift
                d1_ = np.dot(x1-x2_, x1-x2_)
                d1 = min(d1, d1_)
    return d1 

cell = SimulationCell()
cell.lattice.a0 = 3.50579800
cell.add_atom(symbol='Ni', position=[0.0, 0.0, 0.0], atom_id="Ni1")
cell.add_atom(symbol='Ni', position=[0.0, 0.5, 0.5], atom_id="Ni2")
cell.add_atom(symbol='Ni', position=[0.5, 0.0, 0.5], atom_id="Ni3")
cell.add_atom(symbol='Ni', position=[0.5, 0.5, 0.0], atom_id="Ni4")
supercell = make_super_cell(cell, sc = [3,3,3])

bond_lengths = []
for i in range(supercell.n_atoms):
    for j in range(supercell.n_atoms):
        if i < j:
            d = get_distance(supercell, i, j)
            bond_lengths.append(d)
bond_lengths = np.array(bond_lengths)

from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

x = np.linspace(0, max(bond_lengths)*1.25, 1000)
kde = gaussian_kde(bond_lengths, bw_method=0.05)
bond_length_kde = kde(x)
peaks_idx, peak_properties = find_peaks(kde(x), height=0)

import matplotlib.pyplot as plt

fix, ax = plt.subplots()
ax.hist(bond_lengths, bins=100, density=True)
ax.plot(x, kde(x), color='black')
ax.plot(x[peaks_idx], bond_length_kde[peaks_idx], 'x', color='black')
ax.set_xlim(0, max(bond_lengths)*1.25)
ax.set_xlabel('r')
ax.set_ylabel('g(r)')
plt.show()