
results_path = 'Results'
rejected_path = 'Rejected'
temperature = 400

cells = []
with open(results_path, 'r') as f:
    lines = f.readlines()
for i_line, line in enumerate(lines):
    if i_line == 0:
        pass
    else:
        tokens = line.strip().split()
        cells.append([tokens[1], tokens[2], tokens[3]])
        cells.append([tokens[4], tokens[5], tokens[6]])
with open(rejected_path, 'r') as f:
    lines = f.readlines()
for i_line, line in enumerate(lines):
    if i_line == 0:
        pass
    else:
        tokens = line.strip().split()
        cells.append([tokens[1], tokens[2], tokens[3]])
        cells.append([tokens[4], tokens[5], tokens[6]])

cell_energies = {}
for cell in cells:
    n_symbol_1 = cell[0]
    n_symbol_2 = cell[1]
    cell_energy = cell[2]
    composition_key = '{}_{}'.format(n_symbol_1, n_symbol_2)
    if composition_key not in cell_energies:
        cell_energies[composition_key] = {}
        cell_energies[composition_key]['composition'] = [int(n_symbol_1), int(n_symbol_2)]
        cell_energies[composition_key]['energies'] = [float(cell_energy)]
    else:
        cell_energies[composition_key]['energies'].append(float(cell_energy))

for k, v in cell_energies.items():
    cell_energies[k]['energies'] = list(set(cell_energies[k]['energies']))

import numpy as np
from pymatmc2 import constants
for k, v in cell_energies.items():
    Q = 0
    kB = constants.BOLTZMANN
    T = temperature
    cell_energies[k]['energy'] = min(cell_energies[k]['energies'])

import matplotlib.pyplot as plt
x = [v['composition'][0] for k, v in cell_energies.items()]
y = [v['energy'] for k, v in cell_energies.items()]
zipped_xy = list(zip(x,y))
zipped_xy = sorted(zipped_xy, key = lambda x: x[0])
print(zipped_xy)
for i, xy in enumerate(zipped_xy):
    if i == 0:
        dU = zipped_xy[i][1] - zipped_xy[i+1][1]
    elif i == len(zipped_xy)-1:
        dU = zipped_xy[i-1][1] - zipped_xy[i][1]
    else:
        du1 = zipped_xy[i][1] - zipped_xy[i+1][1]
        du2 = zipped_xy[i-1][1] - zipped_xy[i][1]
        dU = (du1 + du2)/2
    x = xy[0]
    plt.plot(x, dU, marker='x', color='black')
plt.show()