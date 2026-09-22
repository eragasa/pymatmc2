import os

from mexm.io.vasp import Poscar

from pymatmc2 import Pymatmc2Configuration
from pymatmc2.widom import WidomTest


configuration_path = os.path.join('input', 'pymatmc2.config')
configuration = Pymatmc2Configuration()
configuration.read(path=configuration_path)

widom_n_iterations = 100

initial_cell_path = os.path.join('input', 'fcc1', 'POSCAR')

o_widom = WidomTest()
o_widom.configuration = configuration
o_widom.symbol1 = 'Au'
o_widom.symbol2 = 'Pt'
o_widom.phase = 'fcc1'
o_widom.cell_initial = Poscar()
o_widom.cell_initial.read(path=initial_cell_path)

assert isinstance(o_widom.symbol1, str)
assert isinstance(o_widom.symbol2, str)
assert isinstance(o_widom.cell_initial, Poscar)

assert o_widom.symbol1 in o_widom.cell_initial.symbols

o_widom.generate_structure_files(n_iterations=widom_n_iterations)