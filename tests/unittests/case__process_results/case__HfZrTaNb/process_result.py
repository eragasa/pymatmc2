""" Processing the results from pymatmc2 """

import os
import numpy as np

from typing import List

from pymatmc2 import constants
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell
from pymatmc2 import Pymatmc2Data
from pymatmc2 import Pymatmc2Results

results_path = 'results'
results_tar_path = 'pymatmc2.results.tar'
results_file_path = 'pymatmc2.results'

configuration_file_path = 'pymatmc2.config'
configuration = Pymatmc2Configuration()
configuration.read(path=configuration_file_path)

print(configuration.results_path)
print(configuration.results_tar_path)
print(configuration.results_file_path)

o_results = Pymatmc2Results()
o_results.configuration = configuration

o_data = Pymatmc2Data()
o_data.configuration = configuration
i_iteration = 0


o_data.write_header()
while True:
    i_iteration += 1
    iteration_path = o_results.get_iteration_path(i_iteration)
    print('processing iteration {}: {}'.format(i_iteration, iteration_path))
    if not os.path.isdir(iteration_path):
        print('completed iteration {}'.format(i_iteration))
        break
  
    symbols = list(configuration.symbols)
    n_symbols = len(symbols) 

    #mutate_type = o_results.get_mutate_type(i_iteration)
    #mc_paths = o_results.get_multicell_paths(i_iteration) 
    #mc = o_results.get_multicells(i_iteration)
    iteration_info = o_results.get_iteration_info(i_iteration)
    o_data.write_results(iteration_info=iteration_info)


