import pytest
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import Pymatmc2Data

resources_path = os.path.join('resources')
pymatmc2_config_path = os.path.join('pymatmc2.config')
pymatmc2_data_path = os.path.join('pymatmc2.results')
symbol_colors = {'Hf': 'blue', 'Zr': 'green', 'Ta': 'red', 'Nb': 'cyan'}

if __name__ == "__main__":

    configuration = Pymatmc2Configuration()
    configuration.read(path=pymatmc2_config_path)

    data = Pymatmc2Data()
    data.configuration = configuration
    data.read(path=pymatmc2_data_path)

    
    fig = plt.figure()
    ax = {}
    ax['bcc1'] = fig.add_subplot(411)
    ax['bcc2'] = fig.add_subplot(412)
    ax['bcc3'] = fig.add_subplot(413)
    ax['bcc4'] = fig.add_subplot(414)

    x = data.df['iteration']
    n_cells = len(configuration.cell_names)
    n_symbols = len(configuration.symbols)
    for i_cell, cn in enumerate(configuration.cell_names):
        y = 0
        for s in configuration.symbols:
            col_name_c = 'final.{}.{}'.format(cn,s)
            y += data.df[col_name_c]
            ax[cn].plot(x, y, linewidth=1., color = symbol_colors[s])

    plt.show()
            