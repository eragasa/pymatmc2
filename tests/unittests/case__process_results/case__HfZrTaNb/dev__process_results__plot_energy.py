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
    ax = []
    ax.append(fig.add_subplot(111))

    x = data.df['iteration']
    y = data.df['final.all.toten']
    ax[0].plot(x, y, linewidth=2., color='black')

    # set x limits, no y limits
    ax[0].set_xlim(min(x), max(x))
    plt.show()
            