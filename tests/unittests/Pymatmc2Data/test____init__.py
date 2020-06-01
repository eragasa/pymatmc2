import pytest
import os
import pandas as pd
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import Pymatmc2Data

resources_path = os.path.join('resources', 'test____init__')
pymatmc2_config_path = os.path.join(resources_path, 'pymatmc2.config')
pymatmc2_data_path = os.path.join(resources_path, 'pymatmc2.results')

def test__default_constructor():
    data = Pymatmc2Data()

def test__property__configuration__set():
    configuration = Pymatmc2Configuration()
    configuration.read(path=pymatmc2_config_path)

    data = Pymatmc2Data()
    data.configuration = configuration

    assert isinstance(data.configuration, Pymatmc2Configuration)

def test__read():
    configuration = Pymatmc2Configuration()
    configuration.read(path=pymatmc2_config_path)

    data = Pymatmc2Data()
    data.configuration = configuration
    data.read(path=pymatmc2_data_path)

    assert isinstance(data.df, pd.DataFrame)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    configuration = Pymatmc2Configuration()
    configuration.read(path=pymatmc2_config_path)

    data = Pymatmc2Data()
    data.configuration = configuration
    data.read(path=pymatmc2_data_path)

    print(data.df)

    import numpy as np
    x = data.df['iteration']
    #y = np.zeros(x.shape)
    #for cn in configuration.cell_names:
    #    col_name = 'final.f.{}'.format(cn)
    #    y = y + data.df[col_name]
    #    plt.plot(x,y)
    #    n_symbols = len(configuration.symbols)

    from collections import OrderedDict
    y_phases = OrderedDict()
    for cn in configuration.cell_names:
        col_name = 'final.f.{}'
        y_phases[cn] = data.df[col_name]

    y_compositions = OrderedDict()
    for cn in configuration.cell_names:
        y_compositions[cn] = OrderedDict()
    for cn in configuration.cell_names:
         for s in configuration.symbols:
             col_name = 'final.{}.{}'.format(cn, s)
             y_compositions[cn][s] = data.df[col_name]

    fig = plt.figure()
    ax = []
    ax.append(fig.add_subplot(333))
    
    y_phase = np.zeros(x.shape)
    xlims = [min(x), max(x)]
    ylims = [0,1]
    for cn in configuration.cell_names:
        col_name_f = 'final.f.{}'.format(cn)
        for s in configuration.symbols:
            col_name_c = 'final.{}.{}'.format(cn,s)
            y += data.df[col_name_f] * data.df[col_name_c]
            plt.plot(x,y)
    plt.show()
            