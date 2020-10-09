import numpy as np
from collections import OrderedDict
from pymatmc2 import constants
class ConfigurationalEntropy():
    def __init__(self):
        self.N = []

if __name__ == "__main__":
    N_total = 32
    N_phases = OrderedDict()
    N_phases['fcc1'] = {'s1':None, 's2':None}
    N_phases['fcc2'] = {'s1':None, 's2':None}
    N_phases['fcc1']['s1'] = np.linspace(0, N_total, N_total + 1)
    N_phases['fcc1']['s2'] = N_total * np.ones(N_total + 1) - N_phases['fcc1']['s1'] 
    N_phases['fcc2']['s2'] = np.linspace(0, N_total, N_total + 1)
    N_phases['fcc2']['s1'] = N_total * np.ones(N_total + 1) - N_phases['fcc2']['s2']

    # calculate entropy
    kB = constants.BOLTZMANN
    phase_names = ['fcc1', 'fcc2']
    symbols = ['s1', 's2']
    temperatures = [0, 100, 200, 300 ,400]
    S_mix = OrderedDict()
    for phase_name in phase_names:
        sum = 0
        for symbol in symbols:
            N_a = N_phases[phase_name][symbol]
            sum += N_a * np.log(N_a/N_total)
        S_mix[phase_name] = -kB * sum
    
    # calculate derivative of entropy
    dS_mix = OrderedDict()
    for phase_name in phase_names:
        sum = 0
        for symbol in symbols:
            N_a = N_phases[phase_name][symbol]
            sum += 1 + np.log(N_a/N_total)
        dS_mix[phase_name] = -kB * sum

    import matplotlib.pyplot as plt
    for T in temperatures:
        x = N_phases[phase_name]['s1'] / N_total
        y = dS_mix[phase_name] * T
        plt.plot(x,y, label='{}K'.format(T))
    plt.xlabel('Composition Percentage')
    plt.ylabel('TdS [eV]')
    plt.legend()
    plt.show()

