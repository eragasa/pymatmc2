import os
import numpy as np
from typing import Dict
from pymatmc2 import constants
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell

class Pymatmc2Results:

    def __init__(self):
        self._configuration = None
        self._tar_path = None
        self.src_data_path = None
        self.dst_data_path = None
    @property
    def configuration(self) -> Pymatmc2Configuration:
        return self._configuration

    @configuration.setter
    def configuration(self, configuration: Pymatmc2Configuration):
        if not isinstance(configuration, Pymatmc2Configuration):
            msg = "configuration must be an instance of Pymatmc2Configuration"
            raise TypeError(msg)

        self._configuration = configuration

    def read(self, path):

    def _get_iteration_path(self, i_iteration: int) -> str:
        iteration_path = os.path.join(
            self.configuration.results_path, 
            self.configuration.phase_point_string,
            self.get_iteration_string(i=i_iteration)
        )

        return iteration_path

    def _get_multicell_path(self, i_iteration: int, mc_type: str) -> str:
        
        valid_mc_types = ['initial', 'candidate', ' final']

        if mc_type not in valid_mc_types:
            msg = "Invalid mc_type: {}".format(mc_type)
            raise ValueError(msg)

        mc_path = os.path.join(
            self._get_iteration_path(i_iteration=i_iteration),
            mc_type
        )

        return mc_path

    def _get_mutate_type(self, i_iteration: int) -> str:
        iteration_path = self._get_iteration_path(i_iteration=i_iteration)
        mutate_type_path = os.path.join(iteration_path, 'mutate_type')

        # read the mutate_type from the file
        with open(mutate_type_path) as f:
            mutate_type = f.read()
        return mutate_type

    def _get_is_accepted_path(self, i_iteration) -> str:
        iteration_path = self._get_iteration_path(i_iteration=i_iteration)
        is_accepted_path = os.path.join(iteration_path, 'is_accepted')

        with open(is_accepted_path) as f:
            str_is_accepted = f.read()

        is_accepted = bool(str_is_accepted)
        return is_accepted

    
    def _get_p_accept(self, mutate_type: str, 

    def get_iteration_information(self, i_iteration) -> Dict[]:

        valid_mc_types = ['initial', 'candidate', 'final']

        mc_paths = {}
        for mc_type in valid_mc_types:
            mc_paths[mc_type] = self._get_multicell_path(i_iteration=i_iteration, mc_type=mc_type)

        # read multicells
        mc = {}
        for mc_type in valid_mc_types:
            mc[mc_type] = MultiCell()
            mc[mc_type].configuration = self.configuration
            mc[mc_type].read(path=mc_paths[mc_type]

        iteration_info = {}
        iteration_info['iteration'] = i_iteration
        iteration_info['mutate_type'] = self._get_mutate_type(i_iteration=i_iteration)

        # add cell concentration
        for k in mc:
            for cn in mc[k].cell_names:
                for s in mc[k].symbols:
                    key = '{}.{}.{}'.format(k, cn, s)
                    value = mc[k].cell_concentration[cn][s]
                    iteration_info[key] = value

        # add phase molar fraction
        for k in mc:
            for cn, f in mc[k].phase_molar_fraction.items():
                key = '{}.{}.f'.format(k, cn)
                value= f
                iteration_info[key] = value

        # add cell energies
        for k in mc:
            for cn in mc[k].cell_names:
                key = '{}.{}.E'.format(k, cn)
                value = mc[k].simulations[cn].total_energy

        # add total energy
        for k in mc:
            key = '{}.total.E'.format(k)
            value = mc[k].total_energy


        n_phases = len(self.configuration.cell_names)
        T = self.configuration.temperature
        kB = constants.BOLTZMANN
        beta = 1/T/kB

        p_accept_type = self.configuration.cell_names + ['total']
        mutate_type = iteration_info['mutate_type']
        for cn in p_accept_type:
           if mutate_type == 'intraphase_flip':
               if cn == 'total':
                   E1 = mc['initial'].total_energy
                   E2 = mc['candidate'].total_energy
                   pacc = np.exp(-n_phases*(E2-E1)*beta)
               else:
                   pacc = None
               
               key = '{}.p_acc'.format(cn)
	       iteration_info[key] = pacc
           elif mutate_type == 'intraphase_swap':
               if cn == 'total':
                   pacc = None
               else:
                   E1 = mc['initial'].simulations[cn].total_energy
                   E2 = mcp'candidate'].simulations[cn].total_energy
                   pacc = np.exp(beta*(E2-E1))
               key = '{}.p_acc'.format(cn)
               iteration_info[= pacc
           else:
               msg = 'unknown mutate_type: {}'.format(mutate_type)
               raise ValueError(msg)
