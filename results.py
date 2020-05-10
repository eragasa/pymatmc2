import os
import numpy as np
from typing import Dict
from pymatmc2 import constants
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import MultiCell

class Pymatmc2Results():
    valid_mc_types = ['initial', 'candidate', 'final']
    def __init__(self):
        self._configuration = None

    @property
    def configuration(self) -> Pymatmc2Configuration:
        return self._configuration

    @configuration.setter
    def configuration(self, configuration: Pymatmc2Configuration):
        self._configuration = configuration

    def get_iteration_path(self, i_iteration:int) -> str:
        iteration_path = os.path.join(
            self.configuration.results_path,
            self.configuration.phase_point_string,
            self.configuration.get_iteration_string(i=i_iteration)
        )
        return iteration_path

    def get_mc_path(self, i_iteration: int, mc_type: str) -> str:

        if mc_type not in Pymatmc2Results.valid_mc_types:
            msg = 'Invalid mc_type: {}'
            msg = msg.format(mc_type)
            raise ValueError(msg)

        iteration_path = self.get_iteration_path(i_iteration)
        mc_path = os.path.join(iteration_path, mc_type)

        return mc_path 

    def get_mutate_type(self, i_iteration: int) -> str:
        iteration_path = self.get_iteration_path(i_iteration)
        mutate_type_path = os.path.join(iteration_path, 'mutate_type')
        with open(mutate_type_path) as f:
            mutate_type = f.read()
        return mutate_type

    def get_is_accepted_path(self, i_iteration: int) -> str:
        iteration_path = get_iteration_path(i_iteration)
        is_accepted_path = os.path.join(iteration_path, 'is_accepted')
        return is_accepted_path

    def get_multicell_paths(self, i_iteration) -> Dict[str, str]:
        mc_paths = {}
        for mc_type in Pymatmc2Results.valid_mc_types:
            mc_paths[mc_type] = self.get_mc_path(i_iteration, mc_type)
        return mc_paths

    def get_multicells(self, i_iteration: int) -> Dict[str, MultiCell]:
        mc_paths = {}
        for mc_type in Pymatmc2Results.valid_mc_types:
            mc_paths[mc_type] = self.get_mc_path(i_iteration, mc_type)
       
        multicells = {} 
        for mc_type, mc_path in mc_paths.items():
            multicells[mc_type] = MultiCell()
            multicells[mc_type].configuration = self.configuration
            multicells[mc_type].read(path=mc_path)

        return multicells

    def get_iteration_info(self, i_iteration: int):
        mc = self.get_multicells(i_iteration)

        iteration_info = {}
        iteration_info['iteration'] = i_iteration
        iteration_info['mutate_type'] = self.get_mutate_type(i_iteration)
        for k in mc:
            for cn in mc[k].cell_names:
                for s in mc[k].symbols:
                    iteration_info['{}.{}.{}'.format(k,cn,s)] = mc[k].cell_concentration[cn][s]
        for k in mc:
            for cn, f in mc[k].phase_molar_fraction.items():
                iteration_info['{}.f.{}'.format(k,cn)] = f
        for k in mc:
            for cn in mc[k].cell_names:
                iteration_info['{}.{}.toten'.format(k,cn)] = mc[k].simulations[cn].total_energy
            iteration_info['{}.all.toten'.format(k)] = mc[k].total_energy
        n_phases = len(self.configuration.cell_names)
        T = self.configuration.temperature
        kB = constants.BOLTZMANN
        beta = 1/T/kB
        
        for cn in self.configuration.cell_names + ['all']:
            if iteration_info['mutate_type'] == 'intraphase_flip':
                if cn == 'all':
                    E1 = mc['initial'].total_energy
                    E2 = mc['candidate'].total_energy
                    pacc = min(1,np.exp(-n_phases*(E2-E1)*beta))
                    iteration_info['{}.p_acc'.format(cn)] = pacc
                else:
                    iteration_info['{}.p_acc'.format(cn)] =  None
            elif iteration_info['mutate_type'] == 'intraphase_swap': 
                if cn == 'total':
                    iteration_info['{}.p_acc'.format(cn)] = None
                else:
                    E1 = mc['initial'].simulations[cn].total_energy
                    E2 = mc['candidate'].simulations[cn].total_energy
                    pacc = min(1,np.exp(-n_phases*(E2-E1)*beta))
                    iteration_info['{}.p_acc'.format(cn)] = pacc
            else:
                msg = 'unknown mutate_type: {}'
                msg = msg.format(iteration_info['mutate_type'])
                raise ValueError(msg)

        return iteration_info

# need to integrate this back in
def print_table(iteration_info):
    total_energy = {}
    total_energy
    for k, v in mc.items():
        total_energy[k] = {}
        for cn in configuration.cell_names:
            total_energy[k][cn] = v.simulations[cn].total_energy
        total_energy[k]['total'] = v.total_energy

        
    cell_concentration = {}
    for k, v in mc.items():
        cell_concentration[k] = v.cell_concentration


    phase_molar_fraction = {}
    for k, v in mc.items():
        phase_molar_fraction[k] = v.phase_molar_fraction 
    
    total_energy_per_atom = {}
    for k, v in mc.items():
        total_energy_per_atom[k] = {}
        for cn in configuration.cell_names:
            E = v.simulations[cn].total_energy
            n_atoms = v.simulations[cn].poscar.n_atoms
            total_energy_per_atom[k][cn] = E/n_atoms

        E = 0
        for cn in configuration.cell_names:
            E += phase_molar_fraction[k][cn] * total_energy_per_atom[k][cn]
        n_phases = len(configuration.cell_names)
        total_energy_per_atom[k]['total'] = E * n_phases

    is_intrinsic = False
    if mutate_type == 'intraphase_flip':
        pacc= {}
        if is_intrinsic:
           n_atoms = sum([v.poscar.n_atoms for v in mc['candidate'].simulations.values()])
           Ec = total_energy['candidate']['total']/n_atoms
           Ei = total_energy['initial']['total']/n_atoms
           pacc['total'] = min(1, np.exp(-n_phases*(Ec-Ei)*beta))
        else:
           Ec = total_energy['candidate']['total']
           Ei = total_energy['initial']['total']
           pacc['total'] = min(1, np.exp(-n_phases*(Ec-Ei)*beta))
    elif mutate_type =='intraphase_swap':
        pacc = {}
        for cn in configuration.cell_names:
            if is_intrinsic:
                n_atoms = mc['candidate'].simulations[cn].poscar.n_atoms
                Ec = total_energy['candidate'][cn]/n_atoms
                Ei = total_energy['initial'][cn]/n_atoms
                pacc[cn] = min(1, np.exp(-beta*(Ec-Ei)))
            else:
                Ec = total_energy['candidate'][cn]
                Ei = total_energy['initial'][cn]
                pacc[cn] = min(1, np.exp(-beta*(Ec-Ei)))

    #print table
    fmt_header = '{:15} {:15}' + n_symbols * '{:^10}'   + '{:^10}'   + '{:^10}'   + '{:^10}' + '{:^10}'
    fmt_row    = '{:15} {:15}' + n_symbols * '{:10.4f}' + '{:10.4f}' + '{:10.4f}' + '{:10.4f}'

    header = fmt_header.format('','',*symbols, 'f', 'E', 'E_per_atom', 'p_acc')
    print(header)
    for k, mc in cell_concentration.items():
         for cn, v in mc.items():
             c = [v[s] for s in symbols]
             row_values = []
             row_values.append(k)
             row_values.append(cn)
             row_values += c
             row_values.append(phase_molar_fraction[k][cn])
             row_values.append(total_energy[k][cn])
             row_values.append(total_energy_per_atom[k][cn])
             row = fmt_row.format(*row_values)
             if cn in pacc.keys() and k == 'candidate':
                 row += '{:10.4f}'.format(pacc[cn])
             print(row)

         # total row
         cn = 'total'
         c = list(configuration.total_concentration.values())
         row_values = []
         row_values.append(k)
         row_values.append(cn)
         row_values += c
         row_values.append(1.0)
         row_values.append(total_energy[k][cn])
         row_values.append(total_energy_per_atom[k][cn])
         row = fmt_row.format(*row_values)
         if cn in pacc.keys() and k == 'candidate':
             row += '{:10.4f}'.format(pacc[cn])
         print(row)
