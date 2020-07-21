from typing import List
import pandas as pd 
from pymatmc2 import MultiCell
from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import ChemicalPotential

class Pymatmc2Data():
    """

    Attributes:
        df (pandas.DataFrame)
    """
    def __init__(self):
        self._configuration = None
        self.df = None
    @property
    def configuration(self) -> Pymatmc2Configuration:
        return self._configuration 

    @configuration.setter
    def configuration(self, configuration: Pymatmc2Configuration):
        if not isinstance(configuration, Pymatmc2Configuration):
            msg = "configuration attribute must be an instance of pymatmc2.Pymatmc2Configuration"
            raise TypeError(msg)
        self._configuration = configuration

    @property
    def column_names(self) -> List[str]:
        if not isinstance(self.configuration, Pymatmc2Configuration):
            msg = 'cannot determine column names if the configuration property is not set first'
            raise ValueError(msg)

        column_names = ['iteration']
        column_names += ['mutate_type']
        column_names += ['{}.accept'.format(k) for k in self.configuration.cell_names]
        column_names += ['{}.accept'.format('all')]
        column_names += ['{}.p_acc'.format(k) for k in self.configuration.cell_names]
        column_names += ['{}.p_acc'.format('all')]
        column_names += ['{}.p_acc_rnd'.format(k) for k in self.configuration.cell_names]
        column_names += ['{}.p_acc_rnd'.format(k) for k in self.configuration.cell_names]
        column_names += ['initial.{}.toten'.format(k) for k in self.configuration.cell_names]
        column_names += ['initial.all.toten']
        column_names += ['candidate.{}.toten'.format(k) for k in self.configuration.cell_names]
        column_names += ['candidate.all.toten']
        column_names += ['final.{}.toten'.format(k) for k in self.configuration.cell_names]
        column_names += ['final.all.toten']

        for cn in self.configuration.cell_names:
            for s in self.configuration.symbols:
                column_name = 'initial.{}.{}'.format(cn, s)
                column_names.append(column_name)

        for cn in self.configuration.cell_names:
            column_name = 'initial.f.{}'.format(cn)
            column_names.append(column_name)

        for cn in self.configuration.cell_names:
            for s in self.configuration.symbols:
                column_name = 'candidate.{}.{}'.format(cn, s)
                column_names.append(column_name)

        for cn in self.configuration.cell_names:
            column_name = 'candidate.f.{}'.format(cn)
            column_names.append(column_name)
        
        for cn in self.configuration.cell_names:
            for s in self.configuration.symbols:
                column_name = 'final.{}.{}'.format(cn, s)
                column_names.append(column_name)

        for cn in self.configuration.cell_names:
            column_name = 'final.f.{}'.format(cn)
            column_names.append(column_name)
        
        for cn in self.configuration.cell_names:
            symbols = self.configuration.symbols
            bond_types = ChemicalPotential.get_bond_types(symbols)
            for bond_type in bond_types:
                column_name = 'final.{}.{}'.format(cn, bond_type)
                column_names.append(column_name)
        return column_names

    def write_header(self):
        results_data_path = self.configuration.results_file_path
        
        s = ",".join(self.column_names) + "\n"
        with open(results_data_path, 'w') as f:
            f.write(s)
        
    def write_results(self, iteration_info):
        results_data_path = self.configuration.results_file_path
        row_values = []
        for k in self.column_names:
            if k in iteration_info:
                row_values.append(iteration_info[k])
            else:
                row_values.append(None)

        s = ",".join([str(k) for k in row_values]) + "\n"

        with open(results_data_path, 'a') as f:
            f.write(s)

    def read(self, path: str):
        self.df = pd.read_csv(path)

