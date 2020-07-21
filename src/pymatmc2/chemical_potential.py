from collections import OrderedDict
from typing import List
from typing import Dict
from typing import Optional
from typing import Tuple
import numpy as np
import scipy.stats
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from mexm.structure import SimulationCell

class ChemicalPotential():

    @staticmethod
    def get_bond_types(symbols:List[str]) -> List[str]:
        """ get a list of bond types from a list of symbols

        Arguments:
            symbols (List[str]): list of symbols

        Returns
            List[str]: list of the bond types, delineated by a `.`
        """

        bond_types = []
        for i1, s1 in enumerate(symbols):
            for i2, s2 in enumerate(symbols):
                if i1 <= i2:
                    bond_types.append('{}{}'.format(s1, s2))
        return bond_types

    @staticmethod
    def get_distance(
        cell: SimulationCell, 
        index_1: int, 
        index_2:int, 
        is_debug=False, 
        max_distance:Optional[float]=1e6
    ):
        """ Calculates distance between two cells using the minimum image convention

        """
        x1 = np.array(cell.atomic_basis[index_1].position) * cell.lattice.a0
        x2 = np.array(cell.atomic_basis[index_2].position) * cell.lattice.a0

        d1 = max_distance
        for ia1 in [-1,0,1]:
            for ia2 in [-1,0,1]:
                for ia3 in [-1,0,1]:
                    shift = ia1 * cell.a1 + ia2 * cell.a2 + ia3 * cell.a3
                    x2_ = x2 + shift
                    d1_ = np.dot(x1-x2_, x1-x2_)
                    d1 = min(d1, d1_)
        return d1 

    @staticmethod
    def calculate_bond_energies(
        cell_bond_counts:List[Dict[str, int]], 
        cell_energies: List[float],
        is_debug: Optional[bool] = False
        ) -> Tuple[Dict[str, float], Dict[str, float]]:

        if is_debug:
            print('cell_energies:')
            print(cell_energies)

        # return object
        cell_bond_count_ = []
        for bond_count in cell_bond_counts:
            cell_bond_count_.append([v for k, v  in bond_count.items()])
        if is_debug:
            print('cell_bond_counts:')
            print(cell_bond_counts)

        matrix_cell_bond_count = np.array(cell_bond_count_)
        if is_debug:
            print('matrix_cell_bond_count:')
            print(matrix_cell_bond_count)

        is_add_constant = False
        X = matrix_cell_bond_count
        Y = np.array(cell_energies)

        bond_energy_model = sm.OLS(Y, X)
        results = bond_energy_model.fit()
        if is_debug:
            print(results.summary())
            print(results.params)
        
        bond_types = [k for k, v in cell_bond_counts[0].items()]
        if is_debug:
            print(bond_types)

        bond_energies = OrderedDict()
        bond_energies_std = OrderedDict()
        for i, v in enumerate(bond_types):
            bond_energies[v] = results.params[i]
            bond_energies_std[v] = results.bse[i]
        return bond_energies, bond_energies_std

    @staticmethod
    def get_bond_lengths(cell: SimulationCell, is_debug: Optional[bool] = False) -> np.ndarray:

        bond_lengths = []
        for i in range(cell.n_atoms):
            for j in range(cell.n_atoms):
                if i < j:
                    d = ChemicalPotential.get_distance(cell, i, j)
                    if is_debug:
                        print('i:{}, j:{}, d:{}'.format(i, j, d))
                    bond_lengths.append(d)
        bond_lengths = np.array(bond_lengths)
        return bond_lengths

    @staticmethod
    def get_bond_lengths_kde(cell: SimulationCell, is_debug:Optional[bool] = False) -> gaussian_kde:

        assert issubclass(type(cell), SimulationCell)

        bond_lengths = ChemicalPotential.get_bond_lengths(cell=cell)
        kde = gaussian_kde(bond_lengths, bw_method=0.05)
        return kde

    @staticmethod
    def get_cell_bond_counts(cell: SimulationCell, is_debug:Optional[bool]=False) -> np.ndarray:

        bond_lengths = ChemicalPotential.get_bond_lengths(cell)
        kde = ChemicalPotential.get_bond_lengths_kde(cell)
        
        # get peaks
        x = np.linspace(0, max(bond_lengths)*1.25, 1000)
        bond_length_kde = kde(x)
        peaks_idx, peak_properties = find_peaks(kde(x), height=0)

        # determine r cutoff
        r1 = x[peaks_idx[0]]
        r2 = x[peaks_idx[1]]
        r_cutoff = (r1 + r2)/2

        if is_debug:
            print('r0:{}'.format(r1))
            print('r1:{}'.format(r2))
            print('r_cutoff:{}'.format(r_cutoff))
            import matplotlib.pyplot as plt
            plt.plot(x, bond_length_kde)
            plt.plot(x[peaks_idx], bond_length_kde[peaks_idx], 'x', color='black')
            plt.show()

        cell_bond_counts = OrderedDict()
        bond_types = ChemicalPotential.get_bond_types(cell.symbols)
        for i in range(cell.n_atoms):
            for j in range(cell.n_atoms):
                d = ChemicalPotential.get_distance(cell, i, j)
                if d < r_cutoff:
                    s1 = cell.atomic_basis[i].symbol
                    s2 = cell.atomic_basis[j].symbol

                    bond_type = '{}{}'.format(s1, s2)
                    if bond_type not in bond_types:
                        bond_type = '{}{}'.format(s2, s1)

                    if bond_type not in cell_bond_counts:
                        cell_bond_counts[bond_type] = 1
                    else:
                        cell_bond_counts[bond_type] += 1
        n_bonds = 0
        for k, v in cell_bond_counts.items():
            n_bonds += v
        coordination_number = n_bonds/cell.n_atoms
        print('n_atoms:{}'.format(cell.n_atoms))
        print('coordination_number:{}'.format(coordination_number))
        return cell_bond_counts
if __name__ == "__main__":
    pass
