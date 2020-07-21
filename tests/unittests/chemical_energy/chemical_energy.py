from collections import OrderedDict
from typing import List
from typing import Dict
import numpy as np
import scipy.stats

symbols = ['Au', 'Pt']

cell_compositions_nominal = [
    OrderedDict([('Au', 4), ('Pt', 0)]),
    OrderedDict([('Au', 0), ('Pt', 4)]),
    OrderedDict([('Au', 0), ('Pt', 4)]),
    OrderedDict([('Au', 0), ('Pt', 4)])
]

cell_bond_counts = [
    OrderedDict([('Au.Au', 24), ('Au.Pt',  0), ('Pt.Pt',  0)]),
    OrderedDict([('Au.Au',  0), ('Au.Pt',  0), ('Pt.Pt', 24)]),
    OrderedDict([('Au.Au', 12), ('Au.Pt', 12), ('Pt.Pt',  0)]),
    OrderedDict([('Au.Au',  0), ('Au.Pt', 12), ('Pt.Pt', 12)])
]

cell_energies = [
    -13.08,
    -24.44,
    -15.71,
    -21.20
]

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
                bond_types.append('{}.{}'.format(s1, s2))
    return bond_types

def calculate_bond_energies(cell_bonds:List[Dict[str, int]], cell_energies: List[float]) -> Dict[str, float]:
    pass
if __name__ == "__main__":

    bond_energies = OrderedDict()
    cell_bond_count_ = []
    for bond_count in cell_bond_counts:
        cell_bond_count_.append([v for k, v  in bond_count.items()])
    print(cell_bond_counts)
    matrix_cell_bond_count = np.array(cell_bond_count_)

    print(cell_energies)
    print(matrix_cell_bond_count)

    import statsmodels.api as sm
    from statsmodels.sandbox.regression.predstd import wls_prediction_std

    is_add_constant = False
    X = matrix_cell_bond_count
    Y = np.array(cell_energies)

    bond_energy_model = sm.OLS(Y, X)
    results = bond_energy_model.fit()
    print(results.summary())
