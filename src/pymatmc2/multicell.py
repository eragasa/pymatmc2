from collections import OrderedDict
from mexm.io.vasp import Poscar
from pymatmc2 import Pymatmc2Configuration

class MultiCell:
    """

    Attributes:
        n_cells (str)
        cells (Dict[str, Poscar])
    """

    def __init__(self):
        """

        Arguments:
            n_cells (int): number of cells
        """
        self.n_cells = None
        self.cells = None

    @staticmethod
    def initialize_from_pymatmc2_configuration(
        configuration: Pymatmc2Configuration
    ):
        """
            configuration (Pymatmc2Configuration)
        """
        if isinstance(configuration, Pymatmc2Configuration):
            obj = MultiCell()
            obj.cells = {}
            for k,v in configuration.simulation_cells.items():
                obj.cells[k] = Poscar()
                obj.cells[k].read(path=v['poscar'])
            return obj
    
    def configure(self, configuration: Pymatmc2Configuration):
        """ configure class from a Pymatmc2Configuration
        
        this method sets up the cells attribute

        Arguments:
            configuration (Pymatmc2Configuration): instance of the 
                Pymatmc2Configuration in which to setup this class.
         """
        self.cells = {}
        for k,v in configuration.simulation_cells.items():
            self.cells[k] = Poscar()
            self.cells[k].read(path=v['poscar'])            

    def add_cells_from_dict(self, simulation_cells):
        """ add cells from dictionary object

        simulation_cells = {
            '<cell_name>' = {'<file_type'>:'<path_to_file>'}
        }
        Arguments:
            simulation_cells (dict) = dictionary of simulation cells
        """

    @property
    def symbols(self):
        symbols = []
        for cell in self.cells:
            for symbol in cell.symbols:
                if symbol not in symbols:
                    symbols.append(symbol)
        return symbols

    def get_number_of_atoms(self, symbol=None):
        n_atoms = 0
        for cell in self.cells:
            n_atoms += cell.get_number_of_atoms(symbol)
        return n_atoms
