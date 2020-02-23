import os
import copy
import yaml
from mexm.io.filesystem import OrderedDictYAMLLoader

class Pymatmc2Configuration():
    """ configuration file for pymatmc2

    Attributes:
        path (str): the path to the directory
        configuration (dict): dictionary of the configuration
    """
    def __init__(self):
        self._path = None
        self.configuration = None
        pass

    @staticmethod
    def initialize_from_dict(config_dict: dict):
        """
        Arguments:
            config_dict (dict): configuration dictionary
        """
        obj = Pymatmc2Configuration()
        obj.configure_from_dict(config_dict=config_dict)

    def configure_from_dict(self, config_dict: dict):
        self.configuration = copy.deepcopy(config_dict)
        
    @property
    def path(self):
        """(str): absolute path of the configuration file"""
        
        return self._path

    @path.setter
    def path(self, path):

        # convert to absolute
        if isinstance(path, str):
            if os.path.isabs(path):
                self._path = path
            else:
                self._path = os.path.abspath(path)
        else:
            msg = "path must be a string"
            raise TypeError(msg)

    @property
    def n_cells(self):
        """:(int) number of simulation_cells """     
        return len(self.configuration['atomic_configuration']['simulation_cells'])

    @property
    def simulation_cells(self):
        """:dict simulations cells"""
        return self.configuration['atomic_configuration']['simulation_cells']
    
    @property
    def temperature(self):
        return self.configuration['environment_variables']['temperature']

    @property
    def pressure(self):
        return self.configuration['environment_variables']['pressure']

    def read(self, path):
        """ read the configuration files

        Arguments:
            path (str): path to the configuration file
        """
        self.path = path

        try:
            with open(self.path, 'r') as f:
                self.configuration= yaml.load(f, OrderedDictYAMLLoader)
        except FileNotFoundError:
            raise

    def write(self, path):
        """ write the configuration files

        Arguments:
            path (str): path to the configuration file
        """
        self.path = path

        with open(self.path, 'w') as f:
            yaml.dump(
                self.configuration, 
                f,
                default_flow_style=False
            )

