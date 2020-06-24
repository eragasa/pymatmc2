from typing import Optional
from typing import Union
from collections import OrderedDict

import matplotlib.pyplot as plt

from pymatmc2 import Pymatmc2Configuration
from pymatmc2 import Pymatmc2Data
from pymatmc2 import Pymatmc2Results

class Pymatmc2Plots:
    def __init__(self,
        configuration:Optional[Union[Pymatmc2Configuration, str]] = None,
        data:Optional[Union[Pymatmc2Data, str]] = None
    ):
        self._configuration = None
        self.configuration = configuration

        self._data = None
        self.data = data

    @property
    def configuration(self):
        return self._configuration
    
    @configuration.setter
    def configuration(self, configuration:Optional[Pymatmc2Configuration]):
        if configuration is None:
            self._configuration = None
        elif isinstance(configuration, str):
            self._configuration = Pymatmc2Configuration()
            self._configuration.read(configuration)
        elif isinstance(configuration, Pymatmc2Configuration):
            self._configuration = configuration
        else:
            m = 'configuration must either be None, '
            raise TypeError(m)

    @property
    def data(self):
        return self._data
    
    @property.setter
    def data(self, data:Union[type(None), Pymatmc2Data, str]):
        if data is None:
            self._data = None
        elif isinstance(data, str):
            self._data = Pymatmc2Data()
            self._data.configuration = self.configuration
            self._data.read(data)
        elif isinstance(data, Pymatmc2Data):
            self._data = data
        else:
            m = 'configuration must either be None, str, or Pymatmc2Data'
            raise TypeError(m)
