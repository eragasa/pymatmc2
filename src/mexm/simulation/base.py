import os
import shutil
from copy import deepcopy
from typing import Dict
from abc import ABC
from abc import abstractmethod
from abc import abstractclassmethod
from collections import OrderedDict

from mexm.structure import SimulationCell
from mexm.exception import MexmException

class MexmSimulationException(MexmException): pass

class Simulation(ABC):

    states = ['INIT','CONFIG','READY','RUNNING','POST','FINISHED','ERROR']
    def __init__(self,
                 name,
                 path):

        # private variables
        self._name = None
        self._path = None

        # constructor arguments
        self.name = name
        self.path = path

        # initialize conditions
        self.conditions = {k:{} for k in Simulation.states}

        # create simulation directory
        self.create_path(path=path)
        self.status = None

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        if not isinstance(name, str):
            raise TypeError("name must be a string")

        bad_symbols = [' ']
        if any([k in name for k in bad_symbols]):
            raise ValueError("name has an illegal character")

        self._name = name

    @property
    def structure(self) -> SimulationCell:
        raise NotImplementedError

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, path):
        if not isinstance(path, str):
            raise TypeError("path must be a string")

        if os.path.isabs(path):
            self._path = path
        else:
            self._path = os.path.abspath(path)

    @abstractmethod
    def get_conditions_init(self): 
        raise NotImplementedError
    
    @abstractmethod
    def get_conditions_config(self): 
        raise NotImplementedError
    
    @abstractmethod
    def get_conditions_ready(self): 
        raise NotImplementedError
    
    @abstractmethod
    def get_conditions_running(self): 
        raise NotImplementedError
    
    @abstractmethod
    def get_conditions_post(self): 
        raise NotImplementedError
    
    @abstractmethod
    def get_conditions_finished(self): 
        raise NotImplementedError

    @abstractmethod
    def on_init(self): 
        raise NotImplementedError
    
    @abstractmethod
    def on_config(self): 
        raise NotImplementedError
    
    @abstractmethod
    def on_ready(self): 
        raise NotImplementedError
    
    @abstractmethod
    def on_running(self): 
        raise NotImplementedError
    
    @abstractmethod
    def on_post(self): 
        raise NotImplementedError
    
    @abstractmethod
    def on_finished(self): 
        raise NotImplementedError
    
    @abstractmethod
    def on_error(self): 
        raise NotImplementedError


    @property
    def conditions_INIT(self):
        return self.conditions['INIT']

    @conditions_INIT.setter
    def conditions_INIT(self, conditions) -> Dict:
        self.conditions['INIT'] = deepcopy(conditions)

    @property
    def conditions_CONFIG(self): 
        return self.conditions['CONFIG']

    @conditions_CONFIG.setter
    def conditions_CONFIG(self, conditions) -> Dict:
        self.conditions['CONFIG'] = deepcopy(conditions)

    @property
    def conditions_READY(self) -> Dict: 
        return self.conditions['READY']

    @conditions_READY.setter
    def conditions_READY(self, conditions):
        self.conditions['READY'] = deepcopy(conditions)

    @property
    def conditions_RUNNING(self): return self.conditions['RUNNING']

    @conditions_RUNNING.setter
    def conditions_RUNNING(self, conditions):
        self.conditions['RUNNING'] = deepcopy(conditions)

    @property
    def conditions_POST(self): return self.conditions['POST']

    @conditions_POST.setter
    def conditions_POST(self, conditions):
        self.conditions['POST'] = deepcopy(conditions)

    @property
    def conditions_FINISHED(self): return self.conditions['FINISHED']

    @conditions_FINISHED.setter
    def conditions_FINISHED(self, conditions):
        self.conditions['FINISHED'] = deepcopy(conditions)

    @property
    def conditions_ERROR(self): return self.conditions['ERROR']

    @conditions_ERROR.setter
    def conditions_ERROR(self, conditions):
        self.conditions['ERROR'] = deepcopy(conditions)

    def read(self):
        raise NotImplementedError

    def write(self, path):
        self.create_path(path=path)

    def run(self): 
        raise NotImplementedError

    def create_path(self, path=None):
        
        if path is not None:
            self.path = path
        path_ = self.path

        cwd = os.path.abspath(os.getcwd())
        if self.path == cwd:
            msg = 'path cannot be set to current working directory'
            raise MexmSimulationException(msg)

        # remove existing directory, if directory exists
        if os.path.isdir(self.path):
            shutil.rmtree(self.path, ignore_errors=True)

        # this is to deal with possible race conditions while deleting directory
        while True:
            try:
                os.mkdir(self.path)
                break
            except FileExistsError as e:
                if os.path.isdir(path_):
                    shutil.rmtree(path_, ignore_errors=True)
                    if e.errno != os.errno.EEXIST:
                        raise
                    pass

    def on_update_status(self):
        status_to_method_map = {
            'INIT':self.on_init,
            'CONFIG':self.on_config,
            'READY':self.on_ready,
            'RUNNING':self.on_running,
            'POST':self.on_post,
            'FINISHED':self.on_finished,
            'ERROR':self.on_error
        }
        status_to_method_map[self.status]()


    def update_status(self):
        assert isinstance(self.conditions_INIT, dict)
        if self.status is None:
            self.get_conditions_init()
            if all([v for k,v in self.conditions_INIT.items()]):
                self.status = 'INIT'
                return
            else:
                for k, v in self.conditions_INIT.items():
                    print(k, v)
                msg = "was not able to initialize"
                raise ValueError(msg)
        elif self.status == 'INIT':
            self.get_conditions_config()
            if all([v for k,v in self.conditions_CONFIG.items()]):
                self.status = 'CONFIG'
                return
            else:
                return
        elif self.status == 'CONFIG':
            self.get_conditions_ready()
            if all([v for k,v in self.conditions_READY.items()]):
                self.status = "READY"
                return
            else:
                return
        elif self.status == 'READY':
            self.get_conditions_running()
            if all([v for k,v in self.conditions_RUNNING.items()]):
                self.status = "RUNNING"
                return
            else:
                return
        elif self.status == 'RUNNING':
            self.get_conditions_post()
            if all([v for k,v in self.conditions_POST.items()]):
                self.status = "POST"
                return
            else:
                return
        elif self.status == 'POST':
            self.get_conditions_finished()
            if all([v for k,v in self.conditions_FINISHED.items()]):
                self.status = "FINISHED"
                return
            else:
                return
        elif self.status == 'ERROR':
            raise ValueError()
        else:
            msg = "unknown status, {}".format(self.status)
            raise ValueError(msg)


class AtomicSimulation(): pass
