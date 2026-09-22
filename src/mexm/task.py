from abc import ABC
from typing import List
import os

class Task():
    def __init__(
        self,
        name: str,
        simulation_type: str,
        simulation_path: str,
        n_cores: int,
        n_nodes: int,
        required: List[str] = {}
    ):
        """

        Arguments:
            task_name (str): the name of the task
            simulation_type (str): the type of simulation
            simulation_path (str): the path to the simulation, will 
                be converted to an absolute path
            n_cores (int): the number of processor cores to use
            n_nodes (int): the number of processor nodes to use
            request (List[str]): a list of task_names that this 
                dependent upon
        """
        self.name = name
        self.simulation_type = simulation_type
        self._simulation_path = simulation_path
        self.n_cores = n_cores
        self.n_nodes = n_nodes
         

    @property
    def simulation_path(self):
        return self._simulation_path

    @simulation_path.setter
    def simulation_path(self, path):
        self._simulation_path = os.path.abspath(path)



        

class TaskManager(ABC):
    def __init__(self, sql_adapter: DbAdapter):
        self.task_repository = {}
        self.sql_adapter = sql_adapter


    def submit_task(self):
        pass


class SerialTaskManager(TaskManager): pass
class MpiTaskManager(TaskManager): pass
