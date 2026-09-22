from abc import ABC
from abc import abstractmethod
import os

class HpcClusterInformation(ABC):
    def __init__(self):
        self.name = None
        self.cores_per_node = None
        self.ram_per_node = None

class HpcJobInformation(ABC):
    pass

class HpcSubmissionScript(ABC):

    def __init__(self):
        self.job_name = None

    @abstractmethod
    def write(self, path):
        raise NotImplementedError


class JobSubmissionManager(ABC):

    @abstractmethod
    def submit_job(self, simulation_path: str, submission_script_path: str):
        raise NotImplementedError

    @abstractmethod
    def get_job_info(self, jobid='all', username=None):
        raise NotImplementedError

from mexm.job.torque import TorqueSubmissionScript
from mexm.job.torque import TorqueJobSubmissionManager

class JobSubmissionManagerFactory():
    obj_submission_scripts = {
        'torque':TorqueSubmissionScript
    }
    obj_job_submission_manager = {
        'torque':TorqueJobSubmissionManager
    }
    
    @staticmethod
    def write_submission_script(hpc_type, script_kwargs, script_path):
        obj_submission_script \
            = JobSubmissionManagerFactory.obj_submission_scripts[hpc_type](
                **script_kwargs
            )
        obj_submission_script.write(path=script_path)

    @staticmethod
    def submit_job(hpc_type, simulation_path, submission_script_path='runjob.sh'):
        obj = JobSubmissionManagerFactory.obj_job_submission_manager[hpc_type]()
        obj.submit_job(
            simulation_path = simulation_path,
            submission_script_path = submission_script_path
        )
