""" module to manage jobs to slurm """
import os, yaml
from copy import copy, deepcopy
from collections import OrderedDict

from mexm.io.hpcutil import HpcCluster
from mexm.io.torque.configuration import TorqueConfiguration
from mexm.io.torque.submissionscript import TorqueSubmissionScript

class TorqueCluster(HpcCluster):
    environment_variables = {
            'version':'$PBS_VERSION',
            'job_name':'$PBS_JOBNAME',
            'job_id':'$PBS_JOBID',
            'submit_path':'$PBS_O_WORKDIR',
            'submit_host':'$PBS_O_HOST',
            'walltime':'$PBS_WALLTIME'
    }

    def __init__(self):
        pass

    def submit_job(self, job_script_path):
        cmd = 'qsub {}'.format(job_script_path)

    def delete_job(self, job_id):
        cmd = 'qdel {}'.format(job_id)

    def get_job_status_all(self):
        cmd = 'qstat'

    def get_job_status_by_job(self, job_id):
        cmd = 'qstat {}'.format(job_id)

    def get_job_status_by_user(self, user=None):
        if user is None:
            user_ = self.user
        cmd = 'qstat -u {}'.format(user_)

    def get_job_details(self, job_id):
        cmd = 'qstat -f {}'.format(job_id)
        cmd = 'checkjob {}'.format)job_id)

class VaspSubmissionScript(TorqueSubmissionScript):
    """ class for creating SLURM submission scripts for VASP

    These classes do require that a SlurmConfiguration file have been already
    created, and provided as the path argument.
    """

    def __init__(self, path=None):
        SlurmSubmissionScript.__init__(self, path)
        self.modules = self.configuration.vasp_configuration['modules']
        self.run_string = self.configuration.vasp_configuration['run_string']

class LammpsSubmissionScript(TorqueSubmissionScript):
    """ class for creating SLURM submission scripts for LAMMPS

    These classes do require that a SlurmConfiguration file have been already
    created, and provided as the path argument.
    """

    def __init__(self, path=None):
        SlurmSubmissionScript.__init__(self, path)
        self.modules = self.configuration.lammps_configuration['modules']
        self.run_string = self.configuration.lammps_configuration['run_string']

class PhontSubmissionScript(TorqueSubmissionScript):
    """ class for creating SLURM submission scripts for PHONTS

    These classes do require that a SlurmConfiguration file have been already
    created, and provided as the path argument.
    """
    def __init__(self, path=None):
        SlurmSubmissionScript.__init__(self,path)
        self.modules = slef.configuration.phonts_configuration['modules']
        self.run_string = self.configuration.phonts_configuration['run_string']
