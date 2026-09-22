import os
import subprocess
import math
from copy import deepcopy
from typing import List
from mexm.simulation import VaspSimulation
from mexm.job import HpcClusterInformation
from mexm.job import HpcJobInformation
from mexm.job import HpcSubmissionScript
from mexm.job import JobSubmissionManager

class TorqueHpcClusterInformation(HpcClusterInformation):
    pass

class TorqueJobInformation(HpcJobInformation):
    def __init__(
        self,
        account: str,
        walltime: int, 
        nodes: int,
        ppn: int,
        jobname: str,
        errpath='job.err',
        stdpath='job.out',
        modules=None
    ):
        """
        Arguments:
            walltime (float): in hours
        """
        super().__init__()
        self.account = account
        self.walltime = walltime
        self.nodes = nodes
        self.ppn = ppn
        self.jobname = jobname
        self.errpath = errpath
        self.stdpath = stdpath
        self.modules = modules

class TorqueSubmissionScript(HpcSubmissionScript):

    def __init__(
        self,
        account=None,
        walltime=None,
        n_nodes=None,
        ppn=None,
        jobname=None,
        errpath=None,
        stdpath=None,
        modules=None,
        cmd=None,
    ):
        super().__init__()
        self.account=account
        self.walltime=walltime
        self.n_nodes=n_nodes
        self.ppn=ppn
        self.jobname=jobname
        self.errpath=errpath
        self.stdpath=stdpath
        self.modules=modules
        self.cmd = cmd
    
    def read(self, path): 
        pass

    def write(self, path):
        submission_script_str = "\n".join([
            self.header_section_to_string(
                account=self.account,
                walltime=self.walltime_to_string(self.walltime),
                n_nodes=self.n_nodes,
                ppn=self.ppn,
                jobname=self.jobname,
                errpath=self.errpath,
                stdpath=self.stdpath
            ),
            self.module_section_to_string(modules=self.modules),
            self.cmd,
            self.footer_section_to_string()
        ])

        with open(path,'w') as f:
            f.write(submission_script_str)
        return submission_script_str
    
    def walltime_to_string(self, walltime: int):
        hr = str(walltime)
        min = '00'
        sec = '00'
        return "{}:{}:{}".format(hr,min,sec)

    def header_section_to_string(
        self,
        account: str, 
        walltime: int, 
        n_nodes: int, ppn, jobname, errpath, stdpath
    ):
        header_str = (
            "#PBS -A {account}\n"
            "#PBS -l walltime={walltime}\n"
            "#PBS -l nodes={n_nodes}:ppn={ppn}\n"
            "#PBS -N {jobname}\n"
            "#PBS -e {errpath}\n"
            "#PBS -o {stdpath}\n"
            "#PBS -S /bin/bash\n\n"
            "cd $PBS_O_WORKDIR\n\n"
        ).format(
            account=account,
            walltime=walltime,
            n_nodes=n_nodes,
            ppn=ppn,
            jobname=jobname,
            errpath=errpath,
            stdpath=stdpath
        )
        return header_str

    def module_section_to_string(self, modules: List[str]) -> str:
        modules_str = "\n".join(
            ["module load {}".format(k) for k in modules]
        ) + "\n\n"
        return modules_str

    def footer_section_to_string(self):
        return "touch jobComplete"

class TorqueJobSubmissionManager(JobSubmissionManager):
    def __init__(
        self,
        qsub_path='/opt/torque/bin/qsub'
    ):
        self.cluster = TorqueHpcClusterInformation()
        self.qsub_path=qsub_path

    def request_job(
        self,
        simulation,
        n_cores,
        n_nodes,
        modules = [],
        ppn=None,
        ram=None
    ):
        if ppn is None:
            ppn_ = ppn
        else:
            ppn_ = self.cluster.cores_per_node

        if n_cores % ppn_ != 0:
            n_nodes_ =  (math.floor(n_cores/ppn_) + 1)
            n_cores_ = n_nodes * ppn_
        else:
            n_nodes_ = n_nodes
            n_cores_ = n_cores

        if isinstance(simulation, VaspSimulation):
            simulation.update_hpc_configuration(
                n_cores,
                n_nodes
            )

        if isinstance(modules, list):
            module_str = [
                'module load {}'.format(k) for k in modules
            ]


    def insert_job_into_database(self,
        simulation_name,
        simulation_scipt
    ):
        pass

    def submit_job(
        self, 
        simulation_path: str,
        submission_script_path: str
    ):
        # setting an initial context, so I can return to it
        initial_path = os.getcwd()

        os.chdir(simulation_path)
        cmd = [self.qsub_path, submission_script_path]
        subprocess_result = subprocess.run(
            cmd, 
            stdout=subprocess.PIPE)
        subprocess_stdout = subprocess_result.stdout.decode('utf-8')
        
        # return to the original context
        os.chdir(initial_path)

        jobid = subprocess_stdout.split('.')[0]
        # return the results
        return jobid


    def get_job_info(
        self,
        jobid='all', 
        username=None
    ):
        if username is None:
            username = os.environ['USER']
        else:
            username_ = username

        jobs_info = self.get_job_info_for_all_jobs(username=username_)
        return jobs_info

    def get_job_info_for_all_jobs(
        self,
        username: str
    ):
        cmd = ['qstat', '-u', username]
        subprocess_result = subprocess.run(
            cmd, 
            stdout=subprocess.PIPE
        )
        subprocess_stdout = subprocess_result.stdout.decode('utf-8')
        if subprocess_stdout == "":
            jobs_info = []
        else:
            jobs_info = {}

            lines = subprocess_stdout.split('\n')
            n_lines = len(lines)

            for i_line in range(2, n_lines):
                job_info = self._process_qstat_line(
                    line=lines[i_line]
                )

                jobs_info[job_info['jobid']] = deepcopy(job_info)

        return jobs_info

    def _process_qstat_line(self, line: str):

            tokens = [k.strip() for k in line.split()]
            for i, v in enumerate(tokens):
                print(i,v)

            # process tokens
            jobid = tokens[0].split('.')[0]
            username = tokens[2]
            jobname = tokens[3]
            jobstatus = self._convert_jobstatus_to_string(
                char_jobstatus=tokens[9]
            )

            #add information about this job
            job_info = {
                'jobid':jobid,
                'username':username,
                'jobname':jobname,
                'jobstatus':jobstatus
                }
            return job_info

    def _convert_jobstatus_to_string(self, char_jobstatus):
        char_to_jobstatus = {
            'C': 'complete',
            'H': 'held',
            'E': 'exiting',
            'Q': 'queued',
            'R': 'running',
            'T': 'moving',
            'W': 'waiting'
        }
        str_jobstatus = char_to_jobstatus[char_jobstatus]
        return str_jobstatus
        
