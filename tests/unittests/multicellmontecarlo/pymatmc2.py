import os
import argparse
from crontab import CronTab

class JobSubmissionSystemOverlord():
    def __init__(self):
        self.username = None
        self.maximum_cores = None
        self.maximum_nodes = None
        self.maximum_cores_per_node = None
        self.maximum_memory_per_node = None
        self.max_cores_per_job = None
        self.max_nodes_per_job = None

    def create_new_project(self, 
                           project_name,
                           project_description):
        pass

    def create_new_job(self,
                       job_name,
                       project_name,
                       project_description):
        pass

    def submit_job(self,
                   job_name,
                   runjob_script='runjob'):
        context_path = os.get_cwd()
        submission_path = self.get_job_path(job_name)

        # switch to the submission path
        os.chdir(submission_path)
        
        cmd  = "qsub {}".format(runjob_script)
    def get_open_jobs(self, project_name):
        pass

    def get_finished_jobs(self, project_name):
        pass

class TorqueSubmissionFile():
    def __init__(self):
        self.src_path = None
        self.dst_path = None
        self.job_name

    def read(self, src_path):
        self.src_path = src_path

    def write(self):
        self.dst_path = dst_path

    def copy(self, src_path, dst_path):
        self.src_path = src_path
        self.dst_path = dst_path

    def header_section(self):
        walltime = "4:00:00"
        jobname = self.job_name
        n_nodes = self.n_nodes
        str_out = "\n".join([
            "#!/bin/bash",
            "PBS -A PAA0028",
            "PBS -l walltime={}".format(walltime),
            "PBS -l nodes={}:ppn={}".format(n_nodes,ppn),
            "PBS -l {}".format(jobname)
        ])

    def load_modules(self):
        modules = [
            'module load intel',
            'module load intelmpi',
            'module load python'
        ]
        for module in modules:
            str_out += "module load {}\n".format(module)

    def execution_lines(self):
        pass

class TorqueSystemOverlord(JobSubmissionSystemOverlord):
    pass

class SlurmSystemOverlord(JobSubmissionSystemOverlord):
    pass

def overlord_turn_on(job_manager_type):
    supported_job_manager_types = ['torque', 'slurm']
    if job_manager_type not in supported_job_manager_types:
        msg = "unsupported job_manager_type: {}".format(job_manager_type)
        raise ValueError(msg)
    python_bin = "/opt/anaconda/bin/python"
    mexm_sbin = os.path.realpath(__file__)
    cronjob_cmd = "{} {} --overlord cronjob --job_manager_type {}".format(
            python_bin, 
            mexm_sbin,
            job_manager_type)
    cronjob_cmt = "mexm job manager"

    cron = CronTab(user=True)
    job = cron.new(command=cronjob_cmd,
             comment=cronjob_cmt)
    job.minute.every(1)
    cron.write_to_user(user=True)

def overlord_turn_off():
    supported_job_manager_types = ['torque', 'slurm']
    
    python_bin = "/opt/anaconda/bin/python"
    mexm_sbin = os.path.realpath(__file__)
    cronjob_cmd = "{} {} --overlord cronjob".format(python_bin, mexm_sbin)
    cronjob_cmt = "mexm job manager"

    cron = CronTab(user=True)
    job = cron.find_comment(cronjob_cmt)
    cron.remove(job)
    cron.write_to_user(user=True)


def overlord_cronjob():
    mexm_sbin = os.path.realpath(__file__)
    print('overlord cronjob')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--overlord', 
                        help = "job submission")
    parser.add_argument('--job_manager_type',
                        help = "job submission manager type")
    args = parser.parse_args()

    job_manager_type_options = ['torque', 'slurm']
    if args.job_manager not in job_manager_options:
        msg = "job_manager must either be torque or slurm"
        print(msg)
        raise ValueError(msg)
    job_manger_type = args.job_manager

    overlord_options = ['on', 'off']
    if args.overlord == 'on':
        overlord_turn_on()
    elif args.overlord == 'off':
        overlord_turn_off()
    elif args.overlord == 'cronjob':
        overlord_cronjob()
    else:
        msg = 'the supported overlord options are {}'.format(','.join(overlord_options))
        print(msg)
