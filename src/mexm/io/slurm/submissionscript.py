from collections import OrderedDict

from mexm.io.slurm import SlurmConfiguration

class SlurmSubmissionScript(object):
    """ class for creating slurm submission scripts """

    def __init__(self, configuration=None):

        if configuration is None:
            self.configuration = SlurmConfiguration()

        elif isinstance(configuration, SlurmConfiguration):
            self.configuration = configuration

        elif isinstance(configuration, OrderedDict):
            self.configuration = SlurmConfiguration.initialize_from_dict(
                    obj_dict=configuration
            )

        elif isinstance(configuration, str):

            if configuration == 'from_environment':
                path_ = os.environ('MEXM_SLURM_CONFIG_PATH')
                self.configuration = SlurmConfiguration(path=path_)

            # otherwise, assume configuration is a path variable
            else:
                path_ = configuration
                self.configuration = SlurmConfiguration(path=path_)

        else:
            raise TypeError()

        self.modules = None
        self.run_string = None

    @property
    def slurm_configuration(self):
        return self.configuration.slurm_configuration

    def write(self, path="runjob.slurm"):
        """ SLURM write slurm submission script

        Args:
            path(str): the path which to write SLURM submission script
        """

        with open(path, 'w') as f:
            f.write(self.slurm_script_to_string())

    def submit(self, path="runjob.slurm"):
        """ SLURM submit job to queue

        Args:
            path(str): the path which to write SLURM submission script
        """

        raise NotImplementedError

    def slurm_script_to_string(self):
        str_out = self.section_header_section_to_string()
        str_out += "#<---------- SLURM debug information\n"
        str_out += self.section_slurm_debug_to_string()
        str_out += '#<---------- load necessary modules\n'
        str_out += self.section_load_modules_to_string()
        str_out += '#<---------- run application\n'
        str_out += self.section_application_to_string()
        str_out += '#<---------- post-processing steps'
        str_out += self.section_postprocessing_to_string()
        return str_out

    def section_header_section_to_string(self):
        job_name_ = self.slurm_configuration['job_name']
        qos_ = self.slurm_configuration['qos']
        mail_type_ = self.slurm_configuration['mail_type']
        mail_name_ = self.slurm_configuration['mail_name']
        ntasks_ = self.slurm_configuration['ntasks']
        time_ = self.slurm_configuration['time']
        output_path_ = self.slurm_configuration['output_path']
        error_path_ = self.slurm_configuration['error_path']
        memory_ = self.slurm_configuration['memory']

        str_out = '#!/bin/bash\n'
        str_out += '#SBATCH --job-name={}\n'.format(job_name_)
        str_out += '#SBATCH --qos={}\n'.format(qos_)
        str_out += '#SBATCH --mail-type={}\n'.format(mail_type_)
        str_out += '#SBATCH --mail-user={}\n'.format(mail_name_)
        str_out += '#SBATCH --ntasks={}\n'.format(ntasks_)
        str_out += '#SBATCH --distribution=cyclic:cyclic\n'
        str_out += '#SBATCH --mem={}\n'.format(memory_)
        str_out += '#SBATCH --time={}\n'.format(time_)
        str_out += '#SBATCH --output={}\n'.format(output_path_)
        str_out += '#SBATCH --error={}\n'.format(error_path_)

        return str_out

    def section_application_to_string(self):

        str_out = None

        if self.run_string is None:
            str_out = "# no application string provided\n"
        else:
            str_out = "{}\n".format(run_string)
        return str_out

    def section_load_modules_to_string(self):

        str_out = None
        fmt_module_load = "module_load_{}\n"

        if self.modules is None:
            str_out = "# no modules loaded \n"
        else:
            str_out = ""
            for module in self.modules:
                str_out += fmt_module_load.format(module)

        return str_out

    def section_slurm_debug_to_string(self):
        str_out = "\n".join([
                'echo slurm_job_id:$SLURM_JOB_ID',
                'echo slurm_job_name:$SLURM_JOB_NAME',
                'echo slurm_job_nodelist:$SLURM_JOB_NODELIST',
                'echo slurm_job_num_nodes:$SLURM_JOB_NUM_NODES',
                'echo slurm_cpus_on_node:$SLURM_CPUS_ON_NODE',
                'echo slurm_ntasks:$SLURM_NTASKS',
                'echo working directory:$(pwd)',
                'echo hostname:$(hostname)',
                'echo start_time:$(date)'
            ]) + "\n"
        return str_out

    def section_postprocessing_to_string(self):
        str_out = "touch jobComplete\n"
        str_out += "echo end_time:$(date)\n"
        return str_out
