
class HpcCluster():
    """ this class modules an HPC cluster """

    def __init__(self):
        pass

    def submit_job(self, job_script_path):
        raise NotImplementedError

    def delete_job(self, job_id):
        raise NotImplementedError

    def get_job_status_all(self):
        raise NotImplementedError

    def get_job_stats_by_job(self, job_id):
        raise NotImplementedError



 
class SlurmCluster(HpcCluster):

class HpcJob():
    def __init__(self):
        self.job_name = None
        self.account = None
        self.queue = None
        self.qos = None
        self.wall_time_limit = None
        self.node_count = None
        self.process_count_per_node = None
        self.memory_limit= = None
        self.memory_limit_per_processor = None
        self.dst_output

class TorqueJob():

    def __init__(self):
        HpcJob.__init__(self)

    def job_name_to_str(self):
        name = self.job_name
        return "#PBS -N {}".format(name)


