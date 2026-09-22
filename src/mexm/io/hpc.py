class HpcClusterConfiguration():

    def __init__(self):
        self.path_ = None
        self.account_ = None
        self.queue_ = None
        self.qos_ = None
        self.wall_time_limit = None
        self.node_count = None

    @property
    def path(self):
        return self.path_
    @path.setter
    def path(self, path):
        if not isinstance(path, str):
            raise TypeError('path must be a string')
        self.path_ = path

class HpcCluster():
    """ this class models an hpc cluster """

    def __init__(self):
        pass

    def submit_job(self, job_script_path):
        raise NotImplementedError

    def delete_job(self, job_id):
        raise NotImplementedError

    def get_job_status_all(self):
        raise NotImplementedError

    def get_job_status_by_job(self):
        raise NotImplementedError


class HpcJob():
    """ this class models an hpc job """

    def __init__(self):
        pass

class SlurmJob():

    def __init__(self):
        self.job_name = None
        self.account = None


