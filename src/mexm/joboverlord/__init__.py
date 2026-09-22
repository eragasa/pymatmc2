from mexm.simulation import AtomicSimulation
from mexm.simulation import VaspSimulation
from mexm.simulation import LammpsSimulation
from mexm.io.vasp import Incar

from mexm.io.torque import TorqueConfiguration
from mexm.io.torque import TorqueSubmissionScript
from mexm.io.slurm import SlurmConfiguration
from mexm.io.slurm import SlurmSubmissionScript

class JobOverlordConfiguration(object):
    supported_cluster_types = {
        'torque':TorqueConfiguration,
        'slurm':SlurmConfiguration
    }
    """

    Attributes:
        mexm_sqlite_path(str)
    """

    def __init__(self):
        try:
            self.mexm_sqlite_path = os.environ['MEXM_SQLITE_PATH']
        except KeyError:
            self.mexm._sqlite_path = None

        self.clusters = {}

    def to_dict(self):
        config_dict = {}
        config_dict['clusters'] = self.clusters
        return config_dict

    def read(self):
        pass

    def write(self):
        pass

    def add_cluster_configuration(self, cluster_name, cluster_configuration):
        if not isinstance(cluster_name, str):
            raise TypeError('cluster_name must be a string')

        if isinstance(cluster_configuration, TorqueConfiguration):
            self.clusters[cluster_name] = cluster_configuration.to_dict()
        elif isinstance(cluster_configuration, SlurmConfiguration):
            self.clusters[cluster_name] = cluster_configuration.to_dict()
        elif isinstance(cluster_configuration, str):
            obj_config = ClusterConfiguration()
            obj_config.read(path=cluster_configuration)
            cluster_type = obj_config.cluster_type
            self.clusters[cluster_name] = self.supported_cluster_types[cluster_type]()
            self.clusters[cluster_name].read(cluster_configuration)
        else:
            msg = (
                'cluster_configuration must either be a path string or '
                'a ClusterConfiguration'
            )
            raise TypeError(msg)


class JobOverlord(object):

    def __init__(self):
        self.configuration_ = JobOverlordConfiguration()

    @property
    def configuration(self):
        return self.configuration_

    @configuration.setter
    def configuration(self, configuration):
        if not isinstance(JobOverlordConfiguration):
            raise TypeError

    def register_job(simulation):
        if not isinstance(simulation, AtomicSimulation):
            raise TypeError('simulation must be an AtomicSimulation')

        if isinstance(simulation, VaspSimulation):
            self.register_vasp_job(simulation)
        elif isinstance(simulation, LammpsSimulation):
            self.register_lammps_job(simulation)

    def register_vasp_job(simulation, n_cores, job_name):
        if not ininstance(simulation, VaspSimulation):
            raise TypeError('simulation must be a VaspSimulation')

        simulation.incar.npar = Incar.determine_npar(n_cores=n_cores)

    def register_lammps_job(simulation, n_cores, job_name):
        if not instance(simulation, LammpsSimulation):
            raise TypeError('simulation must be a LammpsSimulation')