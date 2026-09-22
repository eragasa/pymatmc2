from copy import deepcopy
from collections import OrderedDict
from mexm.io.filesystem import OrderedDictYAMLLoader

class SlurmConfigurationException(Exception): pass

class SlurmConfiguration():
    """
    Args:
        path(str): path to the YAML formatted configuration filew
    Raises:
        KeyError: if the MEXM_SLURM_CONFIG_PATH is not setWhen the path argument is not specified then, the path
            will be retrieved from MEXM_SLURM_CONIG_PATH,
    """

    def __init__(self, path=None):
        self.configuration = OrderedDict()

        if path is not None:
            self.read(path=path)

    @staticmethod
    def initialize_from_dict(obj_dict):
        """ factory method which initializes and object from a dict """
        obj_slurm = SlurmConfiguration()
        obj_slurm.configure_from_dict(obj_dict=obj_dict)
        return obj_slurm

    def read(self, path):
        with open(path, 'r') as f:
            obj_dict = yaml.load(f, OrderedDictYAMLLoader)
            self.configure_from_dict(obj_dict=obj_dict)

    def write(self,path):
        with open(path, 'w') as f:
            yaml_str = yaml.dump(self.to_dict(),
                                 default_flow_style=False
                                 )
            f.write(yaml_str)

    def to_dict(self):
        return self.configuration

    def configure_from_dict(self, obj_dict):
        assert isinstance(obj_dict, dict)

        self.configuration = deepcopy(obj_dict)

    @property
    def slurm_configuration(self):
        try:
            slurm_configuration = self.configuration['slurm']
        except KeyError:
            msg = "SLURM configuration has not been provided"
            raise SlurmConfigurationException(msg)
        return self.configuration['slurm']

    @slurm_configuration.setter
    def slurm_configuration(self, configuration):
        assert isinstance(configuration, OrderedDict)
        self.configuration['slurm'] = deepcopy(configuration)

    @property
    def lammps_configuration(self):
        return self.configuration['lammps']

    @lammps_configuration.setter
    def lammps_configuration(self, configuration):
        assert isinstance(configuration, OrderedDict)
        self.configuration['lammps'] = deepcopy(configuration)

    @property
    def vasp_configuration(self):
        return self.configuration['vasp']

    @vasp_configuration.setter
    def vasp_configuration(self, configuration):
        assert isinstance(configuration, OrderedDict)
        self.configuration['vasp'] = deepcopy(configuration)

    @property
    def phonts_configuration(self):
        return self.configuration['phonts']

    @phonts_configuration.setter
    def phonts_configuration(self, configuration):
        assert isinstance(configuration, OrderedDict)
        self.configuration['phonts'] = deepcopy(configuration)

    def set_slurm_configuration(self,
                                job_name,
                                qos,
                                mail_type,
                                mail_name,
                                ntasks,
                                output_path,
                                error_path,
                                time,
                                memory):
        self.configuration['slurm'] = OrderedDict([
            ('job_name',job_name),
            ('qos', qos),
            ('mail_type', mail_type),
            ('mail_name', mail_name),
            ('ntasks', ntasks),
            ('output_path', output_path),
            ('error_path', error_path),
            ('time', time),
            ('memory', memory)
        ])

    def set_vasp_configuration(self,
                               modules,
                               run_string):
        self.configuration['vasp'] = OrderedDict([
            ('modules', deepcopy(modules)),
            ('run_string', run_string)
        ])

    def set_phonts_configuration(self,
                                 modules,
                                 run_string):
        self.configuration['phonts'] = OrderedDict([
            ('modules', deepcopy(modules)),
            ('run_string', run_string)
        ])

    def set_lammps_configuration(self,
                                 modules,
                                 run_string):
        self.configuration['lammps'] = OrderedDict([
            ('modules', deepcopy(modules)),
            ('run_string', run_string)
        ])
