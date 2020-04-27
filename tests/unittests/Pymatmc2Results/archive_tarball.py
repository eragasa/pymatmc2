import tarfile
import os
from pymatmc2 import Pymatmc2Configuration


class Pymatmc2Tarball():
    def __init__(self, dst_path):
        self._configuration = None
        self.tar_path = dst_path
        self.results_path = 'results'
        self.resource_path = 'resources'
    
    @property
    def configuration(self) -> Pymatmc2Configuration:
        return self._configuration

    @configuration.setter
    def configuration(self, configuration: Pymatmc2Configuration):
        self._configuration = configuration


    def add_resource_dir(self):
        with tarfile.open(self.tar_path, 'w:') as tar:
            src_path = self.resource_path
            tar.add(src_path, arcname = os.path.basename(src_path))

    def get_phase_point_name(self) -> str:
        T = int(self.configuration.temperature)
        P = int(self.configuration.pressure)

        fmt = '{}K_{}GPa'
        return fmt.format(T, P)

    def add_iteration_results_dir(self, i: int):
        src_path = os.path.join(
            self.results_path,
            self.get_phase_point_name(),
            '{:05}'.format(i)
        )
        with tarfile.open(self.tar_path, 'a:') as tar:
            tar.add(src_path, arcname = os.path.basename(src_path))

if __name__ == "__main__":
    
    configuration_path = 'pymatmc2.config'
    configuration = Pymatmc2Configuration()
    configuration.read(configuration_path)
    
    tar_path = configuration.results_tar_path
    print(tar_path)

    tarball = Pymatmc2Tarball(dst_path=tar_path)
    tarball.configuration = configuration
    tarball.add_resource_dir()
    tarball.add_iteration_results_dir(i=1) 
