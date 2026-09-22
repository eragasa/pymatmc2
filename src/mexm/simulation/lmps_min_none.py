import os
from copy import deepcopy
from collections import OrderedDict
from mexm.simulation import LammpsSimulation, StaticCalculation

class LammpsStaticCalculation(LammpsSimulation, StaticCalculation):
    simulation_type = 'lammps_min_none'
    is_base_class = False

    """ Class for LAMMPS structural minimization

    This data class defines additional attributes and methods necessary to
    interact with the Workflow manager.

    Args:
        name (str)
        simulation_path (str)
        structure_path (str)
        bulk_structure_name (str)

    Attributes:
        config
        config_map
    """
    def __init__(self,
                 name,
                 simulation_path,
                 structure_path,
                 bulk_structure_name=None):

        LammpsSimulation.__init__(self,
                                  name=name,
                                  path=name,
                                  structure_path=structure_path,
                                  bulk_structure_name=bulk_structure_name)
    
        self.update_status()
    
    def lammps_input_file_to_string(self):
        str_out = "".join([\
                self._lammps_input_initialization_section(),
                self._lammps_input_create_atoms(),
                self._lammps_input_define_potential(),
                self._lammps_input_run_minimization(),
                self._lammps_input_out_section()])
        return(str_out)

    def get_conditions_config(self):
        config_tests = {
            'potential_initialized':self.is_potential_initialized,
            'parameters_processed':self.is_potential_parameters_processed,
            'bulk_structure_results':self.is_bulk_structure_results
        }

        self.conditions_CONFIG = OrderedDict(
            [(k,v())for k,v in config_tests.items()]
        )
        
        return self.conditions_READY

    def is_bulk_structure_results(self):
        if self.bulk_structure_name is None:
            return True
        else:
            a0_name = '{}.{}.a0'.format(self.bulk_structure_name, 'lmps_min_all')
            if self.results is None:
                return False
            else:
                return a0_name in self.results

    def on_init(self,configuration=None,results=None):
        LammpsSimulation.on_init(self,configuration=configuration)

    def on_config(self,configuration,results=None):
        LammpsSimulation.on_config(self,configuration=None,results=results)

    def on_ready(self,configuration, results=None):
        self.results = deepcopy(results)
        LammpsSimulation.on_ready(self,
                                  configuration=configuration,
                                  results=self.results)

    def modify_structure_file(self, results):
        if self.bulk_structure_name is not None:
            a0_name = '{}.{}.a0'.format(self.bulk_structure_name, 'lmps_min_all')
            self.lammps_structure.a0 = results[a0_name]

    def on_post(self,configuration=None):
        self.__get_results_from_lammps_outputfile()
        LammpsSimulation.on_post(self,configuration=configuration)

    def __get_results_from_lammps_outputfile(self):
        LammpsSimulation.__get_results_from_lammps_outputfile(self)

    def _lammps_input_run_minimization(self):
        str_out = (
            '# ---- define settings\n'
            'compute eng all pe/atom\n'
            'compute eatoms all reduce sum c_eng\n'
            '# ---- run minimization\n'
            'reset_timestep 0\n'
            'thermo 10\n'
            'thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n'
            'run 0\n'
            )
        return str_out
