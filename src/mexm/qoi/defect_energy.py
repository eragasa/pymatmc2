from collections import OrderedDict
from mexm.qoi import Qoi

class DefectFormationEnergy(Qoi):
    qoi_type = 'base_defect'
    calculated_qoi_names = ['E_formation']
    is_base_class = True

    ideal_simulation_type = 'min_all'
    defect_simulation_type = 'min_pos'

    def __init__(self, qoi_name, structures):
        assert isinstance(qoi_name, str)
        assert isinstance(structures, dict)

        Qoi.__init__(self, qoi_name=qoi_name, structures=structures)

    def determine_simulations(self):
        ideal_structure_name = self.structures['ideal']
        ideal_simulation_type = self.ideal_simulation_type
        ideal_simulation_name = '{}.{}'.format(
            ideal_structure_name,
            ideal_simulation_type
        )
        self.add_simulation(
            simulation_id='ideal',
            simulation_name=ideal_simulation_name,
            simulation_type=ideal_simulation_type,
            simulation_structure=ideal_structure_name
        )

        defect_structure_name = self.structures['defect']
        defect_simulation_type = self.defect_simulation_type
        defect_simulation_name = '{}.{}'.format(
            defect_structure_name,
            defect_simulation_type
        )
        self.add_simulation(
            simulation_id='defect',
            simulation_name=defect_simulation_name,
            simulation_type=defect_simulation_type,
            simulation_structure=defect_structure_name,
            bulk_structure=ideal_structure_name
        )

    def calculate_qois(self, results):
        assert isinstance(results, dict)

        ideal_prefix = self.simulation_definitions['ideal']['simulation_name']
        defect_prefix = self.simulation_definition['defect']['simulation_name']

        ideal_ecoh = results['{}.{}'.format(ideal_prefix,'toten')]
        ideal_natoms = results['{}.{}'.format(ideal_prefix,'natoms')]

        defect_ecoh = results['{}.{}'.format(defect_prefix,'toten')]
        defect_natoms = results['{}.{}'.format(ideal_prefix,'natoms')]
