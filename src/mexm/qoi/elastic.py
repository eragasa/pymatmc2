from mexm.qoi import Qoi

class ElasticProperties(Qoi):
    qoi_type = 'base_elastic'
    calculated_qoi_names = [
        'c11','c12','c13'
        ,'c22','c33','c44','c55','c66',
        'bulk', 'shear']
    is_base_class = False
    ideal_simulation_type = 'base_elastic'

    def __init__(self, qoi_name, structures):
        Qoi.__init__(self,
                     qoi_name=qoi_name,
                     structures=structures)

    @staticmethod
    def calculate_bulk_modulus(c11, c12, c13, c22, c33, c44, c55, c66):

        bulk = (c11+2.*c12)/3.
        return bulk

    @staticmethod
    def calculate_shear_modulus(c11, c12, c13, c22, c33, c44, c55, c66):
        shear = (c11-c12)/2
        return shear

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

    def calculate_qois(self, results):
        ideal_prefix = self.simulation_definitions['ideal']['simulation_name']
        toten = results['{}.{}'.format(ideal_prefix,'toten')]
        natoms = results['{}.{}'.format(ideal_prefix,'natoms')]
        totpress = results['{}.{}'.format(ideal_prefix,'totpress')]
