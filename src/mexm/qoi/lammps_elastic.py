from collections import OrderedDict
from mexm.qoi import ElasticProperties

class LammpsElasticProperties(ElasticProperties):
    qoi_type = 'lammps_elastic'
    calculated_qoi_names = [
        'lammps_{}'.format(k) for k in ElasticProperties.calculated_qoi_names
    ]
    ideal_simulation_type = 'lmps_elastic'

    def __init__(self,qoi_name,structures):
        assert isinstance(qoi_name, str)
        assert isinstance(structures, dict)

        ElasticProperties.__init__(self,
                                   qoi_name=qoi_name,
                                   structures=structures)

    def determine_simulations(self):
        super().determine_simulations()

    def calculate_qois(self,simulation_results):
        ideal_sim_name = self.simulation_definitions['ideal']

        sim_result_names = ['c11','c12','c13','c22','c33','c44','c55','c66']
        for i, k in enumerate(qoi_names):
            result_name = '{}.{}'.format(ideal_sim_name,
                                         k)
            qoi_name = '{}.{}'.format(ideal_sim_name,
                                      LammpsElasticProperties.qois_calculated[i])
            self.qois[qoi_name] = simulation_results[result_name]

        c11 = simulation_results['{}.{}'.format(ideal_sim_name,'c11')]
        c12 = simulation_results['{}.{}'.format(ideal_sim_name,'c12')]
        # c13 = simulation_results['{}.{}'.format(_prefix,'c13')]
        # c22 = simulation_results['{}.{}'.format(_prefix,'c22')]
        # c33 = simulation_results['{}.{}'.format(_prefix,'c33')]
        # c44 = simulation_results['{}.{}'.format(_prefix,'c44')]
        # c55 = simulation_results['{}.{}'.format(_prefix,'c55')]
        # c66 = simulation_results['{}.{}'.format(_prefix,'c66')]

        # References:
        # http://homepages.engineering.auckland.ac.nz/~pkel015/SolidMechanicsBooks/Part_I/BookSM_Part_I/06_LinearElasticity/06_Linear_Elasticity_03_Anisotropy.pdf

        bulk = (c11+2*c12)/3
        shear = (c11-c12)/2

        self.qois['{}.bulk_lammps'.format(ideal_sim_name)] = bulk
        self.qois['{}.shear_lammps'.format(ideal_sim_name)] = shear
