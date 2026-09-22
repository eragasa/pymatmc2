from collections import OrderedDict
from mexm.qoi import ElasticProperties

class GulpElasticProperties(ElasticProperties):
    is_base_class = False
    qoi_type = 'gulp_elastic'
    calculated_qoi_names = [
        'gulp_{}'.format(k) for k in ElasticProperties.calculated_qoi_names
    ]
    ideal_simulation_type = 'gulp_elastic'

    def __init__(self, qoi_name, structures):
        super().__init__(qoi_name=qoi_name, structures=structures)
