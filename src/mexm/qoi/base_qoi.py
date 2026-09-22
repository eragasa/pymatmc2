from copy import deepcopy

class Qoi(object):
    qoi_type = 'base_qoi'
    calculated_qoi_names = []
    is_base_class = True

    """ Abstract Quantity of Interest

    Args:
        qoi_name(str)
        qoi_type(str)
        structure_names(list of str)

    Attributes:
        qoi_name(str): this is a unique identifier
        structures(list of str):
        simulations_definitions(dict):
             simulations_definition[task_name]
             simulations_definition[task_name][task_type]
             simulations_definition[task_name][structure]
             simulations_definition[task_name][requires]
    """
    def __init__(self, qoi_name, structures):
        assert isinstance(qoi_name,str)
        assert isinstance(structures,dict)
        assert all([isinstance(k, str) for k, v in structures.items()])
        assert all([isinstance(v, str) for k, v in structures.items()])

        self.qoi_name = qoi_name
        self.structures = deepcopy(structures)
        self.simulation_definitions = {}

    def determine_simulations(self):
        raise NotImplementedError

    def calculate_qoi(self):
        raise NotImplementedError

    def add_simulation(self,
                       simulation_id,
                       simulation_name,
                       simulation_type,
                       simulation_structure,
                       bulk_structure=None):

        self.simulation_definitions[simulation_id] = {
            'simulation_name':simulation_name,
            'simulation_type':simulation_type,
            'simulation_structure':simulation_structure,
            'bulk_structure':bulk_structure
        }

    def determine_required_simulations(self):
        raise NotImplementedError

    def calculate_qois(self, simulation_results):
        raise NotImplementedError
