import yaml
from collections import OrderedDict
from mexm.io.filesystem import OrderedDictYAMLLoader
# from mexm.exceptions import MexmDesignConstraintException
# from mexm.exceptions import MexmObjectiveConstraintException

design_variable_constraint_types = ['=','<','>']
objective_variable_constraint_types = ['=','<','>']
objective_variable_optimization_type = ['min', 'max']

class MultiObjectiveOptimization(object):
    """ base multi objective optimization object

    This object defines the optimization problem which can then be subclassed
    to define a specific optimization problem.  The purpose of this implementation
    is to separate the definition of the problem, from the optimizer that
    solves the problem.

    """

    def __init__(self):
        self.design_variables = OrderedDict()
        self.additional_design_constraints = []

        self.objective_variables = OrderedDict()
        self.additional_objective_constraints = []

    def write(self, path):
        assert isinstance(path, str)

        obj_dict = None
        if str.endswith('.yaml'):
            with open(path, 'w') as f:
                yaml.dump(obj_dict, f, default_flow_style=False)
                f.write(str_out)
        else:
            msg = "unknown file extention"
            raise ValueError(msg)

    def read(self, path):
        assert isinstance(path, str)

        if str.endswith('.yaml'):
            with open(path, OrderedDictYAMLLoader) as f:
                obj_dict = yaml.load(f,OrderedDictYAMLLoader)
                self.configure_from_dict(obj_dict=obj_dict)
        else:
            msg = "unknown file extension"
            raise ValueError(msg)

    def evaluate(self, design_variables):
        raise NotImplementedError

    def __configure_from_dict(self, obj_dict):
        self.design_variables = obj_dict['design_variables']
        self.objective_variables = obj_dict['objective_variables']
        self.additional_design_constraints \
            = obj_dict['additional_design_constraints']
        self.additional_objective_constraints \
            = obj_dict['additional_objective_constraints']

    def initialize_from_dict(obj_dict):
        """ static factory method to initialize from a dictionary object """
        assert isinstance(obj_dict, dict)

        moo = MultiObjectiveOptimization()
        self.__configure_from_dict(obj_dict=obj_dict)
        return moo

    def initialize_from_yaml_file(path):
        """ static factory method to initialize from a yaml file """

        moo = MultiObjectiveOptimization()
        moo.read(path=path)
        return moo

    def to_dict(self):
        """ deserialize object to dictionary """
        obj_dict = OrdereDict()
        obj_dict['design_variables'] = self.design_variables
        obj_dict['additional_design_constraints'] \
            = self.additional_design_constraints
        obj_dict['objective_variables'] = self.objective_variables
        obj_dict['additional_objective_constraints'] \
            = self.additional_objective_constraints
        return obj_dict

    def to_yaml_string(self):
        """ deserialize object to yaml string """

        obj_dict = self.to_dict()
        yaml_str = ""
        raise NotImplementedError()
        return yaml_str

    def _create_design_variable(self, name):
        self.design_variables[name] = OrderedDict()
        self.design_variables[name]['constraints'] = []

    def _create_objective_variable(self, name):
        self.objective_variables[name] = OrderedDict()
        self.objective_variables[name]['optimization_type'] = None
        self.objective_variables[name]['constraints'] = []

    @property
    def design_variable_names(self):
        return [k for k in self.design_variables.items()]

    @property
    def objective_variable_names(self):
        return [k for k in self.objective_variables]

    @property
    def design_constraints(self):
        design_constraints = [
            [k, v['constraints'][0], v['constraints'][1]]
            for k,v in self.design_variables.items()
        ]

        desgin_constraints += self.additional_design_constraints
        return objective_contraints

    @property
    def objective_constraints(self):
        objective_constraints = [
            [k, v['constraints'][0], v['constraints'][1]]
            for k,v in self.objective_variables.items()
        ]

        objective_constraints += self.additional_objective_constraints
        return objective_constraints

    def add_design_variable(self, name, constraints=None):
        assert isinstance(name, str)
        assert isinstance(constraints, list) or constraints is None

        self._create_design_variable(name)
        if constraints is not None:
            for k in constraints:
                constraint_type = k[0]
                assert constraint_type in design_variable_constraint_types

                try:
                    constraint_value = float(k[1])
                except ValueError as e:
                    if isinstance(k[1],str):
                        constraint_value = k[1]
                    else:
                        raise

            self.design_variables[name]['constraints'].append(
                [constraint_type, constraint_value]
            )

    def add_design_constraint(self,
                              constraint_string,
                              constraint_type,
                              constraint_value):
        """ add a design constraint


        Args:
            constraint_string(str): should be a string consisting of design
                variables plus any arithmetic operator


        """
        if constraint_string in self.design_variables:
            self.constraint_variables[constraint_string]['constraints'] = [
                constraint_type,
                constraint_value
            ]
        else:
            self.additional_design_constraints.append([
                constraint_string,
                constraint_type,
                constraint_value
            ])

    def add_objective_variable(self,
                               name,
                               optimization_type,
                               constraints=None):
        """ add objective variable

        Args:
            name(str): name
            optimization_type(str): either min or max.
            constraints(list)
        """
        assert isinstance(name, str)
        assert optimization_type in objective_variable_optimization_type
        assert any([
            isinstance(constraints, list),
            constraints is None
        ])

        self._create_objective_variable(name)

        self.objective_variables[name]['optimization_type'] = optimization_type
        if constraints is not None:
            for k in constraints:
                constraint_type = k[0]
                assert constraint_type in design_variable_constraint_types

                # (1) try to cast into a float
                # (2) if this fails, then it might be a string
                # (3) everything else is an invalid options
                try:
                    constraint_value = float(k[1])
                except ValueError as e:
                    if isinstance(k[1],str):
                        constraint_value = k[1]
                    else:
                        raise

            self.objective_variables[name]['constraints'].append(
                [constraint_type, constraint_value]
            )

    def design_variables_to_string(self):

        return_str = ""
        for k,v in self.design_variables.items():
            variable_name = k
            variable_constraints = v['constraints']

            for i, constraint in enumerate(variable_constraints):
                constraint_type = constraint[0]
                constraint_value = constraint[1]
                if i == 0:
                    return_str += "{:^20}{:^5}{:<60}\n".format(
                            k,
                            constraint_type,
                            constraint_value
                    )
                else:
                    return_str += "{:^20}{:^5}{:<60}\n".format(
                            "",
                            constraint_type,
                            constraint_value
                    )
        return return_str

    def objective_variables_to_string(self):

        return_str = ""
        for k,v in self.objective_variables.items():
            variable_name = k
            variable_type = v['optimization_type']
            variable_constraints = v['constraints']

            for i, constraint in enumerate(variable_constraints):
                constraint_type = constraint[0]
                constraint_value = constraint[1]
                if i == 0:
                    return_str += "{:^20}{:^5}{:^7}{:<60}\n".format(
                            k,
                            variable_type,
                            constraint_type,
                            constraint_value
                    )
                else:
                    return_str += "{:^20}{:^5}{:^7}{:<60}\n".format(
                            "",
                            variable_type,
                            constraint_type,
                            constraint_value
                    )
        return return_str

    def __str__(self):
        return_str = ""
        return_str += 80*"-" + "\n"
        return_str += "{:^80}\n".format("DESIGN VARIABLES")
        return_str += 80*"-" + "\n"
        return_str += self.design_variables_to_string()
        return_str += 80*"-" + "\n"
        return_str += "{:^80}\n".format("OBJECTIVE VARIABLES")
        return_str += 80*"-" + "\n"
        return_str += self.objective_variables_to_string()

        return return_str
