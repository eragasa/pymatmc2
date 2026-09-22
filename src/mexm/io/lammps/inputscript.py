class LammpsScript():

    def __init__(self):
        pass

    def initialize_from_dict(self, obj_dict):
        self.configuration = obj_dict

    def typeid_section_to_string(self):
        typeid_fmt = "group {s} type {i}"

        return_str = ""
        for i, s in enumerate(self.potential.symbols):
            return_st += typeid_fmt.format(s=s, i=str(i)) + "\n"

        return return_str

    def charge_section_to_string(self):

        charge_method_map ={
            'no_charge':self.charge_section_no_charge_to_string(),
            'static_charge':self.charge_section_static_charge_to_string(),
            'qeq':self.charge_section_qeq_charge_to_string_to_string()
        }

        return_str = charge_method_map[cls.potential.charge_type]()

        return return_str

    def charge_section_no_charge_to_string(self):
        raise NotImplementedError

    def charge_section_static_charge_to_string(self):
        raise NotImplementedError

    def charge_section_qeq_charge_to_string_to_string(self):
        qeq_compute_atom_charge_str = "compute charge{s} {s} property/atom q"
        qeq_compute_average_charge_str = "compute q{s} {s} reduce ave c_charge{s}"

        return_str = ""
        for symbol in self.potential.symbols():
            return_str += qeq_compute_atom_charge_str.format(s=symbol) + "\n"

        for symbol in self.potential.symbols():
            return_str += qeq_compute_average_charge_str.format(s=symbol) + "\n"

        qeq_variable_qtotal = "qtotal equal {}".format(
                "+".join(['count({s})*c_q{s}'])
            ) + "\n"


    def charge_summation_section_to_string(self):
        charge_summation_method_map = {
            'ewald': self.charge_summation_secton_ewald_to_string()
        }

        summation_type = self.configuration['charge_summation']['charge_summation_type']
        return_str = charge_summation_method_map[self.charge_summation]

    def charge_summation_section_ewald_to_string(self):
        kwargs = {
            'lammps_fix_id':lammps_fix_id,
            'lammps_group_id':lammps_group_id,
            'qfile_path':qeq_param_path
        }

        return_str = ""
        for potential in self.potentials.items():
            if isinstance(potential, Qeq):
                return_str = potential.lammps_fix_qeq_to_string(self, **kwargs)

        return_str
