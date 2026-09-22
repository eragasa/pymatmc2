from mexm.qoi import StaticStructure

class VaspStaticStructure(StaticStructure):
    qoi_type = 'vasp_min_none'
    is_base_class = False
    qois_calculated = [
        'vasp_Ecoh_min_none',
        'vasp_a1_min_none',
        'vasp_a2_min_none',
        'vasp_a3_min_none',
        'vasp_a11_min_none',
        'vasp_a12_min_none',
        'vasp_a13_min_none',
        'vasp_a21_min_none',
        'vasp_a22_min_none',
        'vasp_a23_min_none',
        'vasp_a31_min_none',
        'vasp_a32_min_none',
        'vasp_a33_min_none',
        'vasp_p11_min_none',
        'vasp_p12_min_none',
        'vasp_p13_min_none',
        'vasp_p21_min_none',
        'vasp_p22_min_none',
        'vasp_p23_min_none',
        'vasp_p31_min_none',
        'vasp_p32_min_none',
        'vasp_p33_min_none'
    ]

    def __init__(self, qoi_name, structures):
        assert isinstance(qoi_name,str)
        assert isinstance(structures, dict)

        StaticStructure.__init__(self, qoi_name=qoi_name, structures=structures)

    def determine_simulations(self):
        ideal_structure_name = self.structure['ideal']
        ideal_simulation_type = 'vasp_min_none'
        ideal_simulation_name = '{}.{}'.format(ideal_structure_name,
                                               ideal_simulation_type)

        self.add_simulation(sim_id='ideal',
                            simulation_name=ideal_simulation_name,
                            simulation_type=ideal_simulation_type,
                            simulation_structure=ideal_structure_name)

    def calculate_qois(self, simulation_results):
        raise NotImplementedError
