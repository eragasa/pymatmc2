from mexm.qoi import StructuralRelaxation

class LammpsStructuralRelaxation(StructuralRelaxation):
    qoi_type = 'lammps_min_all'
    is_base_class = False
    qois_calculated = [
        'lammps_Ecoh_min_all',
        'lammps_a1_min_all',
        'lammps_a2_min_all',
        'lammps_a3_min_all',
        'lammps_a11_min_all',
        'lammps_a12_min_all',
        'lammps_a13_min_all',
        'lammps_a21_min_all',
        'lammps_a22_min_all',
        'lammps_a23_min_all',
        'lammps_a31_min_all',
        'lammps_a32_min_all',
        'lammps_a33_min_all',
        'lammps_p11_min_all',
        'lammps_p12_min_all',
        'lammps_p13_min_all',
        'lammps_p21_min_all',
        'lammps_p22_min_all',
        'lammps_p23_min_all',
        'lammps_p31_min_all',
        'lammps_p32_min_all',
        'lammps_p33_min_all'
    ]

    def __init__(self, qoi_name, structures):
        assert isinstance(qoi_name, str)
        assert isinstance(structures, dict)
        StructuralRelaxation.__init__(self,
                                      qoi_name=qoi_name,
                                      structures=_structures)

    def determine_simulations(self):
        ideal_structure_name = self.structures['ideal']
        ideal_simulation_type = 'lmps_min_all'
        ideal_simulation_name = "{}.{}".format(
                _structure_ideal_name,
                _task_type)
        _task_requires = None
        self.add_simulation(
                simulation_id='ideal',
                simulation_name=ideal_simulation_name,
                simulation_type=ideal_simulation_type,
                simulation_structure=ideal_structure_name)

    def calculate_qois(self,simulation_results):
        ideal_sim_name = self.simulation_definitions['ideal']['simulation_name']

        e_min_none = simulation_results['{}.{}'.format(ideal_sim_name,'toten')]
        n_atoms = simulation_results['{}.{}'.format(ideal_sim_name,'natoms')]
        ecoh = e_min_none / n_atoms

        p_tot = simulation_results['{}.{}'.format(ideal_sim_name,'totpress')]

        # getting the individual components of the H-matrix
        # H = [ [a11,a12,a13]
        #       [a21,a22,a33]
        #       [a31,a32,a33]]

        a11 = simulation_results['{}.{}'.format(ideal_sim_name,'a11')]
        a12 = simulation_results['{}.{}'.format(ideal_sim_name,'a12')]
        a13 = simulation_results['{}.{}'.format(ideal_sim_name,'a13')]
        a21 = simulation_results['{}.{}'.format(ideal_sim_name,'a21')]
        a22 = simulation_results['{}.{}'.format(ideal_sim_name,'a22')]
        a23 = simulation_results['{}.{}'.format(ideal_sim_name,'a23')]
        a31 = simulation_results['{}.{}'.format(ideal_sim_name,'a31')]
        a32 = simulation_results['{}.{}'.format(ideal_sim_name,'a32')]
        a33 = simulation_results['{}.{}'.format(ideal_sim_name,'a33')]

        # calculate the length of cells
        # a = a1
        # b = a2
        # c = a3
        a1 = (a11**2+a12**2+a13**2)**0.5
        a2 = (a21**2+a22**2+a23**2)**0.5
        a3 = (a31**2+a32**2+a33**2)**0.5
        self.qois = {
            '{}.{}'.format(ideal_sim_name,'lammps_Ecoh_min_all'): ecoh,
            '{}.{}'.format(ideal_sim_name, 'lammps_a1_min_all'): a1,
            '{}.{}'.format(ideal_sim_name, 'lammps_a2_min_all'): a2,
            '{}.{}'.format(ideal_sim_name, 'lammps_a3_min_all'): a3
        }
        for i, k in enumerated(StructuralRelaxation.qois_calculated):
            if k not in ['Ecoh', 'a1', 'a2', 'a3']:
                qoi_k = '{}.{}'.format(
                    ideal_sim_name,
                    LammpsStructuralRelaxation.qois_calculated[i]
                )
                qoi_v = simulation_results['{}.{}'.format(ideal_sim_name,k)]
                self.qois[qoi_k] = qoi_v
