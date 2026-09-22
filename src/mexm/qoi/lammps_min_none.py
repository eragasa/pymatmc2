from collections import OrderedDict
from mexm.qoi import StaticStructure

class LammpsStaticStructure(StaticStructure):
    qoi_type = 'lmps_min_none'
    qois_calculated = [
            'lammps_Ecoh_min_none',
            'lammps_a1_min_none',
            'lammps_a2_min_none',
            'lammps_a3_min_none',
            'lammps_a11_min_none',
            'lammps_a12_min_none',
            'lammps_a13_min_none',
            'lammps_a21_min_none',
            'lammps_a22_min_none',
            'lammps_a23_min_none',
            'lammps_a31_min_none',
            'lammps_a32_min_none',
            'lammps_a33_min_none',
            'lammps_p11_min_none',
            'lammps_p12_min_none',
            'lammps_p13_min_none',
            'lammps_p21_min_none',
            'lammps_p22_min_none',
            'lammps_p23_min_none',
            'lammps_p31_min_none',
            'lammps_p32_min_none',
            'lammps_p33_min_none'
    ]

    def __init__(self,qoi_name,structures):
        assert isinstance(qoi_name,str)
        assert isinstance(structures, dict)

        StaticStructure.__init__(self,
                                 qoi_name=qoi_name,
                                 structures=structures)

    def determine_simulations(self):
        ideal_structure_name = self.structures['ideal']
        ideal_simulation_type = 'lmps_min_none'
        ideal_simulation_name = "{}.{}".format(
                ideal_structure_name,
                ideal_simulation_type)
        self.add_simulation(
                simulation_id='ideal',
                simulation_name=ideal_simulation_name,
                simulation_type=ideal_simulation_type,
                simulation_structure=ideal_structure_name)

    def calculate_qois(self, results):
        _prefix = '{}.{}'.format(
                self.structures['ideal'],
                'lmps_min_none')

        _e_min_none = results['{}.{}'.format(_prefix,'toten')]
        _n_atoms = results['{}.{}'.format(_prefix,'natoms')]
        _p_tot = results['{}.{}'.format(_prefix,'totpress')]

        # getting the individual components of the H-matrix
        # H = [ [a11,a12,a13]
        #       [a21,a22,a33]
        #       [a31,a32,a33]]

        _a11 = results['{}.{}'.format(_prefix,'a11')]
        _a12 = results['{}.{}'.format(_prefix,'a12')]
        _a13 = results['{}.{}'.format(_prefix,'a13')]
        _a21 = results['{}.{}'.format(_prefix,'a21')]
        _a22 = results['{}.{}'.format(_prefix,'a22')]
        _a23 = results['{}.{}'.format(_prefix,'a23')]
        _a31 = results['{}.{}'.format(_prefix,'a31')]
        _a32 = results['{}.{}'.format(_prefix,'a32')]
        _a33 = results['{}.{}'.format(_prefix,'a33')]


        # getting the individual components of the pressure tensor
        # H = [ [p11,p12,p13]
        #       [p21,p22,p33]
        #       [p31,p32,p33]]
        _p11 = results['{}.{}'.format(_prefix,'p11')]
        _p12 = results['{}.{}'.format(_prefix,'p12')]
        _p13 = results['{}.{}'.format(_prefix,'p13')]
        _p21 = results['{}.{}'.format(_prefix,'p21')]
        _p22 = results['{}.{}'.format(_prefix,'p22')]
        _p23 = results['{}.{}'.format(_prefix,'p23')]
        _p31 = results['{}.{}'.format(_prefix,'p31')]
        _p32 = results['{}.{}'.format(_prefix,'p32')]
        _p33 = results['{}.{}'.format(_prefix,'p33')]

        # calculate the length of cells
        # a = a1
        # b = a2
        # c = a3
        _a1 = (_a11**2+_a12**2+_a13**2)**0.5
        _a2 = (_a21**2+_a22**2+_a23**2)**0.5
        _a3 = (_a31**2+_a32**2+_a33**2)**0.5

        self.qois = OrderedDict()
        self.qois['{}.Ecoh_min_none'.format(_prefix)] = _e_min_none / _n_atoms
        self.qois['{}.a1_min_none'.format(_prefix)] = _a11
        self.qois['{}.a2_min_none'.format(_prefix)] = _a12
        self.qois['{}.a3_min_none'.format(_prefix)] = _a13
        self.qois['{}.a11_min_none'.format(_prefix)] = _a11
        self.qois['{}.a12_min_none'.format(_prefix)] = _a12
        self.qois['{}.a13_min_none'.format(_prefix)] = _a13
        self.qois['{}.a21_min_none'.format(_prefix)] = _a21
        self.qois['{}.a22_min_none'.format(_prefix)] = _a22
        self.qois['{}.a23_min_none'.format(_prefix)] = _a23
        self.qois['{}.a31_min_none'.format(_prefix)] = _a31
        self.qois['{}.a32_min_none'.format(_prefix)] = _a32
        self.qois['{}.a33_min_none'.format(_prefix)] = _a33
        self.qois['{}.p_11_min_none'.format(_prefix)] = _p11
        self.qois['{}.p_12_min_none'.format(_prefix)] = _p12
        self.qois['{}.p_13_min_none'.format(_prefix)] = _p13
        self.qois['{}.p_21_min_none'.format(_prefix)] = _p21
        self.qois['{}.p_22_min_none'.format(_prefix)] = _p22
        self.qois['{}.p_23_min_none'.format(_prefix)] = _p23
        self.qois['{}.p_31_min_none'.format(_prefix)] = _p31
        self.qois['{}.p_32_min_none'.format(_prefix)] = _p32
        self.qois['{}.p_33_min_none'.format(_prefix)] = _p33

    def get_required_variables(self):
        return list(self.variables.keys())
