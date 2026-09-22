from mexm.qoi import StructuralRelaxation

class VaspStructuralRelaxation(StructuralRelaxation):
    qoi_type = 'vasp_min_all'
    is_base_class = False
    qois_calculated = [
        'vasp_Ecoh_min_all',
        'vasp_a1_min_all',
        'vasp_a2_min_all',
        'vasp_a3_min_all',
        'vasp_a11_min_all',
        'vasp_a12_min_all',
        'vasp_a13_min_all',
        'vasp_a21_min_all',
        'vasp_a22_min_all',
        'vasp_a23_min_all',
        'vasp_a31_min_all',
        'vasp_a32_min_all',
        'vasp_a33_min_all',
        'vasp_p11_min_all',
        'vasp_p12_min_all',
        'vasp_p13_min_all',
        'vasp_p21_min_all',
        'vasp_p22_min_all',
        'vasp_p23_min_all',
        'vasp_p31_min_all',
        'vasp_p32_min_all',
        'vasp_p33_min_all'
    ]

    def __init__(self, qoi_name, structures):
        assert isinstance(qoi_name,str)
        assert isinstance(structures, dict)

        StructuralRelaxation.__init__(self, qoi_name=qoi_name, structures=structures)
