from collections import OrderedDict
from mexm.qoi import Qoi

class StaticStructure(Qoi):
    qoi_type = 'min_none'
    is_base_class = True
    qois_calculated = [
        'Ecoh_min_none',
        'a1_min_none','a2_min_none','a3_min_none',
        'a11_min_none','a12_min_none','a13_min_none',
        'a21_min_none','a22_min_none','a23_min_none',
        'a31_min_none','a32_min_none','a33_min_none',
        'p11_min_none','p12_min_none','p13_min_none',
        'p21_min_none','p22_min_none','p23_min_none',
        'p31_min_none','p32_min_none','p33_min_none'
    ]

    def __init__(self, qoi_name, structures):
        assert isinstance(qoi_name,str)
        assert isinstance(structures, dict)

        Qoi.__init__(self, qoi_name=qoi_name, structures=structures)
