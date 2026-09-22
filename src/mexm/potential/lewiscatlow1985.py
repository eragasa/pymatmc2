from collections import OrderedDict

# Lewis and Catlow, J Phys Solid State Phys 18 (1985) 1149-1165
# Catlow, Proc. B. Lond A. 353, 553-561 (1997)

divalent_cations = ['Mg','Ca','Sr','Ba','Mn','Fe','Co','Ni','Zn']
trivalent_cations = ['Al','Ti','V','Cr','V','Cr','Fe']
tetravalent_cations = ['Ge','Sn','Ti','Ce','U','Th']
all_cations = divalent_cations + trivalent_cations + tetravalent_cations
all_other_parameters = {
    'MgO_A': 821.6, 'MgO_rho':0.3242, 'MgO_C': 0.0,
    'CaO_A':1228.9, 'CaO_rho':0.3372, 'CaO_C': 0.0,
    'SrO_A':1400.0, 'SrO_rho':0.3500, 'SrO_C': 0.0,
    'BaO_A': 931.7, 'BaO_rho':0.3949, 'BaO_C': 0.0,
    'MnO_A': 715.8, 'MnO_rho':0.3464, 'MnO_C': 0.0,
    'FeO_A': 694.1, 'FeO_rho':0.3399, 'FeO_C': 0.0,
    'CoO_A': 696.3, 'CoO_rho':0.3362, 'CoO_C': 0.0,
    'NiO_A': 683.5, 'NiO_rho':0.3332, 'NiO_C': 0.0,
    'ZnO_A': 499.6, 'ZnO_rho':0.3595, 'ZnO_C': 0.0,
    'AlO_A':1114.9, 'AlO_rho':0.3118, 'AlO_C': 0.0,
    'TiO_A':1715.7, 'TiO_rho':0.3069, 'TiO_C': 0.0,
    'VO_A': 1790.2, 'VO_rho': 0.3061, 'VO_C':  0.0,
    'CrO_A':1734.1, 'CrO_rho':0.3010, 'CrO_C': 0.0,
    'FeO_A':1102.4, 'FeO_rho':0.3299, 'FeO_C': 0.0,
    'GeO_A':1035.5, 'GeO_rho':0.3464, 'GeO_C': 0.0,
    'SnO_A': 938.7, 'SnO_rho':0.3813, 'SnO_C': 0.0,
    'TiO_A': 754.2, 'TiO_rho':0.3879, 'TiO_C': 0.0,
    'CeO_A':1013.6, 'CeO_rho':0.3949, 'CeO_C': 0.0,
    'OO_A':22764.3, 'OO_rho': 0.1490, 'OO_C':27.88
}

parameters = {}
for s in divalent_cations:
    parameters['{}_chrg'.format(s)] = 2.0
for s in trivalent_cations:
    parameters['{}_chrg'.format(s)] = 3.0
for s in tetravalent_cations:
    parameters['{}_chrg'.format(s)] = 4.0
parameters['O_chrg'] = -2.0
for i1, s1 in enumerate(all_cations):
    for i2, s2 in enumerate(all_cations):
        if i1 <= i2:
            parameters['{}{}_A'.format(s1,s2)] = 0.0
            parameters['{}{}_rho'.format(s1,s2)] = 0.5
            parameters['{}{}_C'.format(s1,s2)] = 0.0
for parameter_name, parameter_value in all_other_parameters.items():
    parameters[parameter_name] = parameter_value

LewisCatlow1985 = OrderedDict([
    ('potential_type','buckingham'),
    ('symbols',all_cations),
    ('parameters',parameters),
    ('references',['lewiscatlow1985_buck','catlow1977'])
])
