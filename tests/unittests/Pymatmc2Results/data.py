from pymatmc2 import Pymatmc2Configuration

class Pymatmc2Data:
    def __init__(self):
        self._configuration = None
        self.path = 'pymatmc2.data'
    @property
    def configuration(self) -> Pymatmc2Configuration:
        return self._configuration

    @configuration.setter
    def configuration(self, configuration: Pymatmc2Configuration):
        self._configuration = configuration


    def write_header(self):
        
    def write_results(self, i_iteration: int, mc_type, mc_initial: MultiCell, mc_candidate: MultiCell, mc_final_multicell):
        mcs = {}
        mcs['initial'] = mc_initial
        mcs['candidate'] = mc_candidate
        mcs['final'] = mc_final

        symbols = list(self.configuration.total_concentration.values())

        'iteration'
        'mutate_type'
        'initial.fcc1.Au'
        'initial.fcc1.Pt'
	'initial.fcc1.f'
        'initial.fcc1.toten'
        'initial.fcc1.eatom'
        'initial.fcc2.Au'
        'initial.fcc2.Pt'
        'initial.fcc2.f'
        'initial.fcc2.toten'
        'initial.fcc2.eatom'
        'initial.total.toten'
        'initial.total.eatom'
        
