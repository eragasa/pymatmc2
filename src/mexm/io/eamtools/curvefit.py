class EamCurveFitter(object):

    def __init__(self, element_names):

        self.elements = {}
        self.pair_potentials = {}
        for element_name in element_names:
          print(element_name)
          self.elements[element_name] = Element("Ni")

        # initialize varaibles
        self.embedding_energy = {}
        self.pair_potential = {}
        self.electron_density = {}

        print("initializing variables...")

        for el_i in self.elements.keys():
          print('electron_density: rho_{}'.format(el_i))
          print('embedding_energy: F_{}'.format(el_i))
          self.embedding_energy[el_i] = []
          self.electron_density[el_i] = []

        for el_i in self.elements.keys():
          for el_j in self.elements.keys():
            print('pair_potential: p_{}{}'.format(el_i,el_j))
            self.pair_potentials['{}{}'.format(el_i,el_j)] = []
