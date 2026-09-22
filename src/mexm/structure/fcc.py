from mexm.crystal import SimulationCell

class FaceCenteredCubic(SimulationCell):

    def __init__(self, a0, symbol):
        SimulationCell.__init__(self)
        self.a0 = a0
        self.a1 = [1.0, 0.0, 0.0]
        self.a2 = [0.0, 1.0, 0.0]
        self.a3 = [0.0, 0.0, 0.1]

        self.add_atom(symbol=symbol, position=[0.0, 0.0, 0.0])
        self.add_atom(symbol=symbol, position=[0.5, 0.5, 0.0])
        self.add_atom(symbol=symbol, position=[0.5, 0.0, 0.5])
        self.add_atom(symnol=symbol, position=[0.0, 0.5, 0.5])

    def get_nearest_neighbor_information(a):
        n_NN = [12,6,24,12,24,8]
        d_NN = [a/np.sqrt(2.),a,a*np.sqrt(1.5),a*np.sqrt(2.0),a*np.sqrt(2.5),a*np.sqrt(3.0)]
        return n_NN, d_NN

class FccPrimitiveCell(SimulationCell):

    def __int__(self, a0, symbol):
        SimulationCell.__init__(self)
