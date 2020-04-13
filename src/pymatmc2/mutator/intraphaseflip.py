from pymatmc2.mutator import MultiCellMutateAlgorithm
from pymatmc2 import MultiCell

class IntraphaseFlip(MultiCellMutateAlgorithm):
    mutate_type = 'intraphase_flip'

    def mutate_multicell(self, multicell: MultiCell) -> MultiCell:
        """ create a new candidate multicell

        Create a deepcopy of multicell to the attribute multicell_initial.


        Arguments:
            multicell (MultiCell): the initial multicell from which to get a    candidate
        """
        self.multicell_initial = deepcopy(multicell)
        rtn_multicell = MultiCell.initialize_from_obj(multicell=multicell)
        while True:
            try:
                for phase in multicell.simulations:
                    simulation = multicell.simulations[phase]

                    # extract the structure file
                    if isinstance(simulation, VaspSimulation):
                        if isinstance(simulation.contcar, Poscar):
                            cell = deepcopy(simulation.contcar)
                        else:
                            cell = deepcopy(simulation.poscar)
                    else:
                        raise ValueError('unknown simulation type')
                    assert isinstance(cell, SimulationCell)

                    # mutate the cell
                    # select atom at random
                    idx_atoms = list(range(cell.n_atoms))
                    idx_atom = np.random.choice(idx_atoms, 1)[0]
                    old_symbol = cell.atomic_basis[idx_atom].symbol
                    
                    # flip the symbol
                    symbols = cell.symbols
                    symbols.remove(old_symbol)
                    new_symbol = np.random.choice(symbols, 1)[0]
                    
                    # assign new symbol
                    cell.atomic_basis[idx_atom].symbol = new_symbol
                
                    if isinstance(simulation, VaspSimulation):
                        rtn_multicell.simulations[phase].poscar \
                            = Poscar.initialize_from_object(obj=cell)
                    else:
                        raise ValueError("unknown simulation type")
                    
                # this will raise a LinAlg error if rank deficient
                rtn_multicell.phase_molar_fraction
                break
            except linalg.LinAlgError:
                pass
        
        self.multicell_candidate = rtn_multicell
        return rtn_multicell

    def acceptance_probability(
        self, 
        E0: float, 
        E1: float, 
        temperature: float
    ):
        kB = constants.BOLTZMANN
        T = temperature
	beta = 1/kB/T
        n_phases = len(configuration.cell_names)
 
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')


            p_accept = min(1,np.exp(-n_phases*(E1-E0)*beta))

        self.p_accept = {}
        self.p_accept['total'] = p_accept
        return p_accept

    def accept_or_reject(
        self,
        multicell_initial: MultiCell, 
        multicell_candidate: MultiCell,
        temperature: float
    ) -> Tuple[bool, MultiCell]:
        

        E0 = multicell_initial.total_energy
        E1 = multicell_candidate.total_energy

        if E1 < E0:
            is_accept = True
        else:
            p_accept = self.acceptance_probability(
                E0 = E0, 
                E1 = E1, 
                temperature = temperature
            )
            
            self.p = np.random.random()
            if self.p < p_accept:
                is_accept = True
            else:
                is_accept = False

        if is_accept:
            return True, multicell_candidate
        else:
            return False, multicell_initial
