from pymatmc2 import MultiCell

class MulticellComparer():
    """ Compare two multicells """

    def __init__(self, mc0: MultiCell, mc1: MultiCell):
        if not isinstance(mc0, MultiCell):
            msg = "mc0 argument must be an instance of pymatmc2.MultiCell"
            raise TypeError(msg)

        if not isinstance(mc1, MultiCell):
            msg = "mc1 argument must be an instance of pymatmc2.MultiCell"
            raise TypeError(msg)

        self.mc0 = mc0
        self.mc1 = mc1

    
