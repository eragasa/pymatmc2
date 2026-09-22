
class Oszicar():
    str_scf_step = "N       E                     dE             d eps       ncg     rms          rms(c)"

    """ class to read OSZICARs from VASP

    Properties:
        n_scf (int): total_number of self-consistent steps
        n_ionic (int): total number of ionic relaxation
        electric_scf (list of dict): convergence info for SCF calculations
        ionic_relaxation (list of dict): convergernce info for ionic relaxation steps
    """
    def __init__(self):
        self.path_ = None
        self.n_scf = None
        self.scf_info = None
        self.ionic_info = None

    @property
    def path(self):
        return self.path_

    @path.setter
    def path(self, path):
        if not isinstance(path, str):
            msg = 'path must be a string'
            raise TypeError(msg)

        self.path_ = path

    @property
    def n_ionic(self):
        return len(self.ionic_info)

    def write(self, path):
        raise NotImplementedError

    def read(self, path):
        """ read an OSZICAR file

        Args:
           path (str): path to the OSZICAR file
        """
        
        self.path = path

        with open(self.path, 'r') as f:
            lines = [k.strip() for k in f.readlines()]

        n_lines = len(lines)

        n_scf = 0
        electric_scf = []
        ionic_info = []

        i_line = 0
        while i_line < n_lines:

            if lines[i_line] == Oszicar.str_scf_step:
                electric_scf.append([])
                n_scf += 1
            elif (lines[i_line].startswith('DAV:'))or (lines[i_line].startswith('RMM:')):

                tokens = [k.strip() for k in lines[i_line].split()]
                try:
                    E = float(tokens[2])
                except ValueError:
                    E = None

                try:
                    dE = float(tokens[3])
                except ValueError:
                    dE = None

                try:
                    d_eps = float(tokens[4])
                except ValueError:
                    d_eps = None

                try:
                    ncg = int(tokens[5])
                except ValueError:
                    ncg = None

                try:
                    rms = float(tokens[6])
                except ValueError:
                    rms = None

                try:
                    rms_c = float(tokens[7])
                except ValueError:
                    rms_c = None
                except IndexError:
                    rms_c = None
		
                electric_scf[n_scf-1].append({
                    'E':E,
                    'dE':dE,
                    'd_eps':d_eps,
                    'ncg':ncg,
                    'rms':rms,
                    'rms_c':rms_c
                })
            else:
                tokens = [k.strip() for k in lines[i_line].split()]
                ionic_info.append({
                    'N':int(tokens[0]),
                    'F':float(tokens[2]),
                    'E0':float(tokens[4]),
                    'dE0':float(tokens[7][1:])
                })
            i_line += 1

        # put into attribute, eventually change to properties        
        self.n_scf = n_scf
        self.scf_info = electric_scf
        self.ionic_info = ionic_info
            
