import os
import re
import shutil

class VaspPotcarError(Exception):
    def __init__(self,*args,**kwargs):
        """ Error class for reading/writing VASP POTCAR IO issues """
        Exception.__init__(self,*args,**kwargs)

class Potcar(object):
    supported_xc_types = [
        'PAW_LDA',
        'PAW_GGA',
        'PAW_PBE'
    ]
    """ object to deal with POTCAR files for VASP simulations

    Args:
        symbols (:obj:`list` of :obj:`str`): a list of chemical symbols
        path (str): the path of the POTCAR
        xc (str): the exchange correlation model.  Default is 'GGA'

    Attributes:
        potcar_dir (str): the path of the VASP potential directory.
        potcars (dict): a dictionary of key-value pairs of symbols to potcar
            type where the key value is the ISO symbol, and the value is the
            the psuedopotnetial type.

    If an :obj:`str` is passed into symbols, then the string will be
    cast into a :obj:`list` of :obj:`str`

    A directory of potentials is not included in pypospack because it is
    protected by copy right.  If the potcar_dir is not set, then pypospack
    will use the environment variables:

    For linux style sytems:
    >>> export VASP_LDA_DIR=$(cd ~/opt/vasp/potentials/LDA;pwd)
    >>> export VASP_GGA_DIR=$(cd ~/opt/vasp/potentials/GGA;pwd)

    for windows systems:
    >>> setx VASP_LDA_DIR "C:\\vasp\\LDA" /M
    >>> setx VASP_GGA_DIR "C:\\vasp\\GGA" /M
    """
    def __init__(self, symbols = None,
                       path= None,
                       xc_type = 'PAW_GGA'):

        self.symbols_ = None
        self.path_ = None
        self._xc_type = None
        self.encut_min_ = None
        self.encut_max_ = None

        self.symbols = symbols
        self.path = path
        self.xc_type = xc_type
        self.models = None

    @property
    def symbols(self):
        return self.symbols_

    @symbols.setter
    def symbols(self, symbols):
        if isinstance(symbols, list):
            if not all([isinstance(s, str) for s in symbols]):
                raise ValueError('all elements of symbols must be strings')
            self.symbols_ = list(symbols)
        elif isinstance(symbols, str):
            self.symbols_ = [symbols]
        elif symbols is None:
            self.symbols_ = None
        else:
            msg = (
                'symbols must be a list of strings consisting of chemical '
                'symbols'
            )
            raise TypeError(msg)

    @property
    def xc_type(self):
        return self._xc_type

    @xc_type.setter
    def xc_type(self, xc_type):
        if isinstance(xc_type, str):
            if xc_type not in Potcar.supported_xc_types:
                msg = (
                    '{xc} is an unsupported exchange correlation functional'
                ).format(xc=xc_type)
                raise VaspPotcarError(msg)
            self._xc_type = xc_type
        else:
            msg = (
                'xc_type argument must be a string'
            )
            raise TypeError(msg)

    @property
    def path(self):
        return self.path_

    @path.setter
    def path(self, path):
        if path is None:
            self.path_ = None
        elif isinstance(path, str):
            # ensure that the path variable is stored as an absolute path
            if os.path.isabs(path):
                self.path_ = path
            else:
                self.path_ = os.path.abspath(path)
        else:
            raise TypeError('path must be a string to the abspath')

    @property
    def encut_max(self):
        return self.encut_max_

    @encut_max.setter
    def encut_max(self, encut_max):
        if isinstance(encut_max, float):
            self.encut_max_ = encut_max
        else:
            self.encut_max_ = float(encut_max)

    @property
    def encut_min(self):
        return self.encut_min_

    @encut_min.setter
    def encut_min(self, encut_min):
        if isinstance(encut_min, float):
            self.encut_min_ = encut_min
        else:
            self.encut_min_ = float(encut_min)

    def read(self, path = "POTCAR"):
        # initialize arrays
        symbols = []
        encut_min = []
        encut_max = []
        models = []
        xc = []

        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if 'TITEL' in line:
                    
                    symbol = line.split('=')[1].strip().split(' ')[1]
                    symbols.append(symbol)

                    xc.append(
                        line.split('=')[1].strip().split(' ')[0]
                    )
                elif 'ENMIN' in line:
                    obj_re = re.match('ENMAX  =  (.*); ENMIN  =  (.*) eV.*',line,re.M|re.I)
                    enmax = float(obj_re.group(1))
                    enmin = float(obj_re.group(2))
                    encut_min.append(enmin)
                    encut_max.append(enmax)
                elif "VRHFIN" in line:
                    models.append(line.split('=')[1].strip())

        self.symbols = symbols
        self.encut_min = min(encut_min)
        self.encut_max = max(encut_max)

        # check to see if all the exchange correlation functionals are the same        
        assert xc != []
        assert xc.count(xc[0]) == len(xc)

        if xc[0] in self.supported_xc_types:
            self.xc_type = xc[0]
        else:
            msg = 'unknown xc_type: {}'.format(xc[0])
            raise ValueError(msg)
        
    def write(self, 
              path = 'POTCAR', 
              src = None, 
              potcars = None):
        """ write POTCAR file

        Arguments:
            path (str): If none, then the method will use the internal
                attribute.  Other

        Raises:
            pypospack.io.vasp.VaspPotcarError
        """

        self.path = path

        #  if a source potcar is given just copy the file
        if src is not None:
            shutil.copy(src=src, dst=self.path)

        self.__set_xc_directory()

        #  get full path of potcars
        if potcars is not None:
            self.potcars = potcars.copy()
        else:
            self.__find_potcar_files()


        with open(self.path,'w') as f_out:
            for s in self.symbols:
                with open(self.potcars[s],'r') as f_in:
                    # avoid reading large files into memory
                    chunk_sz = 1024*1024*10 # 10MB
                    shutil.copyfileobj(f_in,f_out, chunk_sz)

    def __write_potcar_by_copy(self,src_fn,dst_fn):
        with open(dst_fn) as f_out:
            with open(src_fn) as f_in:
                chunk_sz = 1024*1024*10
                shutil.copyfileobj(f_in,f_out,chunk_sz)

    def __find_potcar_files(self):
        # a more intelligent way of selecting potcars based upon VASP
        # recommendations should be implemented
        # http://cms.mpi.univie.ac.at/vasp/vasp/Recommended_PAW_potentials_DFT_calculations_using_vasp_5_2.html
        self.potcars = {} # initialize
        for s in self.symbols:
            
            potcar_path = os.path.join(self.potcar_dir,s,'POTCAR')
            potcar_new_path = os.path.join(
                self.potcar_dir,"{}_new".format(s),'POTCAR'
            )
            potcar_hard_path = os.path.join(
                self.potcar_dir,"{}_h".format(s),'POTCAR'
            )
            potcar_sv_path = os.path.join(
                self.potcar_dir, '{}_sv'.format(s), 'POTCAR'
            )
            if os.path.isfile(potcar_path):
                # try VASP_XC_DIR/symbol/POTCAR
                self.potcars[s] = potcar_path

            elif os.path.isfile(potcar_new_path):
                # try VASP_XC_DIR/symbol_new/POTCAR
                self.potcars[s] = os.path.join(
                    self.potcar_dir,"{}_new".format(s),'POTCAR'
                )

            elif os.path.isfile(potcar_hard_path):
                # try VASP_XC_DIR/symbols_h/POTCAR
                self.potcars[s] = potcar_hard_path
            elif os.path.isfile(potcar_sv_path):
                self.potcars[s] = potcar_sv_path
            else:
                print(os.path.join(self.potcar_dir,"{}_new".format(s),'POTCAR'))
                msg = 'cannot find a POTCAR file for {}.{}\n'.format(self.xc_type,s)
                msg += 'potcar_dir:{}'.format(self.potcar_dir)
                raise VaspPotcarError(msg)

    @staticmethod
    def get_xc_directory(self, xc_type):
        if xc_type == 'GGA':
            try:
                return os.environ['VASP_GGA_DIR']
            except KeyError:
                msg = (
                    'need to set the environment variable VASP_GGA_DIR'
                )
                raise VaspPotcarError(msg)
        elif xc_type == 'LDA':
            try:
                return os.environ['VASP_LDA_DIR']
            except KeyError:
                msg = (
                    'need to set the environment variable VASP_LDA_DIR'
                )
                raise VaspPotcarError(msg)
        else:
            msg = (
                'unknown xc_type:{}'
            ).format(xc_type)
            raise VaspPotcarError(msg)

    def __set_xc_directory(self) -> str:
        """ sets the directory of the exchange correlation functional

        the exchange correlation functional is set to either 'LDA' or 'GGA'
        based upon environment variables

        Raises:
            (pypospack.io.vasp.VaspPotcarError): if the directory is not set
                 as an environment variable.
        Returns:
            str: the path to the XC directory.
        """
        if self.xc_type in ['PAW_GGA', 'PAW_PBE']:
            try:
                self.potcar_dir = os.environ['VASP_POTPAW_GGA']
            except KeyError:
                msg = ('need to set environment variable VASP_POTPAW_GGA to the '
                       'location of the VASP GGA-PBE potential files' )
                raise VaspPotcarError(msg)
        elif self.xc_type == 'PAW_LDA':
            try:
                self.potcar_dir = os.environ['VASP_POTPAW_LDA']
            except KeyError:
                msg = ('need to set environment variable VASP_POTPAW_LDA to the '
                       'location of the VASP LDA-CA potential files' )
        else:
            raise ValueError('unknown xc_type:{}'.format(self.xc_type))
        return self.potcar_dir
