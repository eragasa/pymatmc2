class VaspIncarError(Exception):
    def __init__(self,*args,**kwargs):
        """Error class for reading/writing VASP INCAR IO issues """
        Exception.__init__(self,*args,**kwargs)