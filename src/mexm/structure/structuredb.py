
class StructureDatabase(object):
    """ structure database

    Attributes:
        filename(str): file to read/write the yaml file
        directory(str): the directory of the structure database
        structures(dict): key is structure name, value is another dict with
            key/value pairs for 'filename' and 'filetype'
    """

    def __init__(self):
        self.filename = 'pypospack.structure.yaml'
        self.directory = None
        self.structures = {}

    def add_structure(self,name,filename,filetype):
        self.structures[name] = {'filename':filename,'filetype':filetype}

    def contains(self,structure):
        """ check to see that the structure in the structure database

        Args:
            structure(str):

        Returns:
            bool: True if structure is in structure database. False if the
                structure is not in the structure database
        """

        return structure in self.structures.keys()

    def read(self,fname = None):
        """ read qoi configuration from yaml file

        Args:
            fname(str): file to read yaml file from.  If no argument is passed
                then use the filename attribute.  If the filename is set, then
                the filename attribute is also set.
        """

        # set the attribute if not none
        if fname is not None:
            self.filename = fname

        try:
            self.structure_db = yaml.load(open(self.filename))
        except:
            raise

        self.directory = self.structure_db['directory']
        self.structures = copy.deepcopy(self.structure_db['structures'])

    def write(self,fname = None):
        if fname is not None:
            self.filename = fname

        # marshall attributes into a dictionary
        self.structure_db = {}
        self.structure_db['directory'] = self.directory
        self.structure_db['structures'] = {}
        self.structure_db['structures'] = copy.deepcopy(self.structures)

        # dump to as yaml file
        with open(fname,'w') as f:
            yaml.dump(self.structure_db, f, default_flow_style=False)

    def check(self):
        """sanity checks for the fitting database

        This method checks the fitting database for the following errors: (1)
        the structure database exists, (2) files for the structure database
        exist

        Returns:
            str: returns a string if there is a problem
        Raises:
            ValueError: if there is a problem with the configuration
        """

        src_dir = self.directory

        # check to see if the source directory exists
        if not os.path.exists(src_dir):
            err_msg = "cannot find simulation directory\n"
            err_msg += "\tcurrent_working_directory:{}\n".format(os.getcwd())
            err_msg += "\tstructure_db_directory:{}\n".format(src_dir)
            return err_msg

        # check to see if the source directory is a directory
        if not os.path.isdir(src_dir):
            err_msg = "path exists, is not a directory\n"
            err_msg += "\tcurrent_working_directory:{}".format(os.getcwd())
            err_msg += "\tstructure_db_directory:{}\n".format(src_dir)
            return err_msg

        # check to see if files exist in the source directory
        files_exist = True
        msg = "structure files are missing:\n"
        for name, v in self.structures.items():
            filename = os.path.join(src_dir,v['filename'])
            if not os.path.isfile(filename):
                files_exist = False
                msg += "\t{}:{}\n".format(name,filename)

        if not files_exist:
            return msg
        else:
            return True

    def get_structure_dict(self,name):
        """

        Args:
            name(str): name of the structure
        """

        structure_db_dir = self.directory
        structure_filename = self.structures[name]['filename']
        structure_dict = {}
        structure_dict['name'] = name
        structure_dict['filename'] = os.path.join(
                structure_db_dir,
                structure_filename)

        return copy.deepcopy(structure_dict)
