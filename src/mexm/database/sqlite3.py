import os
import sqlite3
from mexm.database import DatabaseAdapter

__all__ = ['Sqlite3Adapter']
sqlstrings = {
    "create_job_transaction_table": (
        'CREATE TABLE'
    )
}

class Sqlite3Adapter(DatabaseAdapter):
    sqlstrings = sqlstrings
    def __init__(
        self,
        hostname=None,
        port=None,
        dbname="mexm.db",
        user=None,
        password=None
    ):
        """

        Arguments:
            hostname (str): the hostname to the database
            port (int): the port to the database    
            dbname (str): the name of the database or filename
            user (str): the username to the database
            password (str): the password to the database
        """
        # the following arguments are not used
        # but they must be defined for the abstract class
        self.db_conn = sqlite3.connect(
            os.path.abspath(dbname)
        )