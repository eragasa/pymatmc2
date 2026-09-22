import os
from abc import ABC
from abc import abstractmethod

class DatabaseAdapter(ABC):
    def __init__(self, dbpath: str = None):
        self.db_path = os.path.abspath(dbpath)


from mexm.database.sqlite3 import Sqlite3Adapter

