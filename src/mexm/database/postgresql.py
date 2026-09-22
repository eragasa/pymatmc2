import psycopg2
from mexm.database import DatabaseAdapter
class PsqlAdapter(DatabaseAdapter):
    
    def __init__(
        self,
        hostname='localhost',
        port=5432,
        dbname="mexm",
        user="username",
        password="password"
    ):
        self.db_conn=psycopg2.connect(
            self.get_sql_connection_string(
                hostname=hostname,
                port=port,
                dbname=dbname,
                user=user,
                password=password
            )
        )

    def get_sql_connection_string(
        self, 
        hostname,
        dbname,
        port, 
        user, 
        password
    ):
        str_connection = (
            "host={host} "
            "port={port} "
            "dbname={dbname} "
            "user={user} "
            "password={password}"
        ).format(
            host=hostname, 
            port=port, 
            dbname=dbname,
            user=user,
            password=password
        )
        return str_connection
