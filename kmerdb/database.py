'''
   Copyright 2020 Matthew Ralston

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

'''




import logging
logger = logging.getLogger(__file__)


import io
import os
import json
#import sqlite3
import sys
import hashlib
import random

from datetime import datetime
from sqlalchemy import create_engine
from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey, Sequence, JSON, Boolean

from sqlalchemy.pool import NullPool
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker


import psycopg2
from psycopg2.extensions import parse_dsn
from psycopg2 import sql

rows_per_loading_transaction = 65535
Base = declarative_base()


# ins = kmers.insert()
# conn = engine.connect()
# conn.execute(ins, [{'count': 0}, ...])





def histogram(conn):
    res = conn.execute("SELECT COUNT(count) FROM kmers GROUP BY count ORDER BY count")
    hist = array.array('H')
    for c in res.fetchall():
        hist.append(c[0])
    return hist

class PostgresKdb:
    def __init__(self, k:int, connection_string: str, tablename=None, filename=None):
        
        if type(k) is not int:
            raise TypeError("kmerdb.database.PostgresKdb expects an int as its first positional argument")
        elif type(connection_string) is not str:
            raise TypeError("kmerdb.database.PostgresKdb expects a str as its second positional argument")
        elif filename is not None and type(filename) is not str:
            raise TypeError("kmerdb.database.PostgresKdb expects the keyword argument 'filename' to be a str")
        elif filename is not None and not os.path.exists(filename):
            raise OSError("kmerdb.database.PostgresKdb expects the keyword argument 'filename' to be an existing filepath")
        elif tablename is not None and type(tablename) is not str:
            raise TypeError("kmerdb.database.PostgresKdb expects the keyword argument 'tablename' to be a str")
        self._max_records = 4**k


        if filename is not None and tablename is None:
            m = hashlib.md5()
            #m.update((filename + str(datetime.utcnow()) + str(os.getpid())).encode('utf-8'))

            # Database tablename hash format: filename, random, pid, time
            hash_combination = ("{}_{}_{}_{}".format(filename, random.SystemRandom(), str(os.getpid), str(datetime.utcnow()))).encode('utf-8') # UTF8 encoded text for the hashlib
            m.update(hash_combination)
            
            
        if tablename is not None:
            self._tablename = tablename
        else:
            self._tablename = "kdb_{0}".format(m.hexdigest()[0:15])
        try:
            self._engine = create_engine(connection_string, echo=False, poolclass=NullPool)
        except Exception as e:
            logger.error("Error while establishing connection to the PostgreSQL database...")
            logger.error(e)
            sys.exit(1)


        self.conn = self._engine.connect()
        #self.conn = psycopg2.connect(**parse_dsn(connection_string))

        if tablename is None:
            Session = sessionmaker(bind=self._engine)
            session = Session()
            self.__make_empty_database(session)

        self._id = 1

    def __exit__(self, exc_type, exc_value, traceback):
        if self.conn is not None:
            self.conn.close()
        self._engine.dispose()
        

    def __iter__(self):
        return self

    def __next__(self):
        with self.conn.begin():
            res = self.conn.execute("SELECT * FROM {} WHERE id = %s".format(self._tablename), self._id)
        self._id += 1
        results = res.fetchone()
        if results is None or len(results) == 0:
            raise StopIteration
        else:
            return results


    def _get_total_records(self):
        with self.conn.begin():
            res = self.conn.execute("SELECT COUNT(id) FROM {}".format(self._tablename))
        return res.fetchone()[0]
        
    def _get_nullomers(self):
        with self.conn.begin():
            res = self.conn.execute("SELECT COUNT(id) FROM {} WHERE count = 0".format(self._tablename))
        return res.fetchone()[0]

    def _get_sum_counts(self):
        with self.conn.begin():
            res = self.conn.execute("SELECT SUM(count) from {}".format(self._tablename))
        return res.fetchone()[0]

        
    def __make_empty_database(self, session):
        logger.info("Initializing database...")
        class Kmer(Base):
            __tablename__ = self._tablename

            id = Column(Integer, primary_key=True)
            count = Column(Integer)
            starts = Column(JSON)
            reverses = Column(JSON)
            seqids = Column(JSON)

        # reads = Table('reads', metadata,
        #               Column('read_id', String),
        #               Column('kmer_id', None, ForeignKey('kmers.id'))
        # )
        # Create tables
        Base.metadata.create_all(self._engine)
        logger.info("Loading {0} records at a time into the PostgreSQL database, in case there are limits to memory that are imposed.".format(rows_per_loading_transaction))
        try:
            logger.debug("Populating...")
            for x in range(self._max_records):
                z = Kmer(count=0, starts=[], reverses=[], seqids=[])
                session.add(z)
                if x % rows_per_loading_transaction == 0:
                    session.commit()
        
        except Exception as e:
            session.rollback()
            raise e
        finally:
            session.commit()
            session.close()
            with self.conn.begin():
                self.conn.execute("CREATE INDEX {0}_idx ON {0}(id)".format(self._tablename))
        

class SqliteKdb:
    def __init__(self, filename:str, k:int):
        
        if type(filename) is not str:
            raise TypeError("kmerdb.database.SqliteKdb expects a str as its first positional argument")
        elif not os.path.exists(filename):
            raise OSError("kmerdb.database.SqliteKdb expects an existing Sqlite3 file as its first positional argument")
        elif type(k) is not int:
            raise TypeError("kmerdb.database.SqliteKdb expects an int as its second positional argument")
        self._max_records = 4**k
        self.filepath = filename
        self._engine = create_engine("sqlite:///{0}".format(filename), echo=False, poolclass=NullPool)
        if not os.path.exists(filename) or os.path.getsize(filename) == 0:
            logger.info("No data exists in '{0}', initializing empty database".format(filename))
            ##SQLAlchemy 1.3
            self.__make_empty_database()
        self.conn = self._engine.connect()
        self.__validate_database()
        self._id = 1

        
    def __make_empty_database(self):
        metadata = MetaData()
        kmers = Table('kmers', metadata,
                      Column('id', Integer, Sequence('kmer_id_seq'), primary_key=True),
                      Column('count', Integer),
                      Column('starts', Text),
                      Column('reverses', Text),
                      Column('seqids', Text)
        )
        # reads = Table('reads', metadata,
        #               Column('read_id', String),
        #               Column('kmer_id', None, ForeignKey('kmers.id'))
        # )
        # Create tables
        metadata.create_all(self._engine)
        logger.debug("Populating a length 4^k = {0} profile in Python space for loading into the SQLite3 database".format(self._max_records))
        logger.info("Loading {0} records at a time into the SQLite3 database, in case there are limits to memory that are imposed.".format(rows_per_loading_transaction))
        with self._engine.connect() as conn:
            logger.debug("Populating...")
            for x in range(self._max_records):
                conn.execute(kmers.insert(), {'count': 0, 'starts': '[]', 'reverses': '[]', 'seqids': '[]'})
                # if x % rows_per_loading_transaction == 0:
                #     conn.commit()
                # sys.stderr.write(".")
            logger.info("Deconstructing the sqlite connection")
        logger.info("Finished populating the empty database")
        
        
    def __exit__(self, exc_type, exc_value, traceback):
        self.conn.close()
        self._engine.dispose()

    def __validate_database(self):
        assert self._get_total_records() == self._max_records, "The number of total records in the database should equal 4^k"

    def _get_total_records(self):
        res = self.conn.execute("SELECT COUNT(id) FROM {0}".format(self._tablename))
        return res.fetchone()[0]
        
    def _get_nullomers(self):
        res = self.conn.execute("SELECT COUNT(id) FROM {0} WHERE count = 0".format(self._tablename))
        return res.fetchone()[0]

    def _get_sum_counts(self):
        res = self.conn.execute("SELECT SUM(count) from {0}".format(self._tablename))
        return res.fetchone()[0]

    def _load_metadata(self):
        """ Not implemented """
        return

    
    def __iter__(self):
        return self

    def __next__(self):
        res = self.conn.execute("SELECT * FROM kmers WHERE id = ?", self._id)
        self._id += 1
        results = res.fetchone()
        if results is None or len(results) == 0:
            raise StopIteration
        else:
            return results


