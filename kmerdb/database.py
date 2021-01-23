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
import array
from sqlalchemy import create_engine
from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey, Sequence, Text, Boolean
from sqlalchemy.pool import NullPool
#from sqlalchemy.ext.declarative import declarative_base
#from sqlalchemy.orm import sessionmaker

rows_per_loading_transaction = 1000000
#Base = declarative_base()


# ins = kmers.insert()
# conn = engine.connect()
# conn.execute(ins, [{'count': 0}, ...])


def histogram(conn):
    res = conn.execute("SELECT COUNT(count) FROM kmers GROUP BY count ORDER BY count")
    hist = array.array('H')
    for c in res.fetchall():
        hist.append(c[0])
    return hist


class SqliteKdb:
    def __init__(self, filename:str, k:int):
        
        if type(filename) is not str:
            raise TypeError("kmerdb.database.SqliteKdb expects a str as its first positional argument")
        elif not os.path.exists(filename):
            raise OSError("kmerdb.database.SqliteKdb expects an existing .kmerdb file as its first positional argument")
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
        #null_profile = array.array('B')
        null_profile = []
        # FIXME This won't work for certain values of k
        for x in range(self._max_records):
            null_profile.append({
                'count': 0,
                'starts': '[]',
                'reverses': '[]',
                'seqids': '[]'
            })

        with self._engine.connect() as conn:
            while len(null_profile) > 0:
                temp = []
                for y in range(rows_per_loading_transaction):
                    if len(null_profile) > 0:
                        temp.append( null_profile.pop(0) )
                conn.execute(kmers.insert(), temp)

            
    def __exit__(self, exc_type, exc_value, traceback):
        self.conn.close()
        self._engine.dispose()

    def __validate_database(self):
        assert self._get_total_records() == self._max_records, "The number of total records in the database should equal 4^k"

    def _get_total_records(self):
        res = self.conn.execute("SELECT COUNT(id) FROM kmers")
        return res.fetchone()[0]
        
    def _get_nullomers(self):
        res = self.conn.execute("SELECT COUNT(id) FROM kmers WHERE count = 0")
        return res.fetchone()[0]

    def _get_sum_counts(self):
        res = self.conn.execute("SELECT SUM(count) from kmers")
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


    
# class Kmer(Base):
#     __tablename__ = 'kmers'

#     id = Column(Integer, Sequence('kmer_id_seq'), primary_key=True)
#     count = Column(Integer)


#     def __repr__(self):
#         return "<Kmer(id='%s', count='%s')>" % (self.id, self.count)


    
