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


import os
import sys
import yaml
import json
import time
from datetime import datetime
from math import ceil
from itertools import chain, repeat
from concurrent.futures import ThreadPoolExecutor

from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.attributes import flag_modified
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey, Sequence, JSON, Boolean

from psycopg2 import sql
import tempfile

#import threading
from multiprocessing import Pool


from kmerdb import seqparser, database, kmer

import logging
logger = logging.getLogger(__file__)




Base = declarative_base()




def parsefile(filepath:str, k:int, connection_string:str, p:int=1, rows_per_batch:int=100000, b:int=50000, n:int=1000, stranded:bool=True, all_metadata:bool=False):
    """Parse a single sequence file in blocks/chunks with multiprocessing support

    :param filepath: Path to a fasta or fastq file
    :type filepath: str
    :param k: Choice of k to shred k-mers with
    :type k: int
    :param p: Number of processes
    :type p: int
    :param b: Number of reads (per block) to process in parallel
    :type b: int
    :param stranded: Strand specificity argument for k-mer shredding process
    :type stranded: bool
    :returns: (db, header_dictionary) header_dictionary is the file's metadata for the header block
    :rtype: (kdb.database.SqliteKdb, dict)

    """
    if type(filepath) is not str:
        raise TypeError("kmerdb.parse.parsefile expects a str as its first positional argument")
    elif not os.path.exists(filepath):
        raise OSError("kmerdb.parse.parsefile could not find the file '{0}' on the filesystem".format(filepath))
    elif type(k) is not int:
        raise TypeError("kmerdb.parse.parsefile expects an int as its second positional argument")
    elif type(connection_string) is not str:
        raise TypeError("kmerdb.parse.parsefile expects a str as its third positional argument")
    elif type(p) is not int:
        raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'p' to be an int")
    elif type(b) is not int:
        raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'b' to be an int")
    elif type(n) is not int:
        raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'n' to be an int")
    elif type(stranded) is not bool:
        raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'stranded' to be a bool")
    # Create temporary SQLite3 database file for on-disk k-mer counting
    # temp = tempfile.NamedTemporaryFile(mode="w+", suffix=".sqlite3", delete=False)
    # logger.debug("Creating temporary database to tally k-mers: '{0}'".format(temp.name))
    # temp.close()

    logger.debug("Creating empty database table to store the k-mers counts.")
    db = database.PostgresKdb(k, connection_string, filename=filepath)
    from kmerdb import kmer


    
    data = {} # This is the dictionary of tuples, keyed on k-mer id, and containing 3-tuples ('kmer_id', 'read/start/reverse')
    keys = set()
    rows = []
    nullomers = set()
    try:
        # Build fasta/fastq parser object to stream reads into memory
        logger.debug("Initializing parser...")
        seqprsr = seqparser.SeqParser(filepath, b, k)
        fasta = not seqprsr.fastq # Look inside the seqprsr object for the type of file
        logger.info("Constructing multiprocessing pool with {0} processors".format(p))
        if not fasta: # Cannot process fast files in parallel due to spawning issue.
            pool = Pool(processes=p) # A multiprocessing pool of depth 'p'


        # Instantiate the kmer class
        Kmer = kmer.Kmers(k, strand_specific=stranded, fasta=fasta, all_metadata=all_metadata) # A wrapper class to shred k-mers with

        recs = [r for r in seqprsr] # A block of exactly 'b' reads-per-block to process in parallel
        if not fasta:
            logger.debug("Read exactly b=={0}=={1} records from the {2} seqparser object".format(b, len(recs), s))
            assert len(recs) == b, "The seqparser should return exactly {0} records at a time".format(b)
        else:
            logger.debug("Skipping the block size assertion for fasta files")
        logger.info("Read {0} sequences from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"))

        while len(recs): # While the seqprsr continues to produce blocks of reads
            # Run each read through the shred method
            num_recs = len(recs)

            if fasta:
                list_of_dicts = list(map(Kmer.shred, recs))
            else:
                list_of_dicts = list(map(Kmer.shred, recs))


            logger.info("Shredded up {0} sequences over {1} parallel cores, like a cheesesteak".format(len(list_of_dicts), p))
            #logger.debug("Everything in list_of_dicts is perfect, everything past here is garbage")
            
            # Flatmap to 'kmers', the dictionary of {'id': read_id, 'kmers': [ ... ]}
            kmer_ids = list(chain.from_iterable(map(lambda x: x['kmers'], list_of_dicts)))
            #logger.debug(kmer_ids)
            logger.debug("Initial k-mers as a dictionary from each of {0} sequences".format(num_recs))
            sus = set()
            logger.info("Removing any eroneous k-mers without an id, cause by 'N' content, creating new gaps in the k-mer profile...")
            logger.info("Will substitute IUPAC residues with their respective residues, adding one count for each")
            num_sus = 0
            for i, k in enumerate(kmer_ids):
                if k is None:
                    sus.add(i) # This propagates the N removal, by adding those to the set for removal later.
                    num_sus += 1

            logger.info("Created {0} gaps in the k-mer profile to clean it of N-content".format(num_sus))
            if num_sus > 0:
                logger.warning("{0} uncountable k-mer identified".format(num_sus))
                logger.warning("Likely caused by the presence of the unspecified nucleotide 'N' in a .fasta file")
                logger.warning("Will remove from the final committed profile")

            ########################
            # K-mer cleaning
            #
            # 
            # The following procedure, removes k-mers that don't have a formal k-mer id
            # More specifically, it eliminates entries that were formally loaded into the k-mers array,
            # and thus would be a minor percentage of the total number of k-mers,
            # Specifically, k-mers or subsequences with N can be excluded from the profile
            # To essentially eliminate an area of sequence and create an artificial gap in the "retained" graph.
            # The graph can be explicitly collapsed when you want to, by supplying reads or sequence information, that is contiguous "in your judgement".
            # By introducing these artificial gaps, we expose the errors or unspecified bases through data.
            #
            logger.info("Found {0} invalid or enumerable k-mers, sequences which have no defined index, i".format(num_sus))
            logger.debug("In other words, one or more characters in your .fasta or .fastq file were 'N' or otherwise contained IUPAC characters that need to be substituted.")
            logger.debug("We have chosen to omit k-mer with unspecified nucleotides 'N', and these make gaps in our database explicitly.")
            logger.info("These can be recovered from the reads if they are needed, and the read ids may be found with the experimental --all-metadata flag")
            if all_metadata:

                logger.warning("\n\n\nGENERATING ALL METADATA\n\n\nThis is extremely expensive and experimental. You have been warned.\n\n\n")
                reads = list(chain.from_iterable(map(lambda x: x['seqids'], list_of_dicts)))
                starts = list(chain.from_iterable(map(lambda x: x['starts'], list_of_dicts)))
                reverses = list(chain.from_iterable(map(lambda x: x['reverses'], list_of_dicts)))
                logger.debug("Checking each k-mer for N-content and IUPAC substitution")
                for i, x in enumerate(kmer_ids): # This removes N content
                    if i in sus: # This is where we actually delete the N content, in case that is eventually supported.
                        kmer_ids[i] = None
                        reads[i] = None
                        starts[i] = None
                        reverses[i] = None
                kmer_ids = list(filter(lambda k: k is not None, kmer_ids)) 
                reads = list(filter(lambda r: r is not None, reads))
                starts = list(filter(lambda s: s is not None, starts))
                reverses = list(filter(lambda r: r is not None, reverses))
            else:
                logger.debug("Checking each k-mer for IUPAC substitution or N content.")
                for i, x in enumerate(kmer_ids): # Here we remove the k-mer ids where N-content is detected, in case they are needed, you can use kmer_ids prior to this point to build functionality.
                    if i in sus:
                        kmer_ids[i] = None
                logger.info("Eliminating suspicious 'sus' k-mers, i.e. those with N-content or IUPAC substitutions")
                kmer_ids = list(filter(lambda k: k is not None, kmer_ids))
                reads = [] # I'm keeping this in, just in case for some reason the variable names are needed in the 
                starts = []
                reverses = []
                if None in kmer_ids:
                    logger.debug("In the no-metadata field")
                    # Actually they were just introduced to be filtered out, instead of deleted
                    # Because each deletion whould cange the array index
                    # So instead we set ghtme to None, and filter out
                    raise ValueError("K-mer ids should never actually be none.")

            #logger.debug(kmer_ids)
            logger.debug("{0} 'clean' kmers were identified successfully from {1} input sequences".format(len(kmer_ids), num_recs))
                
            

            logger.debug("Flatmapped {0} kmers for their metadata aggregation".format(len(kmer_ids), len(starts)))
            # Assert that all list lengths are equal before adding metadata to k-mers
            if all_metadata is True and len(kmer_ids) == len(reads) and len(reads) == len(starts) and len(starts) == len(reverses):
                N = len(starts)


                
                data = list(zip(kmer_ids, reads, starts, reverses))
                # Everything is in the right order
                # logger.debug("num k-mer ids: {0}".format(len(kmer_ids)))
                # logger.debug("K-mer id types: {0}".format(type(kmer_ids[0])))
                # logger.debug("Example: {0}".format(kmer_ids[0]))
                # logger.debug("Num reads: {0}".format(len(reads)))
                # logger.debug("reads type: {0}".format(type(reads[0])))
                # logger.debug("Example: {0}".format(reads[0]))
                # logger.debug("Num reverses: {0}".format(len(reverses)))
                # logger.debug("reverse type: {0}".format(type(reverses[0])))
                # logger.debug("Example: {0}".format(reverses[0]))
                # raise RuntimeError("Deciding whether to set a dictionary, or a 4x? array")
                logger.debug("Appended {0} records to rows".format(N))

                rows += data 

                logger.warning("Dumping all metadata into the .kdb file eventually. This could be expensive...")

            elif not all_metadata and len(reads) == 0 and len(starts) == 0 and len(reverses) == 0:
                pass # If we're not doing metadata, don't do it
            else: # Raise an error if the numbers of items per list are not equal
                logger.error("{0} kmer ids".format(len(kmer_ids)))
                logger.error("{0} sequence/read associations".format(len(reads)))
                logger.error("{0} start positions for the associations found".format(len(starts)))
                logger.error("{0} reverse bools for each association".format(len(reverses)))
                
                raise ValueError("Unexpectedly, the number of ids did not match up with the number of other metadata elements per k-mer OR other unknown error")


            # else:
            #     raise RuntimeError("Still have no clue what's going on...")
            # On disk k-mer counting
            # Thank you, this was brilliant
            # https://stackoverflow.com/a/9294062/12855110
            num_kmers = len(kmer_ids)

            if num_kmers == 0:
                raise ValueError("No k-mers to add. Something likely went wrong. Please report to the issue tracker")
            else:
                # for kmer_id in kmer_ids:
                #     db.conn.execute("UPDATE {0} SET count = count + 1 WHERE id = $1".format(db._tablename),
                inserts = list(map(lambda x: (x+1,), kmer_ids))

                batches = ceil(num_kmers/rows_per_batch)
                t0 = datetime.now()
                for i in range(batches):

                    try:
                        with db.conn.begin():

                            db.conn.execute("UPDATE {0} SET count = count + 1 WHERE id IN (%s)".format(db._tablename), inserts[i*rows_per_batch:(i*rows_per_batch+rows_per_batch)])

                    except Exception as e:
                        logger.error(db.conn)
                        raise e

                    t1 = datetime.now()
                    d = (t1 - t0).seconds
                    logger.info("\n\nLoaded {0} k-mers in batch {1} transaction ({2} rows/second)".format(rows_per_batch, i, rows_per_batch/(d+0.1)))
                    t0 = t1
                    t1 = None
                    
            
            recs = [r for r in seqprsr] # The next block of exactly 'b' reads
            # This will be logged redundantly with the sys.stderr.write method calls at line 141 and 166 of seqparser.py (in the _next_fasta() and _next_fastq() methods)
            #sys.stderr("\n")
            logger.info("Read {0} more records from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"))
        if all_metadata:
            sys.stderr.write("Writing all metadata keys for each k-mer's relationships to reads into the SQLite3 database...\n")

            unique_kmer_ids = list(set(list(map(lambda x: x[0], rows))))



            # If we round up to the nearest page, we should have 'pages' number of pages
            # and we'd read n on each page.
            kid = 0
            Session = sessionmaker(bind=db._engine)
            session = Session()
            class Kmer(Base):
                __tablename__ = db._tablename

                id = Column(Integer, primary_key=True)
                count = Column(Integer)
                starts = Column(JSON)
                reverses = Column(JSON)
                seqids = Column(JSON)

            for i in range(len(unique_kmer_ids)):
                kid = kmer_ids[i] # FOr SQLlite 1-based indexing
                logger.debug("Beginning to commit {0} to the database".format(kid))
                sys.stderr.write("\n")
                kmers = [x for x in rows if x[0] == kid]
                logger.debug("Located {0} relationships involving k-mer {1}".format(len(kmers), kid))
                logger.debug("        === M E S S A G E ===")
                logger.debug("=====================================")
                logger.debug("beginning to process all records of the {0} k-mer".format(kid))


                row = session.query(Kmer).filter_by(id=kid).first()
                if row is None:
                    logger.error(kid)
                    raise RuntimeError("Could not locate k-mer with id {0} in the Postgres table".format(kid))


                logger.debug("Before: {0}".format(len(row.reverses)))
                row.seqids = list(map(lambda y: y[1], kmers))
                row.starts = list(map(lambda y: y[2], kmers))
                row.reverses = list(map(lambda y: y[3], kmers))
                #flag_modified(Kmer, 'seqids')
                #flag_modified(Kmer, 'starts')
                #flag_modified(Kmer, 'reverses')


                logger.debug("After: {0}".format(len(row.reverses)))
                session.add(row)

                if i % rows_per_batch == 0:
                    session.commit()
                    logger.info("Transaction completed for the {0} kmer.".format(kid+1))
            session.commit()
            logger.debug("Sleeping...")
            time.sleep(200)
            session.commit()

            logger.debug("===================================")
            logger.debug("Example record with metadata:")
            result = session.query(Kmer).filter_by(id=kid+1).first()
            logger.debug(result.__dict__)
            logger.debug("===================================")
            session.close()
        seqprsr.nullomers = db._get_nullomers() # Calculate nullomers at the end
        seqprsr.total_kmers = db._get_sum_counts() # The total number of k-mers processed
        # Get nullomer ids
        res = list(db.conn.execute("SELECT id FROM {0} WHERE count = 0".format(db._tablename)).fetchall())
        if type(res) is not list or not all(type(x[0]) is int for x in res):
            logger.error("{0}\nSELECT id FROM kmers WHERE count = 0".format(res))
            logger.error("Type of res: {0}".format(type(res)))
            logger.error("Types of values: {0}".format([type(x[0]) for x in res]))
            logger.error("The result of the query was not a list of singleton tuples of ints")
            raise ValueError("PostgreSQL database query returned unexpected data types")
        nullomers = list(map(lambda x: x[0], res))

    finally:
        sys.stderr.write("\n\n\nFinished counting k-mers{0} from '{1}' into the database...\n\n\n".format(' and metadata' if all_metadata else '', filepath))
        if db.conn is not None:
            db.conn.close()

    return db._tablename, seqprsr.header_dict(), nullomers


# class MetadataAppender:
#     def __init__(self, ids:list, reads:list, starts:list, reverses:list, data:dict={}):
#         if type(ids) is not list or not all(type(i) is int for i in ids):
#             raise TypeError("kmerdb.parse.MetadataAppender.__init__() expects a list of ints as its first positional argument")
#         elif type(reads) is not list or not all(type(r) is str for r in reads):
#             raise TypeError("kmerdb.parse.MetadataAppender.__init__() expects a list of strings as its second positional argument")
#         elif type(starts) is not list or not all(type(s) is int for s in starts):
#             raise TypeError("kmerdb.parse.MetadataAppender.__init__() expects a list of ints as its third positional argument")
#         elif type(reverses) is not list or not all(type(r) is bool for r in reverses):
#             raise TypeError("kmerdb.parse.MetadataAppender.__init__() expects a list of bools as its fourth positional argument")
#         elif type(data) is not dict:
#             raise TypeError("kmerdb.parse.MetadataAppender.__init__{} expects a data dictionary as its fifth positional argument")

#         self.data = data
#         self.ids = ids
#         self.reads = reads
#         self.starts = starts
#         self.reverses = reverses

#     def appendRow(self, i):
#         """
#         :param i: 
#         :type i: int
#         :returns: 
#         :rtype: tuple

#         """
#         if type(i) is not int:
#             raise TypeError("kmerdb.MetadataAppender.appendRow expects an int as its first positional argument")
#         logger.info("Processing the {0} k-mer by appending read relations, start offsets, and reverse 3-tuples to the metadata column".format(i))
#         row = (self.ids[i], self.reads[i], self.starts[i], self.reverses[i])
#         def helper(row, data, j):
#             data = yaml.safe_load(data)
#             if type(data) is str:
#                 logger.error(data)
#                 raise ValueError("data: '{0}' with j={1} parsed as a string".format(data, j))
#             elif type(data) is list:
#                 data.append(row[j])
#             else:
#                 logger.debug(data)
#                 raise RuntimeError("data cannot parse as something other than list...")
#             return data
        
#         # with self.conn.begin():
#         #     seqids_data = self.conn.execute("SELECT seqids FROM kmers WHERE id = ?", self.ids[i])
#         #     seqids = helper(row, seqids_data, 1)
#         #     self.conn.execute("UPDATE kmers SET seqids = ? WHERE id = ?", json.dumps(seqids), self.ids[i])


#         # with self.conn.begin():
#         #     starts_data = self.conn.execute("SELECT starts FROM kmers WHERE id = ?", self.ids[i])
#         #     starts = helper(row, starts_data, 2)
#         #     self.conn.execute("UPDATE kmers SET starts = ? WHERE id = ?", json.dumps(starts), self.ids[i])


#         # with self.conn.begin():
#         #     reverses_data = self.conn.execute("SELECT reverses FROM kmers WHERE id = ?", self.ids[i])
#         #     reverses = helper(row, reverses_data, 3)
#         #     self.conn.execute("UPDATE kmers SET reverses = ? WHERE id = ?", json.dumps(reverses), self.ids[i])

#         return row



class Parseable:
    def __init__(self, arguments):
        self.arguments = arguments
            
        
    def parsefile(self, filename):
        """Wrapper function for parse.parsefile to keep arguments succinct for deployment through multiprocessing.Pool
            
        :param data: { 'filename': ..., 'args': { argparse } }
        :type data: dict
        :returns: (db, m, n)
        """
        return parsefile(filename, self.arguments.k, self.arguments.postgres_connection, p=self.arguments.parallel_fastq, rows_per_batch=self.arguments.batch_size, b=self.arguments.fastq_block_size, n=self.arguments.n, stranded=self.arguments.strand_specific, all_metadata=self.arguments.all_metadata)


