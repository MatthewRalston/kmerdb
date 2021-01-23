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
from itertools import chain, repeat
from concurrent.futures import ThreadPoolExecutor


import tempfile

#import threading
from multiprocessing import Pool


from kmerdb import seqparser, database, kmer

import logging
logger = logging.getLogger(__file__)


def parsefile(filepath:str, k:int, p:int=1, b:int=50000, n:int=1000, stranded:bool=True, all_metadata:bool=False):
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
    elif type(p) is not int:
        raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'p' to be an int")
    elif type(b) is not int:
        raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'b' to be an int")
    elif type(n) is not int:
        raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'n' to be an int")
    elif type(stranded) is not bool:
        raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'stranded' to be a bool")
    # Create temporary SQLite3 database file for on-disk k-mer counting
    temp = tempfile.NamedTemporaryFile(mode="w+", suffix=".sqlite3", delete=False)
    logger.debug("Creating temporary database to tally k-mers: '{0}'".format(temp.name))
    temp.close()
    db = database.SqliteKdb(temp.name, k)
    from kmerdb import kmer


    
    data = {} # This is the dictionary of tuples, keyed on k-mer id, and containing 3-tuples ('kmer_id', 'read/start/reverse')
    keys = set()
    rows = []
    try:
        # Build fasta/fastq parser object to stream reads into memory
        logger.debug("Constructing SeqParser object...")
        seqprsr = seqparser.SeqParser(filepath, b, k)
        fasta = not seqprsr.fastq
        logger.debug("Constructing multiprocessing pool with {0} processors".format(p))
        pool = Pool(processes=p) # A multiprocessing pool of depth 'p'
        Kmer = kmer.Kmers(k, strand_specific=stranded, fasta=fasta) # A wrapper class to shred k-mers with
        # Look inside the seqprsr object for the type of file


        recs = [r for r in seqprsr] # A block of exactly 'b' reads-per-block to process in parallel
        if not fasta:
            logger.debug("Read exactly b=={0}=={1} records from the {2} seqparser object".format(b, len(recs), s))
            assert len(recs) == b, "The seqparser should return exactly {0} records at a time".format(b)
        else:
            logger.debug("Skipping the block size assertion for fasta files")
        logger.info("Read {0} sequences from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"))

        while len(recs): # While the seqprsr continues to produce blocks of reads
            # Run each read through the shred method
            list_of_dicts = pool.map(Kmer.shred, recs)

            logger.info("Shredded up {0} sequences over {1} parallel cores, like a cheesesteak".format(len(list_of_dicts), p))
            #logger.debug("Everything in list_of_dicts is perfect, everything past here is garbage")
            
            # Flatmap to 'kmers', the dictionary of {'id': read_id, 'kmers': [ ... ]}
            kmer_ids = list(chain.from_iterable(map(lambda x: x['kmers'], list_of_dicts)))
            reads = list(chain.from_iterable(map(lambda x: x['seqids'], list_of_dicts)))
            starts = list(chain.from_iterable(map(lambda x: x['starts'], list_of_dicts)))
            reverses = list(chain.from_iterable(map(lambda x: x['reverses'], list_of_dicts)))
            logger.debug("Flatmapped {0} kmers for their metadata aggregation".format(len(kmer_ids), len(starts)))
            # Assert that all list lengths are equal before adding metadata to k-mers
            if all_metadata is True and len(kmer_ids) == len(reads) and len(reads) == len(starts) and len(starts) == len(reverses):
                N = len(starts)

                data = list(zip(kmer_ids, reads, starts, reverses))

                logger.debug("Appended {0} records to rows".format(N))

                rows += data 

                logger.warning("Dumping all metadata into the .kdb file eventually. This could be expensive...")

            elif len(kmer_ids) == len(reads) and len(reads) == len(starts) and len(starts) == len(reverses) and not all_metadata:
                pass # If we're not doing metadata, don't do it
            else: # Raise an error if the numbers of items per list are not equal
                logger.error(len(kmer_ids))
                logger.error(len(reads))
                logger.error(len(start))
                logger.error(len(reverses))
                
                raise ValueError("Unexpectedly, the number of ids did not match up with the number of other metadata elements per k-mer OR other unknown error")


            # else:
            #     raise RuntimeError("Still have no clue what's going on...")
            # On disk k-mer counting
            # Thank you, this was brilliant
            # https://stackoverflow.com/a/9294062/12855110
            logger.debug("Parsed/mapped the remainder or the data from list_of_dicts")
            logger.info("Updating the k-mer counts in the SQLAlchemy connection to SQLite3 for {0} k-mer ids in one transation".format(len(kmer_ids)))
            with db.conn.begin():
                db.conn.execute("UPDATE kmers SET count = count + 1 WHERE id = ?",
                                list(map(lambda x: (x+1,), kmer_ids)))
            
            recs = [r for r in seqprsr] # The next block of exactly 'b' reads
            # This will be logged redundantly with the sys.stderr.write method calls at line 141 and 166 of seqparser.py (in the _next_fasta() and _next_fastq() methods)
            #sys.stderr("\n")
            logger.info("Read {0} more records from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"))
        if all_metadata:
            sys.stderr.write("Writing all metadata keys for each k-mer's relationships to reads into the SQLite3 database...\n")

            unique_kmer_ids = list(set(list(map(lambda x: x[0], rows))))

            def submit(conn, kmers, kid):
                logger.debug("        === M E S S A G E ===")
                logger.debug("=====================================")
                logger.debug("beginning to process all records of the {0} k-mer".format(kid))
                with conn.begin():
                    conn.execute("UPDATE kmers SET seqids = ?, starts = ?, reverses = ? WHERE id = ?", json.dumps(list(map(lambda y: y[1], kmers))), json.dumps(list(map(lambda y: y[2], kmers))), json.dumps(list(map(lambda y: y[2], kmers))), kid+1)

                return kid

            # If we round up to the nearest page, we should have 'pages' number of pages
            # and we'd read n on each page.
            kid = 0
            for i in range(len(unique_kmer_ids)):
                kid = unique_kmer_ids[i] # FOr SQLlite 1-based indexing
                logger.debug("Beginning to commit {0} to the database".format(kid))
                sys.stderr.write("\n")
                kmers = [x for x in rows if x[0] == kid]
                logger.debug("Located {0} relationships involving k-mer {1}".format(len(kmers), kid))
                submit(db.conn, kmers, kid)

                    
                logger.info("Transaction completed for the {0} kmer.".format(kid))

            logger.debug("===================================")
            logger.debug("Example record with metadata:")
            result = db.conn.execute("SELECT * FROM kmers WHERE id = ?", kid).fetchone()
            logger.debug(result)
            logger.debug("===================================")

        seqprsr.nullomers = db._get_nullomers() # Calculate nullomers at the end
        seqprsr.total_kmers = db._get_sum_counts() # The total number of k-mers processed
    finally:

        sys.stderr.write("\n\n\nFinished counting k-mers{0} from '{1}' into '{2}'...\n\n\n".format(' and metadata' if all_metadata else '', filepath, temp.name))

    return db, seqprsr.header_dict()


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
