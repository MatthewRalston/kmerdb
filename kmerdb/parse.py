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

import tempfile

#import threading
from multiprocessing import Pool


from kmerdb import seqparser, database, kmer

import logging
logger = logging.getLogger(__file__)


def parsefile(filepath:str, k:int, p:int=1, b:int=50000, stranded:bool=True, all_metadata:bool=False):
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
    try:
        # Build fasta/fastq parser object to stream reads into memory
        logger.debug("Constructing SeqParser object...")
        seqprsr = seqparser.SeqParser(filepath, b, k)
        logger.debug("Constructing multiprocessing pool with {0} processors".format(p))
        pool = Pool(processes=p) # A multiprocessing pool of depth 'p'
        Kmer = kmer.Kmers(k, strand_specific=stranded) # A wrapper class to shred k-mers with
        # Look inside the seqprsr object for the type of file
        s = "fastq" if seqprsr.fastq else "fasta"

        recs = [r for r in seqprsr] # A block of exactly 'b' reads-per-block to process in parallel
        if s == "fastq":
            logger.debug("Read exactly b=={0}=={1} records from the {2} seqparser object".format(b, len(recs), s))
            assert len(recs) == b, "The seqparser should return exactly {0} records at a time".format(b)
        else:
            logger.debug("Read {0} sequences from the {1} seqparser object".format(len(recs), s))
            logger.debug("Skipping the block size assertion for fasta files")
        while len(recs): # While the seqprsr continues to produce blocks of reads
            # Run each read through the shred method
            list_of_dicts = pool.map(Kmer.shred, recs)

            logger.info("Shredding up {0} sequences over {1} parallel cores, like a cheesesteak".format(len(list_of_dicts), p))




            
            # Flatmap to 'kmers', the dictionary of {'id': read_id, 'kmers': [ ... ]}
            kmer_ids = list(chain.from_iterable(map(lambda x: x['kmers'], list_of_dicts)))
            logger.debug("Flatmapped {0} kmer ids for these {1} sequence ids".format(len(kmer_ids), len(list_of_dicts)))
            read_kmer_relations = list(chain.from_iterable(map(lambda x: x['seqids'], list_of_dicts)))
                                       
            logger.debug("Flatmapped and matched {0} kmers with these {1} sequence ids".format(len(read_kmer_relations), len(list_of_dicts)))


            read_kmer_start_offsets = list(chain.from_iterable(map(lambda x: x['starts'], list_of_dicts)))
            logger.debug("Flatmapped {0} kmers for their {1} offsets".format(len(kmer_ids), len(read_kmer_start_offsets)))

            kmer_reverses = list(chain.from_iterable(map(lambda x: x['reverses'], list_of_dicts)))
            logger.debug("Flatmapped {0} reverse? bools for these {1} k-mers that were shredded".format(len(kmer_reverses), len(kmer_ids)))
            if all_metadata is True and len(kmer_ids) == len(read_kmer_relations) and len(read_kmer_relations) == len(read_kmer_start_offsets) and len(read_kmer_start_offsets) == len(kmer_reverses):
                N = len(read_kmer_start_offsets)
                for i in range(N):
                    logger.info("Processing the {0} k-mer by appending read relations, start offsets, and reverse 3-tuples to the metadata column".format(i))                    
                    if kmer_ids[i] in keys:
                        
                        data[kmer_ids[i]].append((kmer_ids[i], read_kmer_relations[i], read_kmer_start_offsets[i], kmer_reverses[i]))
                    else:
                        data[kmer_ids[i]] = [(kmer_ids[i], read_kmer_relations[i], read_kmer_start_offsets[i], kmer_reverses[i])]

                    logger.debug("============================")
                    logger.debug("id")
                    logger.debug(kmer_ids[i])
                    logger.debug("============================")
                    logger.debug("all metadata for this k-mer")
                    logger.debug(data[kmer_ids[i]])
                    logger.debug("read_kmer_relations")
                    logger.debug(read_kmer_relations[i])
                    logger.debug("start offsets")
                    logger.debug(read_kmer_start_offsets[i])
                    logger.debug("reverses")
                    logger.debug(kmer_reverses[i])
                        
                keys = keys.update(set(data.keys()))
                recs = []
                logger.warning("Dumping all metadata into the .kdb file eventually")
                logger.warning("Deferring to dump metadata to SQLite3 later...")
            elif len(kmer_ids) == len(read_kmer_relations) and len(read_kmer_relations) == len(read_kmer_start_offsets) and len(read_kmer_start_offsets) == len(kmer_reverses) and not all_metadata:
                pass
            else:
                logger.error(len(kmer_ids))
                logger.error(len(read_kmer_relations))
                logger.error(len(read_kmer_start_offsets))
                logger.error(len(kmer_reverses))
                
                raise ValueError("Unexpectedly, the number of ids did not match up with the number of other metadata elements per k-mer OR other unknown error")
            # else:
            #     raise RuntimeError("Still have no clue what's going on...")
            # On disk k-mer counting
            # Thank you, this was brilliant
            # https://stackoverflow.com/a/9294062/12855110
            logger.debug("Parsed/mapped the remainder or the data from list_of_dicts")
            logger.info("Updating the k-mer counts in the SQLAlchemy connection to SQLite3 for {0} k-mer ids in one transation".format(len(kmer_ids)))
            db.conn.execute("UPDATE kmers SET count = count + 1 WHERE id = ?",
                            list(map(lambda x: (x+1,), kmer_ids)))
            
            recs = [r for r in seqprsr] # The next block of exactly 'b' reads
            # This will be logged redundantly with the sys.stderr.write method calls at line 141 and 166 of seqparser.py (in the _next_fasta() and _next_fastq() methods)
            #sys.stderr("\n")
            logger.debug("Read {0} more records from the {1} seqparser object".format(len(recs), s))
        if all_metadata:
            sys.stderr.write("Writing all metadata keys for each k-mer's relationships to reads into the SQLite3 database...")
            for kmer_id in data.keys():
                kmer = data[kmer_id]
                # Tuple positions assigned from the tuple in the previous for loop

                # Merge previous entry in database
                reads = db.conn.execute("SELECT seqids FROM kmers WHERE ID = ?", kmer_id)
                # with tuple positions from the previous loop in the global data hash for this k-mer id
                reads = [r["seqids"] for r in reads if r["seqids"] is not None] + data[kmer_id]

                print(reads)
                logger.info("Exiting...")
                logger.debug("Exiting...")
                sys.exit(1)
                
                # Repeat the same strategy for the following
                starts = db.conn.execute("SELECT starts FROM kmers WHERE ID = ?", kmer_id)
                starts = [s["starts"] for s in starts if s["starts"] is not None] + data[kmer_id]
                reverses = db.conn.execute("SELECT reverses FROM kmers WHERE ID = ?", kmer_id)
                reverses = [r["reverses"] for r in reverses if r["reverses"] is not None] + data[kmer_id]
                db.conn.execute("UPDATE kmers SET seqids = ?, starts = ?, reverses = ? WHERE ID = ?", json.dumps(reads), json.dumps(starts), json.dumps(reverses), kmer_id)
        seqprsr.nullomers = db._get_nullomers() # Calculate nullomers at the end
        seqprsr.total_kmers = db._get_sum_counts() # The total number of k-mers processed
    finally:
        sys.stderr.write("\n")
        logger.info("Finished loading records from '{0}' into '{1}'...".format(filepath, temp.name))
        #db._engine.dispose()
    return db, seqprsr.header_dict()


