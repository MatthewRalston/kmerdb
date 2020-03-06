import os
import sys
import yaml
from itertools import chain
import tempfile

#import threading
from multiprocessing import Pool


from kdb import seqparser, database, kmer

import logging
logger = logging.getLogger(__file__)


def parsefile(filepath, k, p=1, b=50000, stranded=True):
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
        raise TypeError("kdb.parse.parsefile expects a str as its first positional argument")
    elif not os.path.exists(filepath):
        raise OSError("kdb.parse.parsefile could not find the file '{0}' on the filesystem".format(filepath))
    elif type(k) is not int:
        raise TypeError("kdb.parse.parsefile expects an int as its second positional argument")
    elif type(p) is not int:
        raise TypeError("kdb.parse.parsefile expects the keyword argument 'p' to be an int")
    elif type(b) is not int:
        raise TypeError("kdb.parse.parsefile expects the keyword argument 'b' to be an int")
    elif type(stranded) is not bool:
        raise TypeError("kdb.parse.parsefile expects the keyword argument 'stranded' to be a bool")
    # Create temporary SQLite3 database file for on-disk k-mer counting
    temp = tempfile.NamedTemporaryFile(mode="w+", suffix=".sqlite3", delete=False)
    logger.debug("Creating temporary database to tally k-mers: '{0}'".format(temp.name))
    temp.close()
    db = database.SqliteKdb(temp.name, k)

    try:
        # Build fasta/fastq parser object to stream reads into memory
        logger.debug("Constructing fastq parser object...")
        seqprsr = seqparser.SeqParser(filepath, b, k)
        
        pool = Pool(processes=p) # A multiprocessing pool of depth 'p'
        Kmer = kmer.Kmers(k, strand_specific=stranded) # A wrapper class to shred k-mers with
        recs = [r for r in seqprsr] # A block of exactly 'b' reads-per-block to process in parallel
        while len(recs): # While the seqprsr continues to produce blocks of reads
            # Run each read through the shred method
            list_of_dicts = pool.map(Kmer.shred, recs)
            # Flatmap to 'kmers', the dictionary of {'id': read_id, 'kmers': [ ... ]}
            kmer_ids = list(chain.from_iterable(map(lambda x: x['kmers'], list_of_dicts)))
            #read_kmer_relations = list(chain.from_iterable(map(lambda x: list(zip(itertools.repeat(x['id']), x['kmers'])), list_of_dicts)))
            # On disk k-mer counting
            # https://stackoverflow.com/a/9294062/12855110
            db.conn.execute("UPDATE kmers SET count = count + 1 WHERE id = ?",
                            list(map(lambda x: (x+1,), kmer_ids)))
            #db.conn.execute("INSERT INTO reads(read_id, kmer_id) VALUES (?, ?)", read_kmer_relations)
            recs = [r for r in seqprsr] # The next block of exactly 'b' reads

        seqprsr.nullomers = db._get_nullomers() # Calculate nullomers at the end
        seqprsr.total_kmers = db._get_sum_counts() # The total number of k-mers processed
    finally:
        sys.stderr.write("\n")
        logger.info("Finished loading records from '{0}' into '{1}'...".format(filepath, temp.name))
        #db._engine.dispose()
    return db, seqprsr.header_dict()


