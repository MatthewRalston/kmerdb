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

import io
import os
import sys
import logging
import gzip
import hashlib
import yaml
import json
import time
from datetime import datetime
from math import ceil
from itertools import chain, repeat


from Bio import SeqIO


logger = logging.getLogger(__file__)

# from sqlalchemy.orm import sessionmaker
# from sqlalchemy.orm.attributes import flag_modified
# from sqlalchemy.ext.declarative import declarative_base
# from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey, Sequence, JSON, Boolean

#from psycopg2 import sql
import tempfile
import numpy as np


from kmerdb import kmer, util


def parse_sequence_file(seq_filepath:str, return_tuple:bool=True):
    """
    Return 2-tuples of (seq_id, seq) strings from a fastafile (NA or AA)
    
    """
    if type(seq_filepath) is not str:
        raise TypeError("kmerdb.graph.parse_sequence_file() expects a fasta/fastq sequence filepath as a str as its only positional_regument")

    if os.path.exists(seq_filepath) is False or os.access(seq_filepath, os.R_OK) is False:
        raise ValueError("kmerdb.graph.parse_sequence_file() expects the filepath to be be readable on the filesystem")

    try:
        logger.debug("Beginning to process sequence filepath '{0}'".format(seq_filepath))
        if util.is_gz_file is True:
            seqhandle = gzip.open(seq_filepath, "rb")
        else:
            seqhandle = open(seq_filepath, 'r')
        if util.is_fasta(seq_filepath) is True:
            parser = SeqIO.parse(seqhandle, "fasta")
        elif util.is_fastq(seq_filepath) is True:
            parser = SeqIO.parse(seqhandle, "fastq")
        else:
            raise ValueError("Could not determine the format of file '{0}'".format(seq_filepath))
        for s in parser: # s is a Bio.SeqRecord
            seq = str(s.seq)
            seq_id = str(s.id)
            if return_tuple is True:
                yield (seq_id, seq)
            else:
                yield s
    except Exception as e:
        raise e
    finally:
        seqhandle.close()




def parsefile(filepath:str, k:int, replace_with_none:bool=True): 
    """Parse a single sequence file in blocks/chunks with multiprocessing support

    :param filepath: Path to a fasta or fastq file
    :type filepath: str
    :param k: Choice of k to shred k-mers with
    :type k: int
    :raise TypeError: filepath was invalid
    :raise OSError: filepath was invalid
    :raise TypeError: k was invalid
    :raise ValueError: invalid (None) kmer id detected
    :raise ValueError: mismatched number of kmer_ids, associated sequence/read ids, starting locations, and reverse bools
    :raise AssertionError: Error in nullomer count estimation
    :returns: (counts, header_dictionary, nullomer_array, all_metadata) header_dictionary is the file's metadata for the header block
    :rtype: (numpy.ndarray, dict, list, list)

    """

    from kmerdb import kmer
    if filepath is None or type(filepath) is not str:
        raise TypeError("kmerdb.parse.parsefile expects a str as its first positional argument")
    elif not os.path.exists(filepath):
        raise OSError("kmerdb.parse.parsefile could not find the file '{0}' on the filesystem".format(filepath))
    elif k is None or type(k) is not int:
        raise TypeError("kmerdb.parse.parsefile expects an int as its second positional argument")
    elif type(replace_with_none) is not bool:
        raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'replace_with_none' to be a bool")
    N = 4**k
    seq_lengths = []
    total_kmers = 0
    counts = np.zeros(N, dtype="uint64")
    nullomers = set()
    
    final_kmer_ids = []
    md5, sha256 = util.checksum(filepath)



    for seq in parse_sequence_file(filepath, return_tuple=False):
        seq_id = seq.id
        seqlen = len(seq)        
        kmer_ids, seq_ids, pos = kmer.shred(seq, k, replace_with_none=replace_with_none, quiet_iupac_warning=False)

        for kmer_id in kmer_ids:
            if kmer_id is not None:
                counts[kmer_id] += 1
                total_kmers += 1
        seq_lengths.append(seqlen)

    is_nullomer = np.where(counts == 0)
    nullomer_array = np.array(range(N), dtype="uint64")[is_nullomer]
    unique_kmers = int(np.count_nonzero(counts))
    
    num_nullomers = N - unique_kmers
    max_read_length = max(seq_lengths)
    min_read_length = min(seq_lengths)
    avg_read_length = int(np.mean(np.array(seq_lengths)))
    assert num_nullomers == len(nullomer_array), "Internal Error: kmerdb.parse.parsefile() found inconsistencies between two ways of counting nullomers. Internal error."
    
    file_metadata = {
        "filename": filepath,
        "md5": md5,
        "sha256": sha256,
        "total_reads": len(seq_lengths),
        "total_kmers": total_kmers,
        "unique_kmers": unique_kmers, # or N - len(nullomers),
        "nullomers": num_nullomers,
        "min_read_length": min_read_length,
        "max_read_length": max_read_length,
        "avg_read_length": avg_read_length
    }

    logger.info("\n\n\nFinished counting k-mers from '{0}'...\n\n\n".format(filepath))
    return counts, file_metadata, nullomer_array


