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
import gzip
import hashlib
import yaml
import json
import time
from datetime import datetime
from math import ceil
from itertools import chain, repeat


from Bio import SeqIO


# from sqlalchemy.orm import sessionmaker
# from sqlalchemy.orm.attributes import flag_modified
# from sqlalchemy.ext.declarative import declarative_base
# from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey, Sequence, JSON, Boolean

#from psycopg2 import sql
import tempfile
import numpy as np


from kmerdb import kmer





def parsefile(filepath:str, k:int, logger=None, quiet=True): #rows_per_batch:int=100000, b:int=50000, n:int=1000, both_strands:bool=False, all_metadata:bool=False):
    """Parse a single sequence file in blocks/chunks with multiprocessing support

    :param filepath: Path to a fasta or fastq file
    :type filepath: str
    :param k: Choice of k to shred k-mers with
    :type k: int
    :param b: Number of reads (per block) to process in parallel
    :type b: int
    :param stranded: Strand specificity argument for k-mer shredding process
    :type stranded: bool
    :raise TypeError: filepath was invalid
    :raise OSError: filepath was invalid
    :raise TypeError: k was invalid
    :raise TypeError: invalid dtype
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
    # elif rows_per_batch is None or type(rows_per_batch) is not int:
    #     raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'rows_per_batch' to be an int")
    # elif b is None or type(b) is not int:
    #     raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'b' to be an int")
    # elif n is None or type(n) is not int:
    #     raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'n' to be an int")
    # elif both_strands is None or type(both_strands) is not bool:
    #     raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'both_strands' to be a bool")
    # elif all_metadata is None or type(all_metadata) is not bool:
    #     raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'all_metadata' to be a bool")
    # if logger is None:
    #     raise ValueError("Invalid logger")

    _loggable = logger is not None
    
    data = {} # This is the dictionary of tuples, keyed on k-mer id, and containing 3-tuples ('kmer_id', 'read/start/reverse')
    keys = set()
    rows = []
    N = 4**k
    nullomers = set()

    
    # Build fasta/fastq parser object to stream reads into memory
    if _loggable:
        logger.log_it("Initializing parser...", "INFO")


    
    seqprsr = SeqParser(filepath, k, logger)
    fasta = not seqprsr.fastq # Look inside the seqprsr object for the type of file


    
    # Initialize the kmer array
    try:
        kmer_ids = []
        counts = np.zeros(N, dtype="uint64")
        read_lengths = []        
    except TypeError as e:
        if _loggable:
            logger.log_it("Invalid dtype for numpy array instantiation", "ERROR")
            logger.log_it(e.__str__(), "ERROR")
        raise e

    Kmer = kmer.Kmers(k, verbose=not quiet) # A wrapper class to shred k-mers with

    recs = [r for r in seqprsr] # A block of exactly 'b' reads-per-block to process in parallel

    if _loggable:
        logger.log_it("Read {0} sequences from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"), "INFO")

    while len(recs): # While the seqprsr continues to produce blocks of reads

        num_recs = len(recs)

        if _loggable:
            logger.log_it("Processing a block of {0} reads/sequences".format(num_recs), "INFO")

        # Initialize an all_kmer_metadata list


        
        list_of_dicts = list(map(Kmer.shred, recs))

        
        min_read_length = 50
        max_read_length = 0

        if quiet is False:
            print("Generating k-mer id array")
        # Flatmap to 'kmers', the dictionary of {'id': read_id, 'kmers': [ ... ]}
            logger.log_it("Aggregating k-mers, counts, read_lengths across reads from a fastq file, or sequences in a fasta from file '{0}'... (NOTE: this may be done several times..)".format(filepath), "INFO")
        
        _kmer_ids = []
        for d in list_of_dicts:

            if d["kmers"] is not None:
                
                _kmer_ids += list(map(kmer.kmer_to_id, d["kmers"]))
                read_lengths.append(d["read_length"])
            if _kmer_ids is None:
                raise RuntimeError("placeholder cannot be none")

            
        if _kmer_ids is not None:
            logger.log_it("The _kmer_ids placeholder has {0} elements.".format(len(_kmer_ids)), "INFO")
            kmer_ids = kmer_ids + _kmer_ids

                    
        seqprsr.max_read_length = max(read_lengths)
        seqprsr.min_read_length = min(read_lengths)

        
        seqprsr.avg_read_length = int(np.mean(np.array(read_lengths, dtype="uint32")))





        
        sus = set()

        
        if _loggable:
            logger.log_it("Removing any eroneous k-mers without an id, cause by 'N' content, creating new gaps in the k-mer profile...", "WARNING")
            logger.log_it("Will substitute IUPAC residues with their respective residues, adding one count for each", "WARNING")
        num_sus = 0
        for i, k in enumerate(kmer_ids):
            if k is None:
                sus.add(i) # This propagates the N removal, by adding those to the set for removal later.
                num_sus += 1

        logger.log_it("Created {0} gaps in the k-mer profile to clean it of N-content".format(num_sus), "INFO")
        if num_sus > 0 and _loggable:
            logger.log_it("{0} uncountable k-mer identified".format(num_sus), "WARNING")
            logger.log_it("Likely caused by the presence of the unspecified nucleotide 'N' in a .fasta file", "WARNING")
            logger.log_it("Will remove from the final committed profile", "WARNING")

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

        if _loggable:
            logger.log_it("Found {0} invalid or un-enumerable k-mers, sequences which have no defined index due to N-conent, i".format(num_sus), "INFO")
            logger.log_it("In other words, one or more characters in your .fasta or .fastq file were 'N' or otherwise contained IUPAC characters that need to be substituted, which would have been given independent k-mer ids.", "DEBUG")
            logger.log_it("We have chosen to omit k-mer with unspecified nucleotides 'N', and these make gaps in our database explicitly.", "DEBUG")
            logger.log_it("Checking each k-mer for IUPAC substitution or N content.", "DEBUG")

            logger.log_it("{0} 'clean' kmers were identified successfully after adding valid k-mer ids from {1} input sequences".format(len(kmer_ids), num_recs), "DEBUG")
            logger.log_it("Flatmapped {0} kmers for their metadata aggregation".format(len(kmer_ids)), "DEBUG")
        # Assert that all list lengths are equal before adding metadata to k-mers
        # else:
        #     raise RuntimeError("Still have no clue what's going on...")
        # On disk k-mer counting
        # Thank you, this was brilliant
        # https://stackoverflow.com/a/9294062/12855110
        # 01-01-2022 This has been removed in favor of in-memory counting
        num_kmers = len(kmer_ids)

        if num_kmers == 0:
            raise ValueError("No k-mers available to add. Please report to the issue tracker")
        else:
            sys.stderr.write("\nAccumulating all k-mers from this set of records...")
            for _ in kmer_ids:
                counts[_] += 1

        recs = [r for r in seqprsr] # The next block of exactly 'b' reads


        if _loggable:
            logger.log_it("Read {0} more records from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"), "INFO")
    # Get nullomers
    # only nullomer counts
    unique_kmers = int(np.count_nonzero(counts))


    all_theoretical_kmer_ids = list(range(N))

    
    
    # FIXME
    num_nullomers = N - unique_kmers
    is_nullomer = np.where(counts == 0)


    assert num_nullomers == len(is_nullomer[0]), "kmerdb.parse module find inconsistencies between two ways of counting nullomers. Internal error."
    
    seqprsr.total_kmers = int(np.sum(counts))
    seqprsr.unique_kmers = unique_kmers
    seqprsr.nullomers = num_nullomers
    header = seqprsr.header_dict()
    header["min_read_length"] = min_read_length
    header["max_read_length"] = max_read_length
    header["avg_read_length"] = int(np.mean(read_lengths))

    seqprsr.nullomer_array = np.array(all_theoretical_kmer_ids, dtype="uint64")[is_nullomer]


    if _loggable:
        logger.log_it("\n\n\nFinished counting k-mers from '{0}'...\n\n\n".format(filepath), "INFO")

    return counts, header, seqprsr.nullomer_array, logger.logs





class SeqParser:
    """
    Largely useless module, needs 3 pieces of information passed back in from the outside.
    This performs ugly decompression of fasta and fastq files, patching the __next__ methods, effectively.
    It allows you to read either fasta or fastq data in blocks, obviously useful for the latter.


    :ivar filepath: The .fastq, .fastq.gz, .fasta, .fasta.gz, .fna, .fna.gz, .fa, or .fa.gz file.
    :ivar num: The number of records to read in from a .fastq
    :ivar k: The choice of k to initialize the calculation of kmer/nullomer counts.
    """
    def __init__(self, filepath:str, k:int, logger=None, num:int=5000000):
        """
        The SeqParser class wraps up some functionality  
        """
        
        if type(filepath) is not str:
            raise TypeError("kmerdb.parse.SeqParser expects a str as its first positional argument")
        elif not os.access(filepath, os.R_OK):
            raise TypeError("kmerdb.parse.SeqParser expects an existing path on the operating system as its first argument")
        elif type(num) is not int:
            raise TypeError("kmerdb.parse.SeqParser expects an int as its second positional argument")
        elif type(k) is not int:
            raise TypeError("kmerdb.parse.SeqParser expects an int as its third positional argument")
        elif logger is None:
            raise ValueError("Invalid logger")

        
        self.k = k
        self.num = num
        self.reads = []
        # Header items
        self.filepath = filepath
        self.md5 = None
        self.sha256 = None
        self.total_reads = 0
        self.total_kmers = 0
        self.unique_kmers = 0
        self.nullomers = 0
        self.nullomer_array = []
        self.compressed = False
        self.fastq = False
        exts = os.path.splitext(filepath)

        if exts[-1] == ".gz":
            self.compressed = True
            nogzexts = os.path.splitext(exts[0])
            if nogzexts[-1] == ".fq" or nogzexts[-1] == ".fastq":
                self.fastq = True
            elif nogzexts[-1] == ".fna" or nogzexts[-1] == ".fasta" or exts[-1] == ".fa":
                self.fastq = False
            else:
                raise ValueError("Cannot parse files of extension '{0}'.\n\nRequires fasta (.fna, .fasta, .fa), fastq (.fq, .fastq), or their gzipped equivalents")
        else: # Must be fasta or fastq uncompressed
            if exts[-1] == ".fq" or exts[-1] == ".fastq":
                self.fastq = True
            elif exts[-1] == ".fna" or exts[-1] == ".fasta" or exts[-1] == ".fa":
                self.fastq = False
            else:
                raise ValueError("Cannot parse files of extension '{0}'.\n\nRequires fasta (.fna, .fasta, .fa), fastq (.fq, .fastq), or their gzipped equivalents")

        # This is a really ugly patch to add appropriate fastq and fasta next behavior.
        if self.fastq:
            self.__class__.__iter__ = self._iter_fastq
            self.__class__.__next__ = self._next_fastq
        else:
            self.__class__.__iter__ = self._iter_fasta
            self.__class__.__next__ = self._next_fasta
                
        if self.compressed:
            if self.fastq:
                sys.stderr.write("Opening gzipped fastq file '{0}'...\n".format(filepath))
                self._handle = SeqIO.parse(gzip.open(self.filepath, 'rt'), "fastq")
            else:
                sys.stderr.write("Opening gzipped fasta file '{0}'...\n".format(filepath))
                self._handle = SeqIO.parse(gzip.open(self.filepath, 'rt'), "fasta")
        else:
            if self.fastq:
                sys.stderr.write("Opening uncompressed fastq file '{0}'...\n".format(filepath))
                self._handle = SeqIO.parse(open(self.filepath, 'r'), "fastq")
            else:
                sys.stderr.write("Opening uncompressed fasta file '{0}'...\n".format(filepath))
                self._handle = SeqIO.parse(open(self.filepath, 'r'), "fasta")
        # Get checksums
        self.md5, self.sha256 = self.__checksum()


    def __checksum(self):
        """Generates md5 and sha256 checksums of a file
        :returns: (md5, sha256)
        :rtype: tuple
        """
        if not os.path.exists(self.filepath):
            raise IOError("kmerdb.parse.SeqParser.__checksum could not find '{}' on the filesystem".format(filepath))
        hash_md5 = hashlib.md5()
        hash_sha256 = hashlib.sha256()
        with open(self.filepath, 'rb') as ifile:
            for chunk in iter(lambda: ifile.read(4096), b""):
                hash_md5.update(chunk)
                hash_sha256.update(chunk)
        return (hash_md5.hexdigest(), hash_sha256.hexdigest())

    def header_dict(self):
        """ Create a header dictionary to convert into YAML to go in the header block(s) of the compression header. Has a schema to be validated, defined in config.py

        :returns: dict
        :rtype: dict

        """
        return {
            "filename": self.filepath,
            "md5": self.md5,
            "sha256": self.sha256,
            "total_reads": self.total_reads,
            "total_kmers": self.total_kmers,
            "unique_kmers": self.unique_kmers or 4**self.k - self.nullomers,
            "nullomers": self.nullomers,
            "min_read_length": 30,
            "max_read_length": 500,
            "avg_read_length": 0
        }
            
    def __exit__(self, exc_type, exc_value, traceback):
        self._handle.close()

    def __enter__(self):
        return self
        
        
    def _iter_fastq(self):
        """A custom iterator method to add to the 'reads' array as iterated upon.
        """
        try:
            for i in range(self.num):
                self.reads.append(next(self._handle))
        except StopIteration as e:
            pass
        except ValueError as e:
            raise e
        return self

    def _next_fastq(self):
        """
        A custom mononucleotide counter
        
        """
        if not len(self.reads):
            raise StopIteration
        else:
            self.total_reads += 1
            read = self.reads.pop()
            return read

    def _iter_fasta(self):
        return self

    def _next_fasta(self):

        seq = next(self._handle)
        self.total_reads += 1
        sys.stderr.write("Read {0} sequences from '{1}'...\n".format(self.total_reads, self.filepath))
        return seq
        







class Parseable:
    def __init__(self, arguments, logger=None):
        self.arguments = arguments
        self.logger = logger
        self.avg_read_length = None
        self.max_read_length = None
        self.min_read_length = None

        
    def parsefile(self, filename):
        """Wrapper function for parse.parsefile to keep arguments succinct for deployment through multiprocessing.Pool
            
        :param filename: the filepath of the fasta(.gz)/fastq(.gz) to process with parsefile -> parse.SeqParser
        :type filename: str
        :returns: (db, m, n)
        """

        counts, header, nullomer_array, logs = parsefile(filename, self.arguments.k, logger=self.logger)


        self.max_read_length = header["max_read_length"]
        self.min_read_length = header["min_read_length"]
        self.avg_read_length = header["avg_read_length"]
        
        
        return counts, header, nullomer_array, logs


