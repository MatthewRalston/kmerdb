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

import logging
logger = logging.getLogger(__file__)





def parsefile(filepath:str, k:int, ): #rows_per_batch:int=100000, b:int=50000, n:int=1000, both_strands:bool=False, all_metadata:bool=False):
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
    :raise TypeError: rows_per_batch was invalid
    :raise TypeError: b was invalid
    :raise TypeError: n was invalid
    :raise TypeError: both_strands was invalid
    :raise TypeError: all_metadata was invalid
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
    data = {} # This is the dictionary of tuples, keyed on k-mer id, and containing 3-tuples ('kmer_id', 'read/start/reverse')
    keys = set()
    rows = []
    N = 4**k
    nullomers = set()
    # Build fasta/fastq parser object to stream reads into memory
    logger.debug("Initializing parser...")
    seqprsr = SeqParser(filepath, k)
    fasta = not seqprsr.fastq # Look inside the seqprsr object for the type of file

    # Initialize the kmer array
    try:
        counts = np.zeros(N, dtype="uint64")
    except TypeError as e:
        logger.error("Invalid dtype for numpy array instantiation")
        logger.error(e)
        raise e

    logger.info("Successfully allocated space for {0} unsigned integers: {1} bytes".format(N, counts.nbytes))
        # Instantiate the kmer class
    Kmer = kmer.Kmers(k, verbose=fasta) # A wrapper class to shred k-mers with

    recs = [r for r in seqprsr] # A block of exactly 'b' reads-per-block to process in parallel
    if not fasta:
        logger.debug("Read exactly {0} records from the seqparser object".format(len(recs)))
        assert len(recs) <= b, "The seqparser should return exactly {0} records at a time".format(b)
    else:
        logger.debug("Skipping the block size assertion for fasta files")
    logger.info("Read {0} sequences from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"))
    all_kmer_metadata = list([] for x in range(N))
    while len(recs): # While the seqprsr continues to produce blocks of reads

        num_recs = len(recs)
        logger.info("Processing a block of {0} reads/sequences".format(num_recs))

        # Initialize an all_kmer_metadata list
        kmer_metadata = []
        list_of_dicts = list(map(Kmer.shred, recs))

        # Flatmap to 'kmers', the dictionary of {'id': read_id, 'kmers': [ ... ]}
        logger.debug("\n\nAcquiring linked-list of all k-mer ids from {0} sequence records...\n\n".format(num_recs))
        kmer_ids = list(chain.from_iterable(map(lambda x: x['kmers'], list_of_dicts)))
        logger.debug("Successfully allocated linked-list for all k-mer ids read.")
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
                raise ValueError("kmerdb.parse.parsefile encountered an invalid kmer_id. Internal error.")

        #logger.debug(kmer_ids)
        logger.debug("{0} 'clean' kmers were identified successfully from {1} input sequences".format(len(kmer_ids), num_recs))
        logger.debug("Flatmapped {0} kmers for their metadata aggregation".format(len(kmer_ids), len(starts)))
        # Assert that all list lengths are equal before adding metadata to k-mers
        if all_metadata is True and len(kmer_ids) == len(reads) and len(reads) == len(starts) and len(starts) == len(reverses):
            N = len(starts)
                
            kmer_metadata = list(zip(kmer_ids, reads, starts, reverses))
            # Everything is in the right order
            logger.warning("Dumping all metadata into the .kdb file eventually. This could be expensive...")

        elif not all_metadata and len(reads) == 0 and len(starts) == 0 and len(reverses) == 0:
            logger.debug("Skipping metadata allocation")
            pass # If we're not doing metadata, don't do it
        else: # Raise an error if the numbers of items per list are not equal
            logger.error("{0} kmer ids".format(len(kmer_ids)))
            logger.error("{0} sequence/read associations".format(len(reads)))
            logger.error("{0} start positions for the associations found".format(len(starts)))
            logger.error("{0} reverse bools for each association".format(len(reverses)))
                
            raise ValueError("Unexpectedly, the number of ids did not match up with the number of other metadata elements per k-mer OR other unknown error. Internal error.")

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
            for kmer in kmer_ids:
                counts[kmer] += 1
        # all_kmer_metadata
        if all_metadata:
            for single_kmer_id, read, start, reverse in kmer_metadata:
                all_kmer_metadata[single_kmer_id].append((read, start, reverse))

        recs = [r for r in seqprsr] # The next block of exactly 'b' reads
        logger.info("Read {0} more records from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"))
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


    
    seqprsr.nullomer_array = np.array(all_theoretical_kmer_ids, dtype="uint64")[is_nullomer]
    sys.stderr.write("\n\n\nFinished counting k-mers{0} from '{1}'...\n\n\n".format(' and metadata' if all_metadata else '', filepath))

    return counts, seqprsr.header_dict(), seqprsr.nullomer_array, all_kmer_metadata





class SeqParser:
    """
    Largely useless module, needs 3 pieces of information passed back in from the outside.
    This performs ugly decompression of fasta and fastq files, patching the __next__ methods, effectively.
    It allows you to read either fasta or fastq data in blocks, obviously useful for the latter.


    :ivar filepath: The .fastq, .fastq.gz, .fasta, .fasta.gz, .fna, .fna.gz, .fa, or .fa.gz file.
    :ivar num: The number of records to read in from a .fastq
    :ivar k: The choice of k to initialize the calculation of kmer/nullomer counts.
    """
    def __init__(self, filepath:str, k:int, num:int=100000, ):
        """
        The SeqParser class wraps up some functionality  
        """
        
        if type(filepath) is not str:
            raise TypeError("kmerdb.parse.SeqParser expects a str as its first positional argument")
        elif not os.access(filepath):
            raise TypeError("kmerdb.parse.SeqParser expects an existing path on the operating system as its first argument")
        elif type(num) is not int:
            raise TypeError("kmerdb.parse.SeqParser expects an int as its second positional argument")
        elif type(k) is not int:
            raise TypeError("kmerdb.parse.SeqParser expects an int as its third positional argument")
        
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
                logger.info("Opening gzipped fastq file '{0}'...".format(filepath))
                self._handle = SeqIO.parse(gzip.open(self.filepath, 'rt'), "fastq")
            else:
                logger.info("Opening gzipped fasta file '{0}'...".format(filepath))
                self._handle = SeqIO.parse(gzip.open(self.filepath, 'rt'), "fasta")
        else:
            if self.fastq:
                logger.info("Opening uncompressed fastq file '{0}'...".format(filepath))
                self._handle = SeqIO.parse(open(self.filepath, 'r'), "fastq")
            else:
                logger.info("Opening uncompressed fasta file '{0}'...".format(filepath))
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
            logger.error(e)
            logger.error("\n\nFastq format error: '{0}' seems to not be fastq format\n\n")
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
    def __init__(self, arguments):
        self.arguments = arguments
            
        
    def parsefile(self, filename):
        """Wrapper function for parse.parsefile to keep arguments succinct for deployment through multiprocessing.Pool
            
        :param filename: the filepath of the fasta(.gz)/fastq(.gz) to process with parsefile -> parse.SeqParser
        :type filename: str
        :returns: (db, m, n)
        """
        return parsefile(filename, self.arguments.k)


