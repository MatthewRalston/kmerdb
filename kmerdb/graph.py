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
import pdb

import sys
import os
import gzip
import tempfile
import yaml, json
from collections import OrderedDict
import math
import re

from builtins import open as _open

import jsonschema
from Bio import SeqIO, Seq, bgzf


import numpy as np

from kmerdb import fileutil, seqparser, kmer, config, util

# Logging configuration
import logging
logger = logging.getLogger(__file__)



def open(filepath, mode="r", metadata=None, slurp:bool=False):
    """
    Opens a file for reading or writing. Valid modes are 'xrwbt'. 
    Returns a lazy-loading KDBGReader object or a KDBGWriter object.
    The data may be force loaded with 'slurp=True'

    :param filepath:
    :type filepath: str
    :param mode:
    :type mode: str
    :param metadata: The file header/metadata dictionary to write to the file.
    :type metadata: dict
    :param slurp: Immediately load all data into KDBGReader
    :type slurp: bool
    :returns: kmerdb.fileutil.KDBGReader/kmerdb.fileutil.KDBGWriter
    :rtype: kmerdb.fileutil.KDBGReader
    """
    if type(filepath) is not str:
        raise TypeError("kmerdb.graph.open expects a str as its first positional argument")
    elif type(mode) is not str:
        raise TypeError("kmerdb.graph.open expects the keyword argument 'mode' to be a str")
    elif ("w" in mode or "x" in mode) and (metadata is None or not isinstance(metadata, OrderedDict)):
        raise TypeError("kmerdb.graph.open expects an additional metadata dictionary")
    elif type(slurp) is not bool:
        raise TypeError("kmerdb.graph.open expects a boolean for the keyword argument 'slurp'")
    modes = set(mode)
    if modes - set("xrwbt") or len(mode) > len(modes):
        raise ValueError("invalid mode: {}".format(mode))


    creating = "x" in modes
    reading  = "r" in modes
    writing  = "w" in modes
    binary   = "b" in modes
    text     = "t" in modes

    if text and binary:
        raise ValueError("can't have text and binary mode at once")
    elif not (creating or reading or writing):
        raise ValueError("must have exactly one or read/write")

    if "r" in mode.lower():
        return KDBGReader(filename=filepath, mode=mode, slurp=slurp)
    elif "w" in mode.lower() or "x" in mode.lower():
        return KDBGWriter(filename=filepath, mode=mode, metadata=metadata)
    else:
        raise ValueError("Bad mode %r" % mode)


    # DO STUFF

class KDBGReader(bgzf.BgzfReader):
    """
    A class to read .kdbg files

    """
    def __init__():

        return

    
class KDBGWriter(bgzf.BgzfWriter):
    """
    A wrapper class around Bio.bgzf.BgzfWriter to write a .kdbg file to disk.

    :ivar metadata: OrderedDict
    :ivar filename: str
    :ivar mode: str
    :ivar fileobj: io.IOBase
    :ivar compresslevel: int
    """
    
    def __init__(self, metadata:OrderedDict, filename=None, both_strands:bool=False, mode="w", fileobj=None, compresslevel=6):
        """Initilize the class."""
        if not isinstance(metadata, OrderedDict):
            raise TypeError("kmerdb.fileutil.KDBWriter expects a valid metadata object as its first positional argument")
        try:
            logger.debug("Validating metadata schema against the config.py header schema")
            jsonschema.validate(instance=dict(metadata), schema=config.graph_schema)
            self.metadata = metadata
            self.k = self.metadata['k']
        except jsonschema.ValidationError as e:
            logger.debug(metadata)
            logger.debug(e)
            logger.error("kmerdb.fileutil.KDBReader couldn't validate the header YAML")
            raise e

        if fileobj:
            assert filename is None, "kmerdb.graph expects filename to be None is fileobj handle is provided"
            handle = fileobj
        else:
            if "w" not in mode.lower() and "a" not in mode.lower():
                raise ValueError("Must use write or append mode, not %r" % mode)
            if "a" in mode.lower():
                raise NotImplementedError("Append mode is not implemented yet")
                # handle = _open(filename, "ab")
            else:
                handle = _open(filename, "wb")
        self._text = "b" not in mode.lower()
        self._handle = handle
        self._buffer = b"" if "b" in mode.lower() else ""
        self.compresslevel = compresslevel

        """
        Write the header to the file
        """

        
        logger.info("Constructing a new .kdbg file '{0}'...".format(self._handle.name))
        yaml.add_representer(OrderedDict, util.represent_ordereddict)
        if "b" in mode.lower():
            metadata_bytes = bytes(yaml.dump(self.metadata, sort_keys=False), 'utf-8')
            metadata_plus_delimiter_in_bytes = metadata_bytes + bytes(config.header_delimiter, 'utf-8')
            self.metadata["metadata_blocks"] = math.ceil( sys.getsizeof(metadata_plus_delimiter_in_bytes) / ( 2**16 ) ) # First estimate
            metadata_bytes = bytes(yaml.dump(self.metadata, sort_keys=False), 'utf-8')
            metadata_bytes = metadata_bytes + bytes(config.header_delimiter, 'utf-8')
            self.metadata["metadata_blocks"] = math.ceil( sys.getsizeof(metadata_bytes) / ( 2**16 ) ) # Second estimate
            metadata_bytes = bytes(yaml.dump(self.metadata, sort_keys=False), 'utf-8')
            metadata_bytes = metadata_bytes + bytes(config.header_delimiter, 'utf-8')
            logger.info("Writing the {0} metadata blocks to the new file".format(self.metadata["metadata_blocks"]))
            logger.debug(self.metadata)
            logger.debug("Header is being written as follows:\n{0}".format(yaml.dump(self.metadata, sort_keys=False)))
            # 03-04-2024 This is still not a completely functional method to write data to bgzf through the Bio.bgzf.BgzfWriter class included in BioPython
            # I've needed to implement a basic block_writer, maintaining compatibility with the Biopython bgzf submodule.
            #self.write(bytes(yaml.dump(metadata, sort_keys=False), 'utf-8'))
        
            for i in range(self.metadata["metadata_blocks"]):
                metadata_slice = metadata_bytes[:65536]
                metadata_bytes = metadata_bytes[65536:]
                self._write_block(metadata_slice)

                #self._write_block
                self._buffer = b""
                self._handle.flush()
        elif "w" == mode.lower() or "x" == mode.lower():
            self.write(yaml.dump(metadata, sort_keys=False))
            self._buffer = ""
            self._handle.flush()
        else:
            logger.error("Mode: {}".format(mode.lower()))
            raise RuntimeError("Could not determine proper encoding for write operations to .kdbg file")


    

def parsefile(filepath:str, k:int, quiet:bool=True, b:int=50000, both_strands:bool=False):
    """
    Parse a single sequence file.
    
    :param filepath: Path to a fasta or fastq file
    :type filepath: str
    :param k: Choice of k to shred k-mers with
    :type k: int
    :param b: Number of reads (per block) to process in parallel
    :type b: int
    :param both_strands: Strand specificity argument for k-mer shredding process
    :type both_strands: bool
    :returns: (graph, header_dictionary) header_dictionary is the file's metadata for the header block
    :rtype: (numpy.ndarray, dict)

    """

    from kmerdb import kmer
    if type(filepath) is not str:
        raise TypeError("kmerdb.graph.parsefile expects a str as its first positional argument")
    elif not os.path.exists(filepath):
        raise OSError("kmerdb.graph.parsefile could not find the file '{0}' on the filesystem".format(filepath))
    elif type(k) is not int:
        raise TypeError("kmerdb.graph.parsefile expects an int as its second positional argument")
    elif type(b) is not int:
        raise TypeError("kmerdb.parse.parsefile expects the keyword argument 'b' to be an int")
    elif type(both_strands) is not bool:
        raise TypeError("kmerdb.graph.parsefile expects the keyword argument 'both_strands' to be a bool")
    
    total_kmers = 4**k
    kmers = []

    
    logger.debug("Initializing graph parser...")
    seqprsr = seqparser.SeqParser(filepath, b, k)
    fasta = not seqprsr.fastq # Look inside the seqprsr object for the type of file
    
    # Initialize the graph structure
    Kmer = kmer.Kmers(k, strand_specific=not both_strands, verbose=fasta, all_metadata=True) # A wrapper class to shred k-mers with
    
    
    recs = [r for r in seqprsr]
    
    if not fasta:
        logger.debug("Read exactly {0} records from the seqparser object".format(len(recs)))
        assert len(recs) <= b, "The seqparser should return exactly {0} records at a time".format(b)
    else:
        logger.debug("Skipping the block size assertion for fasta files")
        
    logger.info("Read {0} sequences from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"))
        
        
    while len(recs):
        num_recs = len(recs)
        logger.info("Processing a block of {0} reads/sequences".format(num_recs))

        kmer_ids = [x for y in list(map(Kmer._shred_for_graph, recs)) for x in y]
        # Flatmap to 'kmer_ids', the dictionary of {'id': read_id, 'kmers': [ ... ]}
        logger.debug("\n\nAcquiring list of all k-mer ids from {0} sequence records...\n\n".format(num_recs))

        kmers = kmers + kmer_ids
        # END WHILE
        recs = [r for r in seqprsr] # The next block of exactly 'b' reads
        # This will be logged redundantly with the sys.stderr.write method calls at line 141 and 166 of seqparser.py (in the _next_fasta() and _next_fastq() methods)
        #sys.stderr("\n")
        logger.info("Read {0} more records from the {1} seqparser object".format(len(recs), "fasta" if fasta else "fastq"))

    logger.info("Read {0} k-mers from the input file".format(len(kmers)))

    sys.stderr.write("\n\n\n")
    logger.info("K-mer counts read from input(s)")
    sys.stderr.write("="*40 + "\n")
    logger.info("K-space dimension: {0}".format(4**k))
    logger.info("Number of k-mers: {0}".format(len(kmers)))
    sys.stderr.write("\n\n\n")

    logger.info("Constructing weighted edge list...")

    edge_list = make_graph(kmers, k, quiet=quiet)
    logger.info("Number of edges: {0}".format(len(edge_list)))
    
    sys.stderr.write("Edge list generation completed...")
    return (edge_list, seqprsr.header_dict())

def make_graph(kmer_ids, k:int=None, quiet:bool=True):
    """
    Make a neighbor graph from a k-mer id list.

    More specifically, this is a edge list

    This particular adjacency list will have redundant relationships:

    """
    if k is None:
        raise TypeError("kmerdb.graph.make_graph expects the keyword argument k to be an int")

    """
    Iterate over all k-mer ids, creating lists of neighbors derived from the kmer.
    Omit the self relationship, as the kmer is not a neighbor of itself.
    Convert to k-mer ids and add the adjacency tuple to the output

    A minimal neighbor graph is 
    """
    edges = {}
    for k_id in set(kmer_ids):
        current_kmer = kmer.id_to_kmer(k_id, k)

        char_goes_on_right = current_kmer[1:] # 11 characters (k - 1) (default 12), remove the left most, so new char goes on right to make the neighbor
        char_goes_on_left = current_kmer[:-1] # 11 characters ( k - 1)
        #common = current_kmer[1:-1] #  10 characters in common with all k-mers
        
        # Create the neighbor lists from the list of chars, prepended or appended to the k-mer suffix/prefix
        char_first = list(map(lambda c: c + char_goes_on_left, kmer.binaryToLetter))
        char_last = list(map(lambda c: char_goes_on_right + c, kmer.binaryToLetter))
        char_first_ids = list(map(kmer.kmer_to_id, char_first))
        char_last_ids = list(map(kmer.kmer_to_id, char_last))
        #r = {"kmer_id": k_id, "left_neighbors": lneighbors_ids, "right_neighbors": rneighbors_ids}
        neighbors_current_to_last = list(map(lambda x: (k_id, x), char_last_ids))
        neighbors_first_to_current = list(map(lambda x: (x, k_id), char_first_ids))
        neighbors = (neighbors_current_to_last + neighbors_first_to_current)

        """
        Hur dur

        at this point, I realized that that I had been removing the character on the wrong side from the new char insertion.
        """
        # if k_id == 15633431 or k_id == 12202077:
        #     logger.debug("-"*11)
        #     logger.debug("K-mer {0} and all (non-self) neighbors:".format(k_id))
        #     logger.debug(current_kmer)
        #     logger.debug("-"*11)
        #     logger.debug("<:" + ",".join(char_first) + " |||| " + common + "  ||||  " + ",".join(char_last) + ">:")
        #     logger.debug(",".join(list(map(str, char_first_ids))) + " |||| " + ",".join(list(map(str, char_last_ids))))
        for pair_ids in neighbors:
            edges[pair_ids] = 0

        """
        Old win condition. I wanted to see that a neighbor pair that was originally not found:
        The 12-mers with the following ids, which had this subsequence in common 'GTGGATAACCT', was not found in the keys of the edge list dictionary.
        """
        # if (15633431, 12202077) in neighbors:
        #     print("Sequence:")
        #     print(current_kmer)
        #     print("SUCCESS")
        #     print(edges)
        #     sys.exit(1)
    """
    At this point, the adjacency structure is a simple list of tuples. Redundant tuples may exist in disparate parts of the graph.

    The structure of the neighbor_graph is [(from_node_id, to_node_id), ...]
    The 'minimal neighbor graph' is provided by the variable 'min_neighbor_graph' from the preceding code block
    """

    """
    Next, the redundant k-mer id pairs are tallied, so that the edge weights are calculated
    """
    num_kmer_ids = len(kmer_ids)
    for i, val in enumerate(kmer_ids):
        if i+1 == num_kmer_ids:
            break
        else:
            """
            The following constitutes a 'neighboring' pair/edge of two k-mers (nodes):
            the k-mers are sequential in the original k-mer id list (implicitly sorted by virtue of the orientation of the read in the .fasta/.fastq file) 
            this list of k-mers actual neighbors can be validated by inspecting debug level output below (-vv).
            """
            pair_ids = (val, kmer_ids[i+1]) # pair of 'neighboring' k-mer ids
            pair = (kmer.id_to_kmer(val, k), kmer.id_to_kmer(kmer_ids[i+1], k))

            k1, k2 = pair
            k1 = k1[1:]
            k2 = k2[:-1]
            try:
                assert pair[0][1:] == pair[1][:-1], "kmerdb.graph expects neighboring k-mers to have at least k-1 residues in common. Must be from a fastq file."
            except AssertionError as e:
                logger.error("must be fastq")
                
                logger.error("mhmm")
                raise e
            if quiet is False:
                logger.debug("Edge: {0} => {1}".format(pair[0], pair[1]))
            try:
                # Accumulate for each edge pair 
                edges[pair_ids] += 1
            except KeyError as e:
                logger.error("Invalid key(pair) used to access edge list")
                sys.stderr.write("PAIR: {0} => {1}".format(pair[0], pair[1]))
                sys.stderr.write("pair ids: {0}".format(pair_ids))
                sys.stderr.write("Valid keys:")
                sys.stderr.write(list(filter(lambda e: e[0] == 15633431, edges.keys())))
                raise e

    return edges
        

class Parseable:
    def __init__(self, arguments):
        self.arguments = arguments
            
        
    def parsefile(self, filename):
        """Wrapper function for graph.parsefile to keep arguments succinct for deployment through multiprocessing.Pool
            
        :param filename: the filepath of the fasta(.gz)/fastq(.gz) to process with parsefile -> seqparser
        :type filename: str
        :returns: (db, m, n)
        """
        return parsefile(filename, self.arguments.k, quiet=not self.arguments.edges, b=self.arguments.fastq_block_size, both_strands=self.arguments.both_strands)



