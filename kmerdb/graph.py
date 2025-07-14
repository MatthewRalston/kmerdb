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
import io
import re
import logging
from builtins import open as _open
from collections import deque, OrderedDict
import psutil
import json
import math

import yaml

import jsonschema
import numpy as np
import networkx as nx


from Bio import bgzf
from kmerdb import fileutil, parse, kmer, config, util


logger = logging.getLogger(__file__)



def make_secondary_graph_from_edge_list(seq_id:str, k:int, edges:list):
    """
    Creates a list of k+1 mer edges from a list of 2-tuples of kmerids (the primary edge list)
    :param seq_id: A sequence id for all corresponding edges
    :type str:
    :param edges: A list of 2-tuples of k-mer ids
    :type list:
    :raises TypeError: if seq_id is not a str
    :raises TypeError: if k is not an int
    :raises TypeError: if edges is not a list of N 2-tuples of kmer id int
    :raises ValueError: if the number of new edges does not have length of N-1
    """
    if type(seq_id) is not str:
        raise TypeError("kmerdb.graph.make_secondary_graph_from_edge_list() expects a str as its first positional argument")
    elif type(k) is not int:
        raise TypeError("kmerdb.graph.make_secondary_graph_from_edge_list() expects an int as its second positional argument")
    elif type(edges) is not list or not all(type(e) is tuple and type(e[0]) is int and type(e[1]) is int for e in edges):
        raise TypeError("kmerdb.graph.make_secondary_graph_from_edge_list() expects a list of integer pairs as its third positional argument")

    new_edges = []
    kplus1mers = convert_edge_list_to_augmented_kmer(k, edges)
    num_edges = len(kplus1mers)
    for i, kmer1 in enumerate(kplus1mers):
        if i < num_edges - 1:
            kmer2 = augmented_kplus1mers[i+1]
            new_edges.append( (kmer.kmer_to_id(kmer1), kmer.kmer_to_id(kmer2)) )
    if len(new_edges) != num_edges - 1:
        raise ValueError("Internal error: kmerdb.graph.make_secondary_graph_from_edge_list() should create N-1 k+1mers from an edges list with length N")
    return new_edges

def convert_edge_list_to_augmented_kmer(k:int, edges:list):
    """
    Creates a list of k+1mer strings from an edge list of 2-tuples of kmer ids
    :param k:
    :type k:
    :param edges:
    :type list:
    :raises TypeError: if k is not an int
    :raises TypeError: if edges is not a list of 2-tuples of kmer id int
    """
    if type(k) is not int:
        raise TypeError("kmerdb.graph.convert_edge_list_to_augmented_kmer() expects an int as its first positional argument")
    elif type(edges) is not list or not all(type(e) is tuple and type(e[0]) is int and type(e[1]) is int for e in edges):
        raise TypeError("kmerdb.graph.convert_edge_list_to_augmented_kmer() expects a list of integer pairs as its third positional argument")

    augmented_kplus1mers = []
    num_edges = len(edges)
    for i, edge in edges:
        kmerid1, kmerid2 = edge
        kmer1 = kmer.id_to_kmer(kmerid1, k)
        kmer2 = kmer.id_to_kmer(kmerid2, k)
        if kmer.is_right_neighbor(kmer1, kmer2):
            bridging_kmer = kmer1 + kmer2[-1]
            augmented_kplus1mers.append(kmer.kmer_to_id(bridging_kmer))
        else:
            logger.error(f"Sequence {seq_id}  |   {kmer1} @ position {i}  -  {kmer2} @ position {i+1}")
            raise ValueError(f"Internal error: lexicographically ordered kmers don't appear to be valid neighbors in the sequence.")
    return augmented_kplus1mers


def make_edges_from_kmerids(kmer_ids:list, positions:list, seqlen:int, seq_id:str, k:int):
    """
    Make edges from a list of lexicographically ordered kmer ids, given their positions and seq_ids
    :param kmer_ids: A list of kmer ids as ints
    :type list:
    :param positions: A list of positions in the 'seq_id' sequence/read
    :type list:
    :param seqlen: The length of the sequence
    :type int:
    :param seq_id: A sequence identifier from the fasta/fastq sequence
    :type str:
    :param k:  The choice of k
    :type int:
    :raises TypeError: if the kmer_ids positional arg is not a list of integers
    :raises TypeError: if the positions positional arg is not a list of integers
    :raises ValueError: if the list of kmer_ids and positions are not equal in length
    :raises TypeError: if the sequence identifier is not a str
    :raises TypeError: if the k param is not an int
    :returns: A list of 5-tuples - (sequence_id, position1, kmerid1, position2, kmerid2)
    :rtype: list
    """
    if type(kmer_ids) is not list and not all(type(kid) is int for kid in kmer_ids):
        raise TypeError("kmerdb.graph.make_edges_from_kmerids() expects a list of kmerids as its first positional argument")
    elif type(positions) is not list and not all(type(p) is int for p in positions):
        raise TypeError("kmerdb.graph.make_edges_from_kmerids() expects a list of int positions in the sequence as its second positional argument")
    elif len(kmer_ids) != len(positions):
        raise ValueError("kmerdb.graph.make_edges_from_kmerids() expects the lists of kmer_ids and positions to be equal in length")
    elif type(seqlen) is not int:
        raise TypeError("kmerdb.graph.make_edges_from_kmerids() expects an int sequence length as its third positional argument")
    elif type(seq_id) is not str:
        raise TypeError("kmerdb.graph.make_edges_from_kmerids() expects a sequence identifier as a str as its fourth positional argument")
    elif type(k) is not int:
        raise TypeError("kmerdb.graph.make_edges_from_kmerids() expects an int as its fifth positional argument")

    num_kmer_ids = len(kmer_ids)
    final = []
    for i in range(num_kmer_ids):
        kmer_id = kmer_ids[i]
        if kmer_id is not None:
            current_kmer = kmer.id_to_kmer(kmer_id, k)
            pos = positions[i]

            if positions[i] == 0:
                # If i = 0, the very first in the list, do nothing.
                continue
            elif pos + k + 1 == seqlen:
                # If the position is at the end of a sequence, append the last k-mer pair
                final.append( (seq_id, positions[i-1], kmer_ids[i-1], pos, kmer_id, ) ) # The final result is an edge tuple (sequence_id, position1, kmerid1, position2, kmerid2)
            else:
                """
                Here, internal (not terminal) k-mer adjacencies/edges are produced from a list of kmer_ids
                Local neighborhoods are created to correctly resolve k-mer adjacencies via lexical ordering while
                accounting for the correct kmer-to-kmer adjacency wrt ambiguous residues
                """
                left_neighbors = [] # Create local neighbor lists in case of substitution of ambiguous residue
                j = i - 1
                tracking_left = True
                while tracking_left is True:
                    if j >= 0 and j + 1 <= seqlen:
                        if positions[j] == pos: # The i-xth position has the same position in the sequence and is not a neighbor
                            pass
                        elif positions[j] == pos - 1 and kmer.is_left_neighbor(current_kmer, kmer.id_to_kmer(kmer_ids[j], k)): # This is considered a left neighbor lexically
                            left_neighbors.append( (kmer_ids[j], kmer_id) )
                        else: # This is an exit condition
                            tracking_left = False
                    else:
                        tracking_left = False
                    j -= 1

                """
                No longer need to check for right neighbors lexically,
                as all kmer edges are kmer_ids[i-1] -> kmer_id
                """
                # right_neighbors = []                    
                # j = i + 1
                # tracking_right = True
                # while tracking_right is True:
                #     if j < num_kmer_ids and seq_ids[j] == seq_id:
                #         #print(f"{i} : {current_kmer} -> {kmer.id_to_kmer(kmer_ids[j], k)}")
                #         if positions[j] == pos: # The i+xth position has the same position in the sequence and is not a neighbor
                #             pass
                #         elif positions[j] == pos + 1 and kmer.is_right_neighbor(current_kmer, kmer.id_to_kmer(kmer_ids[j], k)):
                #             right_neighbors.append( (kmer_id, kmer_ids[j]) )
                #         else:
                #             tracking_right = False
                #     else:
                #         tracking_right = False
                #     j += 1
                # #This may be needed for debugging purposes
                #ln = [tuple(map(lambda l: kmer.id_to_kmer(l, k), ln)) for ln in left_neighbors]
                #rn = [tuple(map(lambda l: kmer.id_to_kmer(l, k), rn)) for rn in right_neighbors]
                #print(f"{ln} | {i}:{current_kmer} ")
                pairs_ = resolve_edges(i, kmer_id, current_kmer, seq_id, pos, seqlen, left_neighbors, k)
                if len(pairs_) > 0:
                    final += [ (seq_id, positions[i-1], p[0], pos, kmer_id) for p in pairs_] # This
    return final

    
def make_edges_from_fasta(filename:str, k:int, quiet:bool=True, replace_with_none:bool=False):
    if type(filename) is not str:
        raise TypeError("kmerdb.graph.make_coo_graph() expects a fasta/fastq sequence filepath as a str as its first positional_regument")
    elif type(k) is not int:
        raise TypeError("kmerdb.graph.make_coo_graph() expects an int for k as the second positional argument")
    elif type(quiet) is not bool:
        raise TypeError("kmerdb.graph.make_coo_graph() expects the keyword argument 'slurp' to be a bool")

    if os.path.exists(filename) is False or os.access(filename, os.R_OK) is False:
        raise ValueError("kmerdb.graph.make_coo_graph() expects the filepath to be be readable on the filesystem")
    N = 4**k
    dtype = "uint32" if k < 17 else "uint64"
    """
    # {'i': tkmer_id, 'kmer_id': kmer_id, 'seq_id': seq_id, 'pos': pos}
    # (i, kmer_id, seq_id, pos)
    # We need a list of tuples with the following attributes:
    #      - the index of the kmer in the TOTAL k-mer sequence space
    #      - the kmer_id for the k-mer starting that position
    #      - the seq_id for the sequence that was shredded
    #      - and an ALMOST identical index of the position of the k-mer id
    """
    counts = np.zeros(N, dtype="uint64")
    kmer_ids = []
    seq_ids = []
    positions = []
    seq_lengths = {}
    data = []
    num_kmer_ids = 0
    for seq in parse.parse_sequence_file(filename, return_tuple=False):
        """
        Shred sequences and get a list of kmer ids in order. We do not count kmers at the moment.
        Each k-mer needs to be converted to pairs of the id and all neighbors, 8 records total, corresponding to the neighboring kmer_ids by taking a single base from the left/right side of the kmer_id of interest.

        The maximum number (maximum k-mer id) that can be stored in a uint32 numpy.ndarray is 429967295 exactly
        therefore unit32 can only be used when k <= 16
        Store the results in numpy uint32 (~4B is the max size of unit32) if k<16 EXACTLY or uint64 otherwise
        FIXME: bonus_kmer_ids not included
        """
        seqlen = len(seq.seq)
        seq_lengths[seq.id] = seqlen 
        kmer_ids_, seq_ids_, pos = kmer.shred(seq, k, replace_with_none=replace_with_none, quiet_iupac_warning=True)
        data_ = make_edges_from_kmerids(kmer_ids_, pos, seqlen, seq.id, k)
        num_kmer_ids += len(kmer_ids_)
        if None in kmer_ids:
            raise RuntimeError("Internal error: An invalid k-mer id was found within sequence '{0}'".format(seq.id))
        # kmer_ids += kmer_ids_
        # seq_ids += seq_ids_
        # positions += pos
        for kmer_id in kmer_ids:
            if kmer_id is not None:
                counts[kmer_id] += 1
        data += data

    """
    Break off functionality here. Need to pass a slice of k-mers and then resolve the k-mer group using conditionals.
    The while-loop functionality should not be included here.
    Instead, lists of k-mers passed to a function should resolve which neighbors get to be involved.
    Completed 7/13/25
    """
    """
    FIXME: This produces all possible edges
    The commented code will not cover read boundaries...
    I actually DO want to not cover read boundaries, and then try to create as many Eulerian walks
    Through a graph of the EDGES as NODES and the edges as valid adjacencies

    MAKE A kmer.are_neighbors function that creates all neighbors of two k-mers and verifies that they are neighbors.

    This is how we take an EDGE as a node, and then take the left k-mer and find its possible neighbors, and the right node and its possible neighbors, and create an EDGE between the right k-mer and its possible neighbors so long as there are neighbors left? I dunno..
    """

    md5, sha256 = util.checksum(filename)

    is_nullomer = np.where(counts == 0)
    nullomer_array = np.array(range(N), dtype="uint64")[is_nullomer]
    unique_kmers = int(np.count_nonzero(counts))

    
    num_nullomers = N - unique_kmers
    num_reads = len(seq_lengths)
    max_read_length = max(seq_lengths)
    min_read_length = min(seq_lengths)
    avg_read_length = int(np.mean(np.array(list(seq_lengths.values()))))

    file_metadata = {
        "filename": filename,
        "md5": md5,
        "sha256": sha256,
        "total_reads": len(seq_lengths),
        "total_kmers": num_kmer_ids,
        "unique_kmers": unique_kmers, # or N - len(nullomers),
        "nullomers": num_nullomers,
        "num_reads": num_reads, # FIXME
        "min_read_length": min_read_length,
        "max_read_length": max_read_length,
        "avg_read_length": avg_read_length
    }

    return data, file_metadata, counts




def resolve_edges(i:int, kmer_id:int, current_kmer:str, seq_id:str, pos:int, seqlen:int, left_neighbors:list, k:int): #kmer_ids:list, seq_ids:list, positions, seq_lengths:dict, left_neighbors:list, right_neighbors:list, k):
    if pos + 1 > seqlen: # when this is the last position in the sequence (seq_len - 1), validate neighbor and return final kmer_edge
        logger.error(f"kmer {i} '{current_kmer}' has {len(left_neighbors)} left neighbors and {len(right_neighbors)} right neighbors")
        raise IndexError(f"Internal error: {i}th kmer {kmer_id}:{current_kmer} has {len(left_neighbors)} left neighbors: {left_neighbors} at position {pos} in a sequence of length {seqlen}")
    # This sequence resolves
    generic_error = ValueError("kmerdb.graph.make_edges_from_fasta() - difficulty deciphering sequence boundary and which neighboring kmers are considered valid neighbors")
    
    if len(left_neighbors) == 1 and kmer.is_left_neighbor(current_kmer, kmer.id_to_kmer(left_neighbors[0][0], k)):
        assert left_neighbors[0] == (left_neighbors[0][0], kmer_id), "kmerdb.graph.resolve_edges() expects a pair of kmerids as the argument 'left_neighbors'"
        return [ left_neighbors[0] ] # The left neighbor is valid, do retrospective edge list append
    elif len(left_neighbors) > 1:
        # This is the condition in which the current k-mer has multiple 'possible' left/right neighbors
        return [ (ln[0], kmer_id) for ln in left_neighbors if kmer.is_left_neighbor(current_kmer, kmer.id_to_kmer(ln[0], k)) and ln[1] == kmer_id ]

    else:
        left_kmers = list(list(map(lambda km: (kmer.id_to_kmer(km[0], k), kmer.id_to_kmer(km[1], k)), left_neighbors)))
        logger.error(f"index: {i}th kmer @ {pos} in {seq_id} : \n\n{left_kmers} -> {current_kmer}\n\n")
        raise ValueError(f"Internal error: could not determine the category of k-mer relationships amongst the {i}th total k-mer '{current_kmer}' and its neighborhood amongst sequence '{seq_id}'")




def open(filename:str, metadata=None, mode="r", slurp:bool=False):
    if type(filename) is not str:
        raise TypeError("kmerdb.graph.open() expects a str as its first positional argument")
    elif type(mode) is not str or mode not in ("r", "w", "wb"):
        raise ValueError("kmerdb.graph.open() expects either 'r' or 'w' for file open 'mode'")
    elif type(slurp) is not bool:
        raise TypeError("kmerdb.graph.open() expects the keyword argument 'slurp' to be a bool")

    if mode == "r":
        return KDBGReader(filename, mode="r", slurp=slurp)
    elif "w" in mode:
        if metadata is None:
            raise ValueError("kmerdb.graph.open() requires the keyword argument 'metadata' to be specified")
            
        return KDBGWriter(filename, mode="w", metadata=metadata)


class KDBGWriter(bgzf.BgzfWriter):
    """
    A wrapper class around Bio.bgzf.BgzfWriter to write a .kdbg file to disk.


    :ivar filename: valid filepath
    :ivar mode: read/write mode
    :ivar metadata: header information
    :ivar both_strands: *unimplemented*
    :ivar fileobj: io.IOBase
    :ivar compresslevel: compression parameter
    """
    
    def __init__(self, filename:str=None, mode:str="w", metadata:OrderedDict=None, fileobj:io.IOBase=None, compresslevel:int=6):
        """
        A wrapper around Bio.bgzf.BgzfWriter

        :param filename: A valid filepath
        :type filename: str
        :param mode: read/write mode
        :type mode: str
        :param metadata: A metadata header for the .kdbg file
        :type metadata: collections.OrderedDict
        :param compresslevel: compression parameter
        :type compresslevel: int
        :raise TypeError: filename was not a str
        :raise TypeError: mode was invalid
        :raise TypeError: fileobj was invalid
        :raise TypeError: compresslevel was invalid
        :raise ValueError: mode was invalid
        :raise NotImplementedError: append mode was invalid
        """
        
        if filename is None or type(filename) is not str:
            raise TypeError("kmerdb.graph.KDBGWriter expects the filename to be a str")
        elif mode is None or type(mode) is not str:
            raise TypeError("kmerdb.graph.KDBGWriter expects the mode to be a str")
        elif metadata is None or (type(metadata) is not OrderedDict and type(metadata) is not dict):
            raise TypeError("kmerdb.graph.KDBGWriter - invalid metadata argument")
        elif fileobj is not None and not isinstance(fileobj, io.IOBase):
            raise TypeError("kmerdb.graph.KDBGWriter expects the keyword argument 'fileobj' to be an instance of io.IOBase")
        elif compresslevel is None or type(compresslevel) is not int:
            raise TypeError("kmerdb.graph.KDBGWriter expects the keyword argument 'compresslevel' to be an int")

        if fileobj:
            assert filename is None, "kmerdb.graph expects filename to be None is fileobj handle is provided"
            handle = fileobj
        else:
            if "w" not in mode.lower() and "a" not in mode.lower():
                raise ValueError("Must use write or append mode, not %r" % mode)
            elif "wb" == mode:
                pass
            elif mode == "w":
                pass
            elif "a" in mode.lower():
                raise NotImplementedError("Append mode is not implemented yet")
                # handle = _open(filename, "ab")
            else:
                raise RuntimeError("Unknown mode for .kdbg file writing class kmerdb.graph.KDBGWriter: {0}".format(mode))
        self._text = "b" not in mode.lower()
        self._handle = _open(filename, "wb")
        self._buffer = b"" if "b" in mode.lower() else ""
        self.compresslevel = compresslevel

        """
        Write the header to the file
        """

        logger.info("Constructing a new .kdbg file '{0}'...".format(self._handle.name))

        # 3-04-2024
        yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)
        self.metadata = metadata

        #self._write_block(metadata_slice)
        metadata_bytes = bytes(yaml.dump(self.metadata, sort_keys=False), 'utf-8')
        metadata_plus_delimiter_in_bytes = metadata_bytes + bytes(config.header_delimiter, 'utf-8')
        self.metadata["metadata_blocks"] = math.ceil( sys.getsizeof(metadata_plus_delimiter_in_bytes) / ( 2**16 ) ) # First estimate
        metadata_bytes = bytes(yaml.dump(self.metadata, sort_keys=False), 'utf-8')
        metadata_bytes = metadata_bytes + bytes(config.header_delimiter, 'utf-8')
        self.metadata["metadata_blocks"] = math.ceil( sys.getsizeof(metadata_bytes) / ( 2**16 ) ) # Second estimate
        metadata_bytes = bytes(yaml.dump(self.metadata, sort_keys=False), 'utf-8')
        metadata_bytes = metadata_bytes + bytes(config.header_delimiter, 'utf-8')

        logger.info("Writing the {0} metadata blocks to the new file".format(self.metadata["metadata_blocks"]))
        logger.info("Header is being written as follows:\n{0}".format(yaml.dump(self.metadata, sort_keys=False)))

            # 01-01-2022 This is still not a completely functional method to write data to bgzf through the Bio.bgzf.BgzfWriter class included in BioPython
            # I've needed to implement a basic block_writer, maintaining compatibility with the Biopython bgzf submodule.
            #self.write(bytes(yaml.dump(metadata, sort_keys=False), 'utf-8'))
        
        for i in range(self.metadata["metadata_blocks"]):
            metadata_slice = metadata_bytes[:65536]
            metadata_bytes = metadata_bytes[65536:]
            self._write_block(metadata_slice)

            #self._write_block
            self._buffer = b""
            self._handle.flush()
        


class KDBGReader(bgzf.BgzfReader):
    """
    A class to read .kdbg files.
    1 args

    8 keyword args

    fileobj
    mode
    n1_dtype
    n2_dtype
    weights_dtype
    sort
    slurp

    :ivar filename: The choice of k to shred with
    :ivar fileobj: An existing fileobject from io.IOBase
    :ivar mode: read/write mode
    :ivar n1_dtype: NumPy dtype
    :ivar n2_dtype: NumPy dtype
    :ivar weights_dtype: NumPy dtype
    :ivar sort: *unimplemented* - uh.. sort the .kdbg edges somehow
    :ivar slurp: Autoload the .kdbg file

    
    """
    def __init__(self, filename:str, fileobj:io.IOBase=None, mode:str="r", slurp:bool=False):
        """
        A wrapper around Bio.bgzf.BgzfReader

        :param filename: A valid filepath
        :type filename: str
        :param fileobj: An existing fileobject from io.IOBase
        :type fileobj: io.IOBase
        :param mode: read/write mode
        :type mode: str
        :param slurp: Autoload the .kdbg file
        :type slurp: bool
        :raise TypeError: fileobj was not a file-type interface, inheriting from io.IOBase
        :raise TypeError: filename was not a str
        :

        """
        """

         fileobj type checking
        """


        if fileobj is not None and not isinstance(fileobj, io.IOBase):
            raise TypeError("kmerdb.graph.KDBGReader expects the keyword argument 'fileobj' to be a file object")
        if filename is not None and type(filename) is not str:
            raise TypeError("kmerdb.graph.KDBGReader expects the keyword argument 'filename' to be a str")
        elif not os.access(filename, os.R_OK):
            raise ValueError("kmerdb.graph.open expects an existing filepath as its first positional argument")
        elif mode is not None and type(mode) is not str:
            raise TypeError("kmerdb.graph.KDBGReader expects the keyword argument 'mode' to be a str")
        elif slurp is not None and type(slurp) is not bool:
            raise TypeError("kmerdb.graph.KDBGReader expects the keyword argument 'slurp' to be a bool")

        """

        Handle fileobj or 
        """
        
        if fileobj:
            assert filename is None
            handle = fileobj
            assert "b" in handle.mode.lower()
        else:
            if "w" in mode.lower() or "a" in mode.lower():
                raise ValueError("Must use read mode (default), not write or append mode")
            handle = _open(filename, "rb")

        self._text = "b" not in mode.lower()
        if self._text:
            self._newline = "\n"
        else:
            self._newline = b"\n"
        self._handle       = handle
        self._filepath     = self._handle.name
        #       Placeholder
        self.max_cache     = 100
        self._buffers      = {}
        self._block_start_offset = None
        self._block_raw_length = None
        """

        np.array defs

        n1
        n2
        weights

        """
        self.read_and_validate_kdbg_header()
        if slurp is True:
            logger.info("Reading .kdbg contents into memory")
            self.slurp()

        if handle is not None:
            handle.close()
        if self._handle is not None:
            self._handle.close()
        self._handle = None
        handle = None
        fileobj=None
        return

    def read_and_validate_kdbg_header(self):

        """

        KDBGReader PyYAML                metadata ingestion

        """

        '''
        Here we want to load the metadata blocks. We want to load the first two lines of the file: the first line is the version, followed by the number of metadata blocks
        '''
        # 0th block

        logger.info("Loading the 0th block from '{0}'...".format(self._filepath))
        self._load_block(self._handle.tell())

        self._buffer = self._buffer.rstrip(config.header_delimiter)
        
        initial_header_data = OrderedDict(yaml.safe_load(self._buffer))
        # Placeholder
        num_header_blocks = None
        if type(initial_header_data) is str:
            logger.error("Inappropriate type for the header data.")
            #logger.info("Um, what the heckin' is this in my metadata block?")
            raise TypeError("kmerdb.graph.KDBGReader could not parse the YAML formatted metadata in the first blocks of the file")
        elif type(initial_header_data) is OrderedDict:
            logger.info("Successfully parsed the 0th block of the file, which is expected to be the first block of YAML formatted metadata")
            logger.info("Assuming YAML blocks until delimiter reached.")

        """

        KDBGReader
        
        YAML metadata spec validation
        """
            
        if "version" not in initial_header_data.keys():
            raise TypeError("kmerdb.graph.KDBGReader couldn't validate the header YAML")
        elif "metadata_blocks" not in initial_header_data.keys():
            raise TypeError("kmerdb.graph.KDBGReader couldn't validate the header YAML")
        logger.info(initial_header_data)
        if initial_header_data["metadata_blocks"] != 1:
            logger.error("Internal error: does not support more than 1 metadata block yet")
            raise IOError("Internal error: Cannot read more than 1 metadata block yet")
        for i in range(initial_header_data["metadata_blocks"] - 1):
            logger.warning("Multiple metadata blocks read, most likely from a composite edge-graph...")
            self._load_block(self._handle.tell())
            addtl_header_data = yaml.safe_load(self._buffer.rstrip(config.header_delimiter))
            if type(addtl_header_data) is str:
                logger.error(addtl_header_data)
                raise TypeError("kmerdb.graph.KDBGReader determined the data in the {0} block of the header data from '{1}' was not YAML formatted".format(i, self._filepath))
            elif type(addtl_header_data) is dict:
                sys.stderr.write("\r")
                sys.stderr.write("Successfully parsed {0} blocks of YAML formatted metadata".format(i))
                initial_header_data.update(addtl_header_data)
                num_header_blocks = i
            else:
                logger.error(addtl_header_data)
                raise RuntimeError("kmerdb.graph.KDBGReader encountered a addtl_header_data type that wasn't expected when parsing the {0} block from the .kdb file '{1}'.".format(i, self._filepath))

        #raise RuntimeError("kmerdb.graph.KDBGReader encountered an unexpected type for the header_dict read from the .kdb header blocks")
        logger.info("Validating the header YAML...")
        try:
            jsonschema.validate(instance=initial_header_data, schema=config.graph_schema)
            self.metadata = dict(initial_header_data)
            self.k = self.metadata['k']

        except jsonschema.ValidationError as e:
            logger.error("kmerdb.graph.KDBGReader couldn't validate the header/metadata YAML from {0} header blocks".format(num_header_blocks))
            logger.error("""


Failed to validate YAML header.

-------------------------------



 You can store unique k-mer counts, total nullomer counts, and other metadata alongside the weighted edge list with the 'kmerdb graph' command


 ...then try again to view the header with 'kmerdb header'



""")

                
            logger.error(e.__str__)
            raise ValueError("Requires kmerdb v{0} format YAML header. Body is .tsv format table, .bgzf file.       - weighted edge list        (idx, node1, node2, weight)")
        self.metadata["header_offset"] = self._handle.tell()
        #logger.debug("Handle set to {0} after reading header, saving as handle offset".format(self.metadata["header_offset"]))
        #self._reader = gzip.open(self._filepath, 'r')
        self._offsets = deque()
        for values in bgzf.BgzfBlocks(self._handle):
            self._offsets.appendleft(values) # raw start, raw length, data start, data length
        if len(self._offsets) == 0:
            raise IOError("kmerdb.graph.KDBGReader opened an empty file")
        # Skip the zeroth block
        self._load_block()

    def read_line(self):
        """
        Read and parse a single line from the .kdbg file
        :returns: node1_id, node2_id, weight
        :rtype: tuple
        """
        # FIXME
        line = self.readline()
        #print("Line: {0}".format(line))
        if line == '':
            return
        else:
            return list(map(int, line.split("\t")))
            #return parse_kdbg_table_line(line, row_dtype="uint64")

    def slurp(self):
        """
        Autoload the .kdbg file into memory
        """
        #vmem = psutil.virtual_memory()
        seq_ids = []
        pos1 = []
        pos2 = []
        kmerid1 = []
        kmerid2 = []

        reading=True
        while reading is True:
            try:
                line = next(self)
            except StopIteration as e:
                logger.info("Finished loading .kdbg through slurp (on init)")
                break
            if line is None:
                raise RuntimeError("Internal error: kmerdb.graph.KDBGReader.slurp errored on new line")
            try:
                if line is not None:
                    line = line.split("\t")
                    if len(line) != 8:
                        raise ValueError("kmerdb.graph.KDBGReader.slurp() expects .kdbg data to have 8 columns")
                    i, seq_id, p1, id1, kmer1str, p2, id2, kmer2str = line
                    try:
                        if re.match(util.is_integer, i) is None:
                            raise ValueError("kmerdb.graph.KDBGReader.slurp() expects the first column to be an integer value")
                        elif re.match(util.is_integer, p1) is None:
                            raise ValueError("kmerdb.graph.KDBGReader.slurp() expects the third column to be an integer value")
                        elif re.match(util.is_integer, id1) is None:
                            raise ValueError("kmerdb.graph.KDBGReader.slurp() expects the fourth column to be an integer value")
                        elif re.match(util.is_integer, p2) is None:
                            raise ValueError("kmerdb.graph.KDBGReader.slurp() expects the sixth column to be an integer value")
                        elif re.match(util.is_integer, id2) is None:
                            raise ValueError("kmerdb.graph.KDBGReader.slurp() expects the seventh column to be an integer value")
                        
                    except Exception as e:
                        raise e
                    seq_ids.append(seq_id)
                    pos1.append(int(p1))
                    pos2.append(int(p2))
                    kmerid1.append(int(id1))
                    kmerid2.append(int(id2))
            except StopIteration as e:
                reading = False
        self.seq_ids = seq_ids
        self.pos1 = pos1 #np.array(pos1, dtype="uint64")
        self.kmer_id1 = kmerid1 #np.arary(kmerid1, dtype="uint64")
        self.pos2 = pos2 #np.array(pos2, dtype="uint64")
        self.kmer_id2 = kmerid2 #np.arary(kmerid2, dtype="uint64")
        self.G = None

        return
    
    def as_networkx(self):
        """
        Creates a NetworkX representation of the graph.
        
        """

        if self.G is not None:
            return self.G
        else:
            N = 4**self.k
            num_edges = len(self.kmer_id1.shape[0])

            unique_kmers = list(set(self.kmer_id1))
            temp = self.kmer_id1 + self.kmer_id2

            nodes = []
            for kmerid_ in temp:
                try:
                    idx = self.kmer_id1.index(kmerid_)
                    pos = self.pos1[idx]
                except ValueError as e:
                    idx = self.kmer_id2.index(kmerid_)
                    pos = self.pos2[idx]
                nodes.append( (kmerid_, {"seq_id": self.seq_ids[idx], "pos": pos, "kmer": kmer.id_to_kmer(kmerid_, self.k)}) )
            #tuples = [(self.seq_ids[i], self.pos1[i], self.kmer_id1[i], self.pos2[i], self.kmer_id2[i]) for i in range(num_edges)]
        
        
            edge_list = list(zip(self.kmer_id1, self.kmer_id2))
            
            G = nx.Graph()


            G.add_nodes_from(nodes) # Adds a list of k-mer ids as nodes to the graph
            G.add_edges_from(edge_list) # Adds edges between nodes

            logger.info(f"Number of nodes: {G.number_of_nodes()}")
            logger.info(f"Number of edges: {G.number_of_edges()}")

            self.G = G

            # Adjacency matrix
            #rcm = list(cuthill_mckee_ordering(self.G, heuristic=_smallest_degree))
            #self.A = nx.adjacency_matrix(self.G, nodelist=rcm)
            return G

        
    #def make_eulerian_circuit(self):

    def _smallest_degree(self):
        return min(self.G, key=self.G.degree)
        
    def write_dot(self, filepath:str):
        from networkx.utils import cuthill_mckee_ordering
        if type(filepath) is not str:
            raise ValueError("kmerdb.graph.KDBGReader.write_dot() expects a filepath as a str")
        #pg = nx.nx_pydot.to_pydot(G)
        nx.write_dot(G, filepath)

    # def read_dot(self, filepath):
    #     if type(filepath) is not str and not (os.path.exists(filepath) and os.access(filepath, os.R_OK)):
    #         raise ValueError("kmerdb.graph.KDBGReader.read_dot() expects a filepath as a str that points to a readable filepath on the filesystem")
                 
    #def write_sparse6(self, filepath:str):


    
