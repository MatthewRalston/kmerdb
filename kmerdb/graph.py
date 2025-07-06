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


def graph_convert(n1:np.ndarray, n2:np.ndarray, weights:np.ndarray, k:int=None):
    """
    :param n1: A numpy.ndarray of node 1 kmer ids from an edge list
    :type numpy.ndarray:
    :param n2: A numpy.ndarray of node 2 kmer ids from an edge list
    :type numpy.ndarray:
    :param weights: A numpy.ndarray of weights corresponding to the edge
    :type numpy.ndarray:
    """

    if type(k) is not int:
        raise TypeError("kmerdb.graphlib.graph_convert() expects the keyword argument 'k' to be an int")
    elif isinstance(n1, np.ndarray) is False:
        raise TypeError("kmerdb.graphlib.graph_convert() expects a numpy.ndarray as its first positional argument")
    elif isinstance(n2, np.ndarray) is False:
        raise TypeError("kmerdb.graphlib.graph_convert() expects a numpy.ndarray as its second positional argument")
    elif isinstance(weights, np.ndarray) is False:
        raise TypeError("kmerdb.graphlib.graph_convert() expects a numpy.ndarray as its third positional argument")
    elif n1.shape == () or n2.shape == () or weights.shape == ():
        raise ValueError("kmerdb.graphlib.graph_convert() found that numpy.ndarrays were empty")
    elif n1.shape != n2.shape:
        raise ValueError("kmerdb.graphlib.graph_convert() found that the node1 and node2 array lengths did not match |  node1: {0} node2: {1}".format(n1.shape, n2.shape))
    elif n2.shape != weights.shape:
        raise ValueError("kmerdb.graphlib.graph_convert() found that the node2 and weights array lengths did not match |  node2: {0} weights: {1}".format(n2.shape, weights.shape))

    N = 4**k
    # edge_list = list(zip(n1.tolist(), n2.tolist(), weights.tolist()))

    
    edge_list_length = n1.shape[0] #
    edge_list = [ (int(n1[i]), int(n2[i]), {'weight': int(weights[i])}) for i in range(edge_list_length)]
    G = nx.Graph()


    G.add_nodes_from(range(N)) # Adds a list of k-mer ids as nodes to the graph
    G.add_edges_from(edge_list)

    print("Number of nodes: {0}".format(G.number_of_nodes()))
    print("Number of edges: {0}".format(G.number_of_edges()))


    
def make_coo_graph(filename:str, k:int, quiet:bool=True, replace_with_none:bool=False):
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
    seq_lengths = []
    pairs = []
    for seq in parse.parse_sequence_file(filename, return_tuple=False):
        """
        Shred sequences and get a list of kmer ids in order. We do not count kmers at the moment.
        Each k-mer needs to be converted to pairs of the id and all neighbors, 8 records total, corresponding to the neighboring kmer_ids by taking a single base from the left/right side of the kmer_id of interest.

        The maximum number (maximum k-mer id) that can be stored in a uint32 numpy.ndarray is 429967295 exactly
        therefore unit32 can only be used when k <= 16
        Store the results in numpy uint32 (~4B is the max size of unit32) if k<16 EXACTLY or uint64 otherwise
        FIXME: bonus_kmer_ids not included
        """
        seq_lengths.append(len(seq.seq))
        kmer_ids_, bonus_kmer_ids = kmer.shred(seq, k, replace_with_none=replace_with_none, quiet_iupac_warning=True)
        kmer_ids.extend(kmer_ids_)
        
        for kmer_id in kmer_ids:
            if kmer_id is not None:
                counts[kmer_id] += 1
        if replace_with_none is False:
            for kmer_id in bonus_kmer_ids:
                if kmer_id is not None:
                    counts[kmer_id] += 1
        pass
    num_kmer_ids = len(kmer_ids)
    pairs_ = []
    for i, k_id in enumerate(kmer_ids):
        current_kmer = kmer.id_to_kmer(k_id, k)
        #common_seq = current_kmer[1:-1]

        if i == 0:
            left_neighbor_ids = kmer.left_neighbors(kmer.id_to_kmer(k_id, k))
            pairs += list(zip( left_neighbor_ids, [k_id]*4))
        elif i == num_kmer_ids - 1:
            right_neighbor_ids = kmer.right_neighbors(kmer.id_to_kmer(k_id, k))
            pairs += list(zip( [k_id]*4, right_neighbor_ids ))
        else: # ALL OTHER KMERS
            left_neighbor = kmer.id_to_kmer(kmer_ids[i-1], k)
            right_neighbor = kmer.id_to_kmer(kmer_ids[i+1], k)
            has_right_neighbor = kmer.is_right_neighbor(current_kmer, right_neighbor)
            has_left_neighbor = kmer.is_left_neighbor(current_kmer, left_neighbor)
            # FIXME: Convert these to kmer_ids instead of appending str pairs to pairs_
            if has_left_neighbor is True and has_right_neighbor is True:
                pairs.append( (kmer_ids[i-1], k_id) ) # left neighbor edge
                pairs_.append( (left_neighbor, current_kmer ) )
            elif has_right_neighbor is True and has_left_neighbor is False:
                pairs.append( (k_id, kmer_ids[i+1]) )
            elif has_right_neighbor is False and has_left_neighbor is True:
                pairs.append( ( kmer_ids[i-1], k_id ) )
            elif has_right_neighbor is False and has_left_neighbor is False:
                raise ValueError("Internal error: Single k-mer found without obvious neighbors.")

        """
        FIXME: This produces all possible edges
        The commented code will not cover read boundaries...
        I actually DO want to not cover read boundaries, and then try to create as many Eulerian walks
        Through a graph of the EDGES as NODES and the edges as valid adjacencies

        MAKE A kmer.are_neighbors function that creates all neighbors of two k-mers and verifies that they are neighbors.

        This is how we take an EDGE as a node, and then take the left k-mer and find its possible neighbors, and the right node and its possible neighbors, and create an EDGE between the right k-mer and its possible neighbors so long as there are neighbors left? I dunno..
        """

        """
        [ (kmer_id, neighbor_kmer_id) , ... ]
        """

    md5, sha256 = util.checksum(filename)

    is_nullomer = np.where(counts == 0)
    nullomer_array = np.array(range(N), dtype="uint64")[is_nullomer]
    unique_kmers = int(np.count_nonzero(counts))

    
    num_nullomers = N - unique_kmers
    num_reads = len(seq_lengths)
    max_read_length = max(seq_lengths)
    min_read_length = min(seq_lengths)
    avg_read_length = int(np.mean(np.array(seq_lengths)))

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

    return pairs, file_metadata, counts



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
        :param n1_dtype: NumPy dtype
        :type n1_dtype: str
        :param n2_dtype: NumPy dtype
        :type n2_dtype: str
        :param weights_dtype: NumPy dtype
        :type weights_dtype: str
        :param sort: *unimplemented* - uh.. sort the .kdbg edges somehow
        :type sort: bool
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
        print("Line: {0}".format(line))
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
        total_kmer_ids = []
        kmer1 = []
        kmer2 = []
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
                    tkid, k1, k2 = list(map(int, line.split("\t")))
                    total_kmer_ids.append(tkid)
                    kmer1.append(k1)
                    kmer2.append(k2)
            except StopIteration as e:
                reading = False
        self.total_kmer_ids   = np.array(total_kmer_ids, dtype="uint64")
        self.kmer1_ids        = np.array(kmer1, dtype="uint64")
        self.kmer2_ids        = np.array(kmer2, dtype="uint64")

        return
    
