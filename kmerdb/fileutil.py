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
import sys
import os
import re
import logging
import gzip
import tempfile
import json
from collections import deque, OrderedDict
import psutil
#import re
from builtins import open as _open

import yaml

import numpy as np
import jsonschema
from Bio import SeqIO, bgzf

from kmerdb import util, config, appmap

# Logging configuration
logger = logging.getLogger(__file__)



def open(filepath, mode="r", metadata=None, sort:bool=False, slurp:bool=False, with_index:bool=False):
    """
    Opens a file for reading or writing. Valid modes are 'xrwbt'. 'metadata=' is needed when writing/creating.
    Returns a lazy-loading KDBReader object or a KDBWriter object.
    The data may be force loaded with 'slurp=True'

    :param filepath:
    :type filepath: str
    :param mode:
    :type mode: str
    :param metadata: The file header/metadata dictionary to write to the file.
    :type metadata: dict
    :param sort: Sort on read the data into KDBReader
    :type sort: bool
    :param slurp: Immediately load all data into KDBReader
    :type slurp: bool
    :param with_index: Create index data dynamically.
    :type with_index: bool
    :returns: kmerdb.fileutil.KDBReader/kmerdb.fileutil.KDBWriter
    :rtype: kmerdb.fileutil.KDBReader
    """
    if type(filepath) is not str:
        raise TypeError("kmerdb.fileutil.open expects a str as its first positional argument")
    ## We don't want this because it won't write a new .kdb file on mode='w'
    # elif os.access(filepath, os.R_OK) is False:
    #     raise ValueError("kmerdb.fileutil.open expects a valid filepath as its first positional argument")
    elif type(mode) is not str:
        raise TypeError("kmerdb.fileutil.open expects the keyword argument 'mode' to be a str")

    if (mode == "w" or mode == "x") and (metadata is not None and (isinstance(metadata, OrderedDict) or type(metadata) is dict)):
        pass
    elif (mode == "w" or mode == "x"):
        raise TypeError("kmerdb.fileutil.open expects an additional metadata dictionary")
    elif sort is None or type(sort) is not bool:
        raise TypeError("kmerdb.fileutil.open expects a boolean for the keyword argument 'sort'")
    elif slurp is None or type(slurp) is not bool:
        raise TypeError("kmerdb.fileutil.open expects a boolean for the keyword argument 'slurp'")
    elif with_index is None or type(with_index) is not bool:
        raise TypeError("kmerdb.fileutil.open expects a boolean for the keyword argument 'with_index'")

    
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
        return KDBReader(filename=filepath, mode=mode, sort=sort, slurp=slurp, with_index=with_index)
    elif "w" in mode.lower() or "x" in mode.lower():
        return KDBWriter(metadata, filename=filepath, mode=mode)
    else:
        raise ValueError("Bad mode %r" % mode)





class KDBReader(bgzf.BgzfReader):
    """
    A class that reads .kdb files, potentially just for accessing header metadata, or for reading the entire contents into numpy arrays.
    

    :ivar filename: str
    :ivar fileobj: io.IOBase
    :ivar mode: str
    :ivar max_cache: int
    :ivar sort: bool
    :ivar slurp: bool
    :ivar with_index: bool
    """
    def __init__(self, filename:str=None, fileobj:io.IOBase=None, mode:str="r", max_cache:int=100, sort:bool=False, slurp:bool=False, with_index:bool=False):
        if fileobj is not None and not isinstance(fileobj, io.IOBase):
            raise TypeError("kmerdb.fileutil.KDBReader expects the keyword argument 'fileobj' to be a file object")
        elif filename is not None and type(filename) is not str:
            raise TypeError("kmerdb.fileutil.KDBReader expects the keyword argument 'filename' to be a str")
        elif os.access(filename, os.R_OK) is False:
            raise ValueError("kmerdb.fileutil.KDBReader expects a valid filepath as its first positional argument")

        elif sort is None or type(sort) is not bool:
            raise TypeError("kmerdb.fileutil.KDBReader expects the keyword argument 'sort' to be a bool")
        elif max_cache is None or type(max_cache) is not int:
            raise TypeError("kmerdb.fileutil.KDBReader expects the keyword argument 'max_cache' to be a int")
        elif mode is None or type(mode) is not str:
            raise TypeError("kmerdb.fileutil.KDBReader expects the keyword argument 'mode' to be a str")
        elif slurp is None or type(slurp) is not bool:
            raise TypeError("kmerdb.fileutil.KDBReader expects the keyword argument 'slurp' to be a bool")
        elif with_index is None or type(with_index) is not bool:
            raise TypeError("kmerdb.fileutil.KDBReader expects the keyword argument 'with_index' to be a bool")
        
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
        self._handle      = handle
        self._filepath    = self._handle.name
        self.max_cache    = max_cache
        self._buffers     = {}
        self._block_start_offset = None
        self._block_raw_length = None
        self.kmer_ids     = None
        self.counts       = None
        self.frequencies  = None
        # Has been slurped
        self.completed = False
        '''
        Here we want to load the metadata blocks. We want to load the first two lines of the file: the first line is the version, followed by the number of metadata blocks
        '''
        # 0th block
        logger.info("Loading the 0th block from '{0}'...".format(self._filepath))
        self._load_block(self._handle.tell())

        self._buffer = self._buffer.rstrip(config.header_delimiter)
        
        initial_header_data = OrderedDict(yaml.safe_load(self._buffer))
        num_header_blocks = None

        if type(initial_header_data) is str:
            logger.info("Inappropriate type for the header data.")
            #logger.info("Um, what the fuck is this in my metadata block?")
            raise TypeError("kmerdb.fileutil.KDBReader could not parse the YAML formatted metadata in the first blocks of the file")
        elif type(initial_header_data) is OrderedDict:
            logger.info("Successfully parsed the 0th block of the file, which is expected to be the first block of YAML formatted metadata")
            logger.info("Assuming YAML blocks until delimiter reached.")
            if "version" not in initial_header_data.keys():
                raise TypeError("kmerdb.fileutil.KDBReader couldn't validate the header YAML")
            elif "metadata_blocks" not in initial_header_data.keys():
                raise TypeError("kmerdb.fileutil.KDBReader couldn't validate the header YAML")
            else:
                logger.info(str(initial_header_data))
                if initial_header_data["metadata_blocks"] == 1:
                    logger.info("1 metadata block: Not loading any additional blocks")
                    pass
                else:
                    for i in range(initial_header_data["metadata_blocks"] - 1):
                        self._load_block(self._handle.tell())
                        logger.debug("Multiple metadata blocks read...")
                        logger.debug("Secondary blocks read")
                        addtl_header_data = yaml.safe_load(self._buffer.rstrip(config.header_delimiter))
                        if type(addtl_header_data) is str:
                            logger.error(addtl_header_data)
                            logger.error("Couldn't determine this block.::/")
                            raise TypeError("kmerdb.fileutil.KDBReader determined the data in the {0} block of the header data from '{1}' was not YAML formatted".format(i, self._filepath))
                        elif type(addtl_header_data) is dict:
                            sys.stderr.write("\r")
                            sys.stderr.write("Successfully parsed {0} blocks of YAML formatted metadata".format(i))
                            initial_header_data.update(addtl_header_data)
                            num_header_blocks = i
                        else:
                            logger.error(addtl_header_data)
                            raise RuntimeError("kmerdb.fileutil.KDBReader encountered a addtl_header_data type that wasn't expected when parsing the {0} block from the .kdb file '{1}'.".format(i, self._filepath))
        else:
            raise RuntimeError("kmerdb.fileutil.KDBReader encountered an unexpected type for the header_dict read from the .kdb header blocks")
        logger.info("Reading additional blocks as YAML...")
        logger.info("Validating the header data against the schema...")
        try:
            jsonschema.validate(instance=initial_header_data, schema=config.kdb_metadata_schema)
            self.metadata = dict(initial_header_data)



            # Metadata values
            self.k = self.metadata['k']
            self.sorted = self.metadata["sorted"]

            N = 4**self.metadata["k"]
            # Logging statements
            logger.debug(f"Allocating a {N} long uint64 NumPy array... probably /s")
            # kmer id, count, frequency vectors
            self.kmer_ids = np.zeros(N, dtype="uint64")
            self.counts = np.zeros(N, dtype="uint64")
            self.frequencies = np.zeros(N, dtype="float64")
            # .kdb index
            self.index = None
            self.with_index = with_index
            
        except jsonschema.ValidationError as e:
            logger.error("kmerdb.fileutil.KDBReader couldn't validate the header/metadata YAML from {0} header blocks".format(num_header_blocks))
            logger.error("""


Failed to validate YAML header.

-------------------------------



 You can store unique k-mer counts, total nullomer counts, and other metadata alongside a k-mer count vector with the 'kmerdb profile' command



{0}

 ...then try again to view the header with 'kmerdb header' or the whole file : 'kmerdb view -H kmer_count_vector.8.kdb'



""".format(appmap.command_1_usage))

            logger.error(e.__str__)
            raise ValueError("Requires kmerdb v{0} format YAML header. Body is .tsv format table, .bgzf file.      - k-mer count vector        (idx, k-mer id, count, frequency)".format(config.VERSION))
        self.metadata["header_offset"] = self._handle.tell()
        logger.debug("Handle set to {0} after reading header, saving as handle offset".format(self.metadata["header_offset"]))
        #self._reader = gzip.open(self._filepath, 'r')
        self._offsets = deque()
        for values in bgzf.BgzfBlocks(self._handle):
            #logger.debug("Raw start %i, raw length %i, data start %i, data length %i" % values)
            self._offsets.appendleft(values) # raw start, raw length, data start, data length
        if len(self._offsets) == 0:
            raise IOError("kmerdb.fileutil.KDBReader opened an empty file")
        # Skip the zeroth block
        self._load_block()
        # print(str(self._buffer)) # 1
        # print(self.readline())
        # self._load_block()
        # print(self._buffer) # 2
        # print(self.readline())

        if sort is True and self.metadata["sorted"] is True:
            sort is True
        else:
            sort is False
        if slurp is True:
            self.slurp(sort=sort, with_index=with_index)
        handle.close()
        self._handle.close()
        if slurp is True:
            self._handle = None
            handle = None
        fileobj=None
        return
    
    def __iter__(self):
        return self

    def __exit__(self, type, value, tb):
        if self._handle is not None:
            self._handle.close()
        return


    def _slurp(self, sort:bool=False, with_index:bool=False):
        if sort is None or type(sort) is not bool:
            raise TypeError("kmerdb.fileutil.KDBReader._slurp expects the 'sort' keyword argument to be a bool")
        elif with_index is None or type(with_index) is not bool:
            raise TypeError("kmerdb.fileutil.KDBReader._slurp expects the 'with_index' keyword argument to be a bool")
        if self._handle is None:
            self._handle = _open(self._filepath, 'rb')
        # First calculate the amount of memory required by the array
        N = 4**self.k # The dimension of the k-space, or the number of elements for the array
        num_bytes = 4 * N
        #logger.info("Approximately {0} bytes".format(num_bytes))
        vmem = psutil.virtual_memory()
        # Initial zero-valued empty id/count vectors
        self.kmer_ids = np.zeros(N, dtype="uint64")
        self.counts = np.zeros(N, dtype="uint64")
        self.frequencies = np.zeros(N, dtype="float64")
        self.all_kmer_metadata = []
        # Initial Python lists
        kmer_ids = []
        profile = []
        counts = []
        frequencies = []
        
        # Block offsets
        offsets = []
        if vmem.available > num_bytes:
            # Do the slurp
            logger.info("Reading profile into memory")
            if sort is False:
                i = 0
                try:
                    for j in range(N):
                        #logger.debug("Reading {0}th line...".format(j))
                        try:
                            line = next(self)
                            offset = self._handle.tell()
                        except StopIteration as e:
                            logger.error("Finished loading initial profile through slurp-on-init")
                            logger.error("Read profile from: {0}".format(self._filepath))
                            pass
                        if line is None:
                            logger.error("Internal error: 'next' was None")
                            raise RuntimeError("Internal error: Could not complete parsing of file.")
                        # Don't forget to not parse the metadata column [:-1]
                        l = line.rstrip().split("\t")
                        try:
                            assert len(l) == config.KDB_COLUMN_NUMBER, "kmerdb.fileutil.KDBReader requires .kdb files to have {0} columns.".format(config.KDB_COLUMN_NUMBER)
                        except AssertionError as e:
                            logger.error("Found .kdb table line with {0} columns. Requires {1} columns.".format(l, config.KDB_COLUMN_NUMBER))
                            logger.error(line)
                            raise e
                        x, kmer_id, _count, _frequency = l
                        x = int(x)
                        kmer_id = int(kmer_id)                            
                        assert j == x, "Line number doesn't match file row index"
                        kmer_ids.append(kmer_id)
                        count = int(_count)
                        frequency = float(count)/N
                        # K-mer id and count NumPy index-based assignment
                        self.kmer_ids[j] = kmer_id
                        self.counts[kmer_id] = count
                        self.frequencies[kmer_id] = frequency
                        # Block offsets
                        if with_index is True:
                            offsets.append((kmer_id, offset))
                        i += 1
                except StopIteration as e:
                    if i == N:
                        raise e
                    else:
                        logger.debug("Number of lines read: {0}".format(i))
                        logger.debug("Number of k-mers (N) set globally.".format(N))
                        logger.debug("Read only {0} lines from the file...".format(i))
                        logger.error(e.__str__())
                        raise e
                except Exception as e:
                    raise e
                assert self.counts.size == self.kmer_ids.size, "length of k-mer id column has diverged from the counts column"
                assert self.counts.size == self.frequencies.size, "length of frequencies column has diverged from counts column"
            elif sort is True:
                i = 0
                try:
                    for j in range(N):
                        try:
                            line = next(self)
                            offset = self._handle.tell()
                        except StopIteration as e:
                            logger.error("Finished loading initial profile through slurp-on-init")
                            raise e
                        if line is None:
                            logger.warning("Next was None... profile was sparse, breaking", "WARNING")
                            raise IOError("next method was None on sorted import. Internal Error.")
                        # Don't forget to not parse the metadata column [:-1]
                        _, kmer_id, _count, _frequency = line.rstrip().split("\t")
                        if re.match(util.is_integer, kmer_id) is None:
                            raise ValueError("kmerdb.fileutil.KDBReader._slurp() expects the second column to be an integer value")
                        elif re.match(util.is_integer, _count) is None:
                            raise ValueError("kmerdb.fileutil.KDBReader._slurp() expects the third column to be an integer value")
                        elif re.match(util.findall_float, _frequency) is None:
                            raise ValueError("kmerdb.fileutil.KDBReader._slrup() expects the fourth column to be a floating point value")
                        kmer_id = int(kmer_id)
                        count = int(_count)
                        _frequency = float(_frequency)
                        # Trying to get the needle threaded here.
                        # The goal is for the profile to be a running index, similar to the row number
                        # But this could maintain lexical sort order
                        self.kmer_ids[j] = kmer_id
                        self.counts[kmer_id] = count
                        self.frequencies[kmer_id] = _frequency
                        # Block offsets
                        if with_index is True:
                            offsets.append((kmer_id, offset))
                        # The index data, just a list of id and gzip file-offsets
                        i+=1
                        #sys.stderr.write("::DEBUG::   |\-)(||||..... KMER_ID: {0} COUNT: {1}".format(kmer_id, count))
                    assert self.kmer_ids.size == self.counts.size, "Counts size has diverged from kmer_ids size/shape"
                    assert self.counts.size == self.frequencies.size, "Frequencies size has diverged from counts shape"
                    # Index tuple  data
                    #offsets.append((kmer_id, offset))
                except StopIteration as e:
                    if i == N:
                        raise e
                    else:
                        raise e
                except Exception as e:
                    raise e
                logger.info("Read {0} lines from the file...".format(i))
                assert self.kmer_ids.size == self.counts.size, "Counts size is mismatched in with kmer_ids size"
                assert self.frequencies.size == self.counts.size, "Frequencies size is mismatched in count from counts size"
                logger.debug("Correct array sizes match...")
                self._handle.seek(0)
                self._load_block()
                #logger.debug("Dammit, why can't i reset the Bio.bgzf filehandle...")
                if sort is True:
                    # If the file is sorted, do not sort
                    kmer_ids_sorted_by_count = np.lexsort((self.kmer_ids, self.counts))
                    reverse_kmer_ids_sorted_by_count = np.flipud(kmer_ids_sorted_by_count)
                    for i, idx in enumerate(kmer_ids_sorted_by_count): # This is right, not fixing this.
                        kmer_id = self.kmer_ids[idx]
                        count = self.counts[idx]
                        ## I stand corrected. m
                        #self.kmer_ids[i] = kmer_ids[profile[i]]
                        #self.counts[idx] = counts[idx]
                        ###self.frequencies[idx] = frequencies[idx]
                        #logger.debug("Just in casey eggs and bakey...")
            else:
                raise RuntimeError("Internal error. Failed to initialize array when determining sort strategy")
        else:
            raise OSError("The dimensionality at k={0} or 4^k = {1} exceeds the available amount of available memory (bytes) {2}".format(self.k, N, vmem.available))
        if with_index is True:
            self.index = offsets
        self._handle.seek(0)
        self._load_block()
        assert self.kmer_ids.size == self.counts.size, "Number of Counts mismatched in size from kmer_id array size"
        assert self.frequencies.size == self.counts.size, "Number of Frequencies is mismatched in size from Count array size"
        if with_index is True:
            assert len(self.index) == self.kmer_ids.size, "Number of index elements is mismatched in size from kmer_id array size"
        self.completed = True
        return self.counts

    def slurp(self, sort:bool=False, with_index:bool=False):
        """
        A function to lazy-load an entire .kdb file into memory. 
        :param sort: Whether or not to sort the columns?
        :type sort: bool
        :param with_index: dynamically create index offsets while loading
        :type with_index: bool
        """

        if np.sum(self.counts) != 0:
            return self.counts
        else:
            return self._slurp(sort=sort, with_index=with_index)
        


    
    
class KDBWriter(bgzf.BgzfWriter):
    """
    A wrapper class around Bio.bgzf.BgzfWriter to write a .kdb file to disk.

    :ivar metadata: OrderedDict
    :ivar filename: str
    :ivar mode: str
    :ivar fileobj: io.IOBase
    :ivar compresslevel: int
    """
    
    def __init__(self, metadata:OrderedDict, filename=None, mode="w", fileobj=None, compresslevel=6):
        """Initilize the class."""
        import math
        if type(metadata) is not dict and not isinstance(metadata, OrderedDict):
            raise TypeError("kmerdb.fileutil.KDBWriter expects a valid metadata dictionary as its first positional argument")
        try:
            logger.debug("Validating .kdb metadata schema")
            jsonschema.validate(instance=dict(metadata), schema=config.kdb_metadata_schema)
            self.metadata = metadata
            self.k = self.metadata['k']
        except jsonschema.ValidationError as e:
            logger.error("kmerdb.fileutil.KDBReader couldn't validate the header YAML")
            raise e

        if fileobj:
            assert filename is None
            handle = fileobj
        else:
            if "w" not in mode.lower() and "a" not in mode.lower():
                raise ValueError("Must use write or append mode, not %r" % mode)
            if "a" in mode.lower():
                raise NotImplementedError("Append mode is not implemented yet")
            else:
                handle = _open(filename, "wb")
        self._text = "b" not in mode.lower()
        self._handle = handle
        self._buffer = b"" if "b" in mode.lower() else ""
        self.compresslevel = compresslevel
        """
        Write the header to the file
        """
        logger.info("Constructing a new .kdb file '{0}'...".format(self._handle.name))
        yaml.add_representer(OrderedDict, util.represent_yaml_from_collections_dot_OrderedDict)
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

            # 01-01-2022 This is still not a completely functional method to write data to bgzf through the Bio.bgzf.BgzfWriter class included in BioPython
            # 04-10-2023 This is still disgusting to me. I understand a limited amount of the Bio.bgzf source or the nuances of the format specification
            # That said, it's producing files, with data in the correct order. The metadata_blocks/offset calculation is still rudimentary
            # But I'm able to generate my file format reasonably.
            # I've needed to implement a basic block_writer, maintaining compatibility with the Biopython bgzf submodule.
            # 7/7/25 I like .bgzf in concept... how I'm going to do something with it is unknown. Or maybe not possible.
            # I actually thought this could be a funny little memory saving shortcut.
            # I'm not made about this, it just looks hella bloated.
            #self.write(bytes(yaml.dump(metadata, sort_keys=False), 'utf-8'))
            for i in range(self.metadata["metadata_blocks"]):
                metadata_slice = metadata_bytes[:65536]
                metadata_bytes = metadata_bytes[65536:]
                self._write_block(metadata_slice)
                #self._write_block
                self._buffer = b""
                self._handle.flush()
        elif "w" == mode.lower() or "x" == mode.lower():
            metadata_yaml_str = yaml.dump(metadata, sort_keys=False)
            metadata_bytes = bytes(metadata_yaml_str, 'utf-8')
            mode = "b"
            self._buffer = b""            
            self.write(metadata_bytes)
            self._handle.flush()
        else:
            raise RuntimeError("Could not determine proper encoding for write operations to .kdb file")

