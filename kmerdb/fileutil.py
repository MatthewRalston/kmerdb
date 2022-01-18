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
import gzip
import tempfile
import yaml, json
from collections import deque, OrderedDict
import psutil
import numpy as np
import math

#import pdb

from builtins import open as _open

import jsonschema
from Bio import SeqIO, bgzf

#sys.path.append('..')

from kmerdb import kmer, util, config

# Logging configuration
import logging
logger = logging.getLogger(__file__)
# S3 configuration
s3prefix = "s3://"


def is_all_fasta(filenames):
    """Tests if all the strings in a list are fasta format.
    
    :param filenames:
    :type list:
    :returns bool:
    :rtype bool:
    """
    if type(filenames) is not list:
        raise TypeError("kmerdb.fileutil.is_all_fasta() expects a list as its first positional argument")
    elif not all(type(s) is str for s in filenames):
        raise TypeError("kmerdb.fileutil.is_all_fasta() expects a list of str as its first positional argument")
    return all((f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.fa.gz') or f.endswith('.fasta.gz') or f.endswith('.fna') or f.endswith('.fna.gz')) for f in filenames)


def isfloat(num):
    """
    Thanks to the author of:
    https://www.programiz.com/python-programming/examples/check-string-number
    """
    try:
        float(num)
        return True
    except ValueError:
        return False


def _s3_file_download(self, seqpath, temporary=True):
    """
    Note: the file will be downloaded into a temporary file that needs to be deleted afterwards
    It will create the temporary file with respect to the TMP bash variable 'export TMP=/some/temporary/location'
    :param seqpath: The s3 identifier of a object. 's3://bucket/example.fasta'
    :type seqpath: str
    :returns: The location of a downloaded gennomic Fasta file
    :rtype: str
    """
    import boto3
    s3 = boto3.resource('s3')
    s3client = boto3.client('s3')



    if type(seqpath) is not str:
        raise TypeError("kmerdb.fileutil.SeqReader.__s3_file_download expects a str 'seqpath' as its first positional argument")
    elif seqpath[0:5] != s3prefix:
        raise TypeError("kmerdb.fileutil.SeqReader.__s3_file_download expects a s3 object reference its first positional argument. e.g. 's3://bucket/example.txt'")
    elif type(temporary) is not bool:
        raise TypeError("kmerdb.fileutil.SeqReader.__s3_file_download expects the keyword argument temporary to be a bool")
    seqpath = seqpath.lstrip(s3prefix)
    pathsegs =seqpath.split('/')
    bucket = pathsegs.pop(0)
    fname = os.path.basename(seqpath)
    key = '/'.join(pathsegs)
    if seqpath[-3:] == ".gz":
        suffix = '.' + '.'.join(seqpath.split('.')[-2:])
    else:
        suffix = path.splitext(seqpath)[1]
    if temporary is True:
        filepath = tempfile.NamedTemporaryFile(mode='w+b', suffix=suffix, delete=False)
        logger.info("Downloading '{0}' => '{1}'...".format(seqpath, filepath.name))
    else:
        filepath = open(fname, 'w+b')
        logger.info("Downloading '{0}' => '{1}'...".format(seqpath, fname))
    obj = s3.Object(bucket, key)
    obj.download_fileobj(filepath)
    filepath.close()
    return filepath.name

def parse_line(line):
    """
    Parses a line according to the expected syntax, and returns the python data types expected as a tuple.

    :param line:
    :type line: str
    :returns: kmer_id, count, metadata
    :rtype: tuple
    """
    
    if type(line) is not str:
        raise TypeError("kmerdb.fileutil.parse_line expects to a str as its first positional argument")
    else:
        linesplit = line.rstrip().split("\t")
        if len(linesplit) != 3:
            logger.error("Full line:\n{0}".format(line))
            raise ValueError("kmerdb.fileutil.parse_line() encountered a .kdb line without 3 columns, a violation of the format")
        else:
            kmer_id, count, kmer_metadata = linesplit
            if isfloat(count):
                kmer_id, count = int(kmer_id), float(count)
            else:
                kmer_id, count = int(kmer_id), int(count)
            kmer_metadata = yaml.safe_load(kmer_metadata)
            if type(kmer_metadata) is dict:
                return kmer_id, count, kmer_metadata
            else:
                logger.error("Improperly formatted k-mer metadata field")
                logger.error(line)
                raise ValueError("kmerdb.fileutil.parse_line(): Improperly formatted k-mer metadata field")



def open(filepath, mode="r", metadata=None):
    """
    Opens a file for reading or writing. Valid modes are 'xrwbt'. 'metadata=' is needed when writing/creating.

    :param filepath:
    :type filepath: str
    :param mode:
    :type mode: str
    :param metadata: The file header/metadata dictionary to write to the file.
    :type metadata: dict
    :returns: kmerdb.fileutil.KDBReader/kmerdb.fileutil.KDBWriter
    :rtype: kmerdb.fileutil.KDBReader
    """
    if type(filepath) is not str:
        raise TypeError("kmerdb.fileutil.open expects a str as its first positional argument")
    elif type(mode) is not str:
        raise TypeError("kmerdb.fileutil.open expects the keyword argument 'mode' to be a str")
    elif ("w" in mode or "x" in mode) and (metadata is None or not isinstance(metadata, OrderedDict)):
        raise TypeError("kmerdb.fileutil.open expects an additional metadata dictionary")
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
        return KDBReader(filename=filepath, mode=mode)
    elif "w" in mode.lower() or "x" in mode.lower():
        return KDBWriter(metadata, filename=filepath, mode=mode)
    else:
        raise ValueError("Bad mode %r" % mode)





class KDBReader(bgzf.BgzfReader):
    def __init__(self, filename:str=None, fileobj:io.IOBase=None, mode:str="r", max_cache:int=100, dtype:str="uint32"):
        if fileobj is not None and not isinstance(fileobj, io.IOBase):
            raise TypeError("kmerdb.fileutil.KDBReader expects the keyword argument 'fileobj' to be a file object")
        elif filename is not None and type(filename) is not str:
            raise TypeError("kmerdb.fileutil.KDBReader expects the keyword argument 'filename' to be a str")
        elif type(dtype) is not str:
            raise TypeError("kmerdb.fileutil.KDBReader expects the keyword argument 'dtype' to be a str")

        try:
            np.dtype(dtype)
        except TypeError as e:
            raise TypeError("kmerdb.fileutil.KDBReader expects a valid NumPy data type as an optional argument")
        
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
        self.profile      = None
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
            elif "dtype" not in initial_header_data.keys() or initial_header_data["dtype"] not in ["uint32", "uint64", "float32", "float64"]:
                logger.error(initial_header_data)
                raise TypeError("kmerdb.fileutil.KDBReader couldn't determine the type of count/frequency stored in the file.")
            else:
                logger.debug(initial_header_data)
                if initial_header_data["metadata_blocks"] == 1:
                    logger.info("1 metadata block: Not loading any additional blocks")
                else:
                    for i in range(initial_header_data["metadata_blocks"] - 1):
                        logger.info("Multiple metadata blocks read...")
                        self._load_block(self._handle.tell())
                        logger.debug("Secondary blocks read")
                        addtl_header_data = yaml.safe_load(self._buffer.rstrip(config.header_delimiter))
                        if type(addtl_header_data) is str:
                            logger.error(addtl_header_data)
                            logger.error("Couldn't determine yo' block.::/")
                            #logger.error("Couldn't at you goldy.")
                            logger.error("Sorry - Chris Brown.")
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
        logger.debug("Reading additional blocks as YAML...")
        sys.stderr.write("\n")
        logger.info("\n\n\nFly\n\n\n")
        logger.info("Validating the header data against the schema...")
        try:
            jsonschema.validate(instance=initial_header_data, schema=config.metadata_schema)
            self.metadata = initial_header_data
            self.k = self.metadata['k']
            self.dtype = self.metadata["dtype"]
            #logger.info("Squash target my nuts...")
            logger.info("Self assigning dtype to uint32")
            #logger.debug("Checking for metadata inference...")
            #print(self.dtype)
            #sys.exit(1)

            
        except jsonschema.ValidationError as e:
            logger.debug(e)
            logger.error("kmerdb.fileutil.KDBReader couldn't validate the header/metadata YAML from {0} header blocks".format(num_header_blocks))
            raise e
        self.metadata["header_offset"] = self._handle.tell()
        logger.debug("Handle set to {0} after reading header, saving as handle offset".format(self.metadata["header_offset"]))
        self._reader = gzip.open(self._filepath, 'r')
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

        if dtype != self.dtype:
            raise TypeError("kmerdb.fileutil.KDBReader._slurp expects the dtype keyword argument to be equal to the dtype in the file. Got {0} was {1}".format(dtype, self.dtype))

        self.slurp(dtype=dtype)
        self.is_int = False
        self.dtype = dtype
        

        #
    def read_line(self):
        """
        Returns in order, a parsed line from the .kdb file as follows:

        :returns: (kmer_id:int, count:int, metadata:dict)
        :rtype: tuple
        """

        line = self.readline()
        return parse_line(line)

    # def _load_block(self, start_offset=None):
    #     if start_offset is None:
    #         # If the file is being read sequentially, then _handle.tell()
    #         # should be pointing at the start of the next block.
    #         # However, if seek has been used, we can't assume that.
    #         start_offset = self._block_start_offset + self._block_raw_length
    #     if start_offset == self._block_start_offset:
    #         self._within_block_offset = 0
    #         return
    #     elif start_offset in self._buffers:
    #         # Already in cache
    #         self._buffer, self._block_raw_length = self._buffers[start_offset]
    #         self._within_block_offset = 0
    #         self._block_start_offset = start_offset
    #         return
    #     # Must hit the disk... first check cache limits,
    #     while len(self._buffers) >= self.max_cache:
    #         # TODO - Implemente LRU cache removal?
    #         self._buffers.popitem()
    #     # Now load the block
    #     handle = self._handle
    #     if start_offset is not None:
    #         handle.seek(start_offset)
    #     self._block_start_offset = handle.tell()
        
    #     block_size, self._buffer = bgzf._load_bgzf_block(handle, self._text)
    #     self._within_block_offset = 0
    #     self._block_raw_length = block_size
    #     # Finally save the block in our cache,
    #     self._buffers[self._block_start_offset] = self._buffer, block_size
        
    
    def __iter__(self):
        return self

    def __next__(self):
        line = self.readline()
        if len(line):
            return line
        else:
            raise StopIteration

    def __exit__(self, type, value, tb):
        self._handle.close()
        self._reader.close()
        return
        
        
    # def __next__(self):
    #     try:
    #         self._load_block()
    #     except StopIteration as e:
    #         logger.warning(e)
    #         self._handle.close()
    #         raise StopIteration
    #     return self._buffer

    def _slurp(self, dtype:str="uint32"):
        """
        A function to read an entire .kdb file into memory
        """
        if type(dtype) is not str:
            raise TypeError("kmerdb.fileutil.KDBReader.slurp expects the dtype keyword argument to be a str")
        try:
            np.dtype(dtype)
        except TypeError as e:
            logger.error(e)
            logger.error("kmerdb.fileutil.KDBReader.slurp encountered a TypeError while assessing the numpy datatype '{0}'...".format(dtype))
            raise TypeError("kmerdb.fileutil.KDBReader._slurp expects the dtype keyword argument to be a valid numpy data type")
        if dtype != self.dtype:
            raise TypeError("kmerdb.fileutil.KDBReader._slurp expects the dtype keyword argument to be equal to the dtype in the file. Got {0} was {1}".format(dtype, self.dtype))
        # First calculate the amount of memory required by the array
        N = 4**self.k # The dimension of the k-space, or the number of elements for the array
        num_bytes = 4 * N
        logger.info("Approximately {0} bytes".format(num_bytes))
        logger.info("Fly.")
        vmem = psutil.virtual_memory()
        if vmem.available > num_bytes:
            if self.profile is None:
                # Do the slurp
                i = 0
                try:
                    self.profile = np.zeros(4**self.k, dtype=dtype)

                    for j in range(N):
                        #logger.debug("Reading {0}th line...".format(j))
                        line = next(self)
                        if line is None:
                            logger.warning("Next was None... profile was sparse, breaking")
                            sys.exit(1)
                            break
                        # Don't forget to not parse the metadata column [:-1]

                        kmer_id, _count, metadata = line.rstrip().split("\t")
                        kmer_id = int(kmer_id)
                        if isfloat(_count):
                            count = float(_count)
                            self.is_int32 = False
                            self.is_float32 = True
                        else:
                            count = int(_count)
                            self.is_int32 = True
                            self.is_float32 = False
                        #logger.debug("The {0}th line was kmer-id: {1} with an abundance of {2}".format(j, kmer_id, count))
                        i += 1
                        self.profile[kmer_id] = count
                    logger.info("Read {0} lines from the file...".format(i))
                    return self.profile
                except StopIteration as e:
                    if i == N:
                        logger.debug("Read {0} lines from the file...".format(i))
                        logger.warning("StopIteration was raised!!!")
                        raise e
                    else:
                        logger.debug("Read only {0} lines from the file...".format(i))
                        logger.debug("Profile must have been sparse...")
                        raise e
            else:
                logger.warning("Profile, dickhead, it's already in memory. Dont rerun.")
                logger.warning("Profile is already loaded in memory! Did you remember to deallocate it when you were done?")
                return self.profile
        else:
            raise OSError("The dimensionality at k={0} or 4^k = {1} exceeds the available amount of available memory (bytes) {2}".format(self.k, N, vmem.available))
        if self.profile.dtype != dtype:
            raise TypeError("kmerdb.fileutil.KDBReader._slurp expects the dtype keyword argument to be the inferred numpy data type")
        self.dtype = dtype
        return self.profile

    def slurp(self, dtype:str="uint32"):
        self._slurp(dtype=dtype)

    
    
class KDBWriter(bgzf.BgzfWriter):
    def __init__(self, metadata:OrderedDict, filename=None, mode="w", fileobj=None, compresslevel=6):
        """Initilize the class."""
        if not isinstance(metadata, OrderedDict):
            raise TypeError("kmerdb.fileutil.KDBWriter expects a valid metadata object as its first positional argument")
        try:
            logger.debug("Validating metadata schema against the config.py header schema")
            jsonschema.validate(instance=dict(metadata), schema=config.metadata_schema)
            self.metadata = metadata
            self.k = self.metadata['k']
        except jsonschema.ValidationError as e:
            logger.debug(e)
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



        
        logger.info("Constructing a new .kdb file '{0}'...".format(self._handle.name))
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
        elif "w" == mode.lower() or "x" == mode.lower():
            self.write(yaml.dump(metadata, sort_keys=False))
            self._buffer = ""
            self._handle.flush()
        else:
            logger.error("Mode: {}".format(mode.lower()))
            raise RuntimeError("Could not determine proper encoding for write operations to .kdb file")

