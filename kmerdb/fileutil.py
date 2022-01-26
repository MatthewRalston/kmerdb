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
import re

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


is_integer = re.compile("^[-+]?[0-9]+$")
findall_float = re.compile(r"[-+]?(?:\d*\.\d+|\d+)")

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


    :param num: Check if a string is a float, or could be converted to a float
    :type num: str
    :returns: Whether or not the string can be parsed as a float
    :rtype: bool
    """
    if type(num) is str:
        logger.debug("Type of number being interpreted through pure Python : {0}".format(type(num)))
        findall_float.match(num)
        logger.debug(re.match(findall_float, num))
        return re.match(findall_float, num) is not None
    elif type(num) is float:
        return True
    elif type(num) is int:
        return False
    else:
        logger.error(type(num))
        #logger.error(num.dtype)
        raise TypeError("kmerdb.fileutil.isfloat expects a single str as a positional argument")


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
    elif type(line) is str and line == "":
        return None
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



def open(filepath, mode="r", metadata=None, sort:bool=False, slurp:bool=False):
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
    :returns: kmerdb.fileutil.KDBReader/kmerdb.fileutil.KDBWriter
    :rtype: kmerdb.fileutil.KDBReader
    """
    if type(filepath) is not str:
        raise TypeError("kmerdb.fileutil.open expects a str as its first positional argument")
    elif type(mode) is not str:
        raise TypeError("kmerdb.fileutil.open expects the keyword argument 'mode' to be a str")
    elif ("w" in mode or "x" in mode) and (metadata is None or not isinstance(metadata, OrderedDict)):
        raise TypeError("kmerdb.fileutil.open expects an additional metadata dictionary")
    elif type(sort) is not bool:
        raise TypeError("kmerdb.fileutil.open expects a boolean for the keyword argument 'sort'")
    elif type(slurp) is not bool:
        raise TypeError("kmerdb.fileutil.open expects a boolean for the keyword argument 'slurp'")
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
        return KDBReader(filename=filepath, mode=mode, sort=sort, slurp=slurp)
    elif "w" in mode.lower() or "x" in mode.lower():
        return KDBWriter(metadata, filename=filepath, mode=mode)
    else:
        raise ValueError("Bad mode %r" % mode)





class KDBReader(bgzf.BgzfReader):
    """
    A class that reads KDB files, potentially just for accessing header metadata, or for reading the entire contents into numpy arrays.
    

    :ivar filename: str
    :ivar fileobj: io.IOBase
    :ivar mode: str
    :ivar max_cache: int
    :ivar column_dtype: NumPy uint datatype
    :ivar count_dtypes: Numpy uint datatype
    :ivar frequencies_dtype: NumPy float datatype
    :ivar sort: bool
    :ivar slurp: bool
    """
    def __init__(self, filename:str=None, fileobj:io.IOBase=None, mode:str="r", max_cache:int=100, column_dtype:str="uint64", count_dtypes:str="uint64", frequencies_dtype:str="float64", sort:bool=False, slurp:bool=False):
        if fileobj is not None and not isinstance(fileobj, io.IOBase):
            raise TypeError("kmerdb.fileutil.KDBReader expects the keyword argument 'fileobj' to be a file object")
        elif filename is not None and type(filename) is not str:
            raise TypeError("kmerdb.fileutil.KDBReader expects the keyword argument 'filename' to be a str")
        elif type(sort) is not bool:
            raise TypeError("kmerdb.fileutil.KDBReader expects the keyword argument 'sort' to be a bool")
        
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
        self.profile      = None
        self.counts       = None
        self.frequencies  = None
        self.count_dtype = None
        self.frequencies_dtype = None
        self.kmer_ids_dtype = None
        self.profile_dtype = None
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
        logger.info("Validating the header data against the schema...")
        try:
            jsonschema.validate(instance=initial_header_data, schema=config.metadata_schema)
            self.metadata = dict(initial_header_data)
            
            self.k = self.metadata['k']
            self.kmer_ids_dtype = self.metadata["kmer_ids_dtype"]
            self.profile_dtype = self.metadata["profile_dtype"]
            self.count_dtype = self.metadata["count_dtype"]
            self.frequencies_dtype = self.metadata["frequencies_dtype"]
            self.sorted = self.metadata["sorted"]
            logger.info("Self assigning dtype to uint64 probably")
            logger.debug("Checking for metadata inference...")
            self.kmer_ids = np.zeros(4**self.metadata["k"], dtype=self.metadata["kmer_ids_dtype"])
            self.profile = np.zeros(4**self.metadata["k"], dtype=self.metadata["profile_dtype"])
            self.counts = np.zeros(4**self.metadata["k"], dtype=self.metadata["count_dtype"])
            self.frequencies = np.zeros(4**self.metadata["k"], dtype=self.metadata["frequencies_dtype"])
            logger.info("Squash target this.")
        except jsonschema.ValidationError as e:
            logger.debug(e)
            logger.error("kmerdb.fileutil.KDBReader couldn't validate the header/metadata YAML from {0} header blocks".format(num_header_blocks))
            raise e
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
            self.slurp(column_dtypes=column_dtype, count_dtypes=count_dtypes, frequencies_dtype=frequencies_dtype, sort=sort)
        self.is_int = True
        self.kmer_ids_dtype = column_dtype
        self.profile_dtype = column_dtype
        self.count_dtypes = count_dtypes
        self.frequencies_dtype = frequencies_dtype
        handle.close()
        self._handle.close()
        self._handle = None
        handle = None
        fileobj=None
        return

        #
    def read_line(self):
        """
        Returns in order, a parsed line from the .kdb file as follows:

        :returns: (kmer_id:int, count:int, metadata:dict)
        :rtype: tuple
        """

        line = self.readline()
        return parse_line(line)

    
    def __iter__(self):
        return self

    def __exit__(self, type, value, tb):
        self._handle.close()
        return
        

    def _slurp(self, column_dtypes:str="uint64", count_dtypes:str="uint64", frequencies_dtype:str="float64", sort:bool=False):
        if type(column_dtypes) is not str:
            raise TypeError("kmerdb.fileutil.KDBReader.slurp expects the column_dtypes keyword argument to be a str")
        elif type(count_dtypes) is not str:
            raise TypeError("kmerdb.fileutil.KDBReader.slurp expects the count_dtypes keyword argument to be a str")
        elif type(frequencies_dtype) is not str:
            raise TypeError("kmerdb.fileutil.KDBReader.slurp expects the frequencies_dtype keyword argument to be a str")
        try:
            np.dtype(column_dtypes)
            np.dtype(count_dtypes)
            np.dtype(frequencies_dtype)
        except TypeError as e:
            logger.error(e)
            logger.error("kmerdb.fileutil.KDBReader.slurp encountered a TypeError while assessing the numpy datatype '{0}'...".format(dtype))
            raise TypeError("kmerdb.fileutil.KDBReader._slurp expects the dtype keyword argument to be a valid numpy data type")
        if type(sort) is not bool:
            raise TypeError("kmerdb.fileutil.KDBReader._slurp expects the sort keyword argument to be a bool")
        if self._handle is None:
            self._handle = _open(self._filepath, 'rb')
        # First calculate the amount of memory required by the array
        N = 4**self.k # The dimension of the k-space, or the number of elements for the array
        num_bytes = 4 * N
        logger.info("Approximately {0} bytes".format(num_bytes))
        logger.info("Fly.")
        vmem = psutil.virtual_memory()
        self.kmer_ids = np.zeros(N, dtype=column_dtypes)
        self.profile = np.zeros(N, dtype=column_dtypes)
        self.counts = np.zeros(N, dtype=count_dtypes)
        self.frequencies = np.zeros(N, dtype=frequencies_dtype)
        self.all_kmer_metadata = []
        kmer_ids = []
        profile = []
        counts = []
        frequencies = []
        if vmem.available > num_bytes:
            if (self.profile > 0) is not all(x is False for x in (self.profile > 0) ):
                # Do the slurp
                logger.info("Reading profile into memory")
                if sort is False:
                    i = 0
                    try:
                        for j in range(N):
                            #logger.debug("Reading {0}th line...".format(j))
                            try:
                                line = next(self)
                            except StopIteration as e:
                                logger.error("Finished loading initial profile through slurp-on-init")
                                raise e
                            if line is None:
                                logger.warning("Next was None... profile was sparse, breaking")
                                sys.exit(1)
                                break
                            # Don't forget to not parse the metadata column [:-1]

                            x, p, kmer_id, _count, _frequency, kmer_metadata = line.rstrip().split("\t")
                            x = int(x)
                            p = int(p)
                            
                            logger.debug(line.rstrip())


                            assert j == x, "Line number did not match"

                            logger.debug("{0}\t{1}\t{2}\t'{3}'\t'{4}'".format(x, p, kmer_id, _count, _frequency))

                            kmer_id = int(kmer_id)
                            kmer_ids.append(kmer_id)
                            profile.append(kmer_id)

                            if isfloat(_count):
                                parsed_integer_string = is_integer.match(_count)
                                if parsed_integer_string is not None:
                                    logger.info("Casting count field to integer...")
                                    logger.info("Starting count | castin' as int '{0}'".format(_count))

                                    count = int(_count)
                                    #logger.warning("Casting count field to integer.")
                                    #raise TypeError("Shouldn't be normal usage. Unsupported")
                                else:
                                    logger.info("Starting count | cast as float '{0}'".format(_count))
                                    count = float(_count)
                                logger.info("Determining the type of data in the file1...")
                                logger.info("File counts at {0}, file frequency = {1}".format(count, _frequency))
                            else:
                                parsed_integer_string = is_integer.match(_count)
                                if parsed_integer_string is not None:
                                    logger.info("Casting count field to integer...")
                                    logger.info("Starting count | castin' as int '{0}'".format(_count))
                                    count = int(_count)
                                    #logger.warning("Casting count field to integer.")
                                    #raise TypeError("Shouldn't be normal usage. Unsupported")
                                else:
                                    logger.info("Starting count | IDGAF castin' as int '{0}'".format(_count))
                                    count = int(_count)
                                logger.info("Determining the type of data in the file2...")
                                logger.info("File counts at {0}, file frequency = '{1}'".format(count, _frequency))
                            assert type(count) is int, "_slurp | type of count field is not int"
                            frequency = float(count)/N
                            self.profile[j] = p


                            self.profile[j] = p
                            self.kmer_ids[p] = kmer_id
                            self.counts[kmer_id] = count
                            self.frequencies[kmer_id] = _frequency
                            self.all_kmer_metadata.append(kmer_metadata)
                            #sys.stderr.write("::DEBUG::   |\-)(||||..... KMER_ID: {0} COUNT: {1}".format(kmer_id, count))
                            try:
                                if isfloat(_frequency):
                                    logger.debug("Frequency: '{0}'".format(_frequency))
                                    try:
                                        logger.debug("type: {0}".format(type(_frequency)))
                                        logger.debug("Using integer matching regex...")
                                        pis = is_integer.match(_frequency)
                                    except TypeError as e:
                                        logger.error(e)
                                        raise e
                                    if pis is not None:
                                        logger.warning("Matched frequency as integer: perhaps an error or a genuine zero...Reseting to zero")
                                        frequency = float(0)
                                        logger.error("Starting frequency cast as fresh float 0.0 | '{0}'".format(_frequency))
                                        #raise TypeError("Matched frequency to integer instead of float. Aborting.")
                                    else:
                                        logger.info("Starting frequency cast as float | '{0}'")
                                        frequency = float(_frequency)
                                    logger.info("Determining type of data in the file3...")
                                    logger.info("It looks like the frequency data is a double. interpretted frequency as {0}, frequency string: '{1}'".format(frequency, _frequency))
                                    assert frequency == float(_frequency), "kmerdb.fileutil.KDBReader._slurp | Frequency did not match expected value based on the count during recalculation of frequencies..."

                                else:
                                    parsed_integer_string = is_integer.match(_frequency)
                                    if parsed_integer_string is not None:
                                        logger.warning("Casting frequency field to integer.")
                                        raise TypeError("Casting frequency field to integer is unsupported.")
                                    else:
                                        logger.info("Starting frequency | cast as int '{0}'".format(_frequency))
                                        logger.info("Castin' frequency field to an integer.")
                                        
                                        frequency = int(count)
                                    logger.info("Determining type of data in the file4...")
                                    logger.info("It looks like the frequency data is an int. python frequency as {0}, frequency string: '{1}'".format(frequency, _frequency))
                                    #logger.info("It looks like the frequency data is a double. Caclulated frequency at {0}, frequency read: {1}".format(frequency, _frequency))
                                    assert float(frequency) == float(_frequency), "kmerdb.fileutil.KDBReader._slurp | Frequency did not match expected value based on the count during recalculation of frequencies..."


                            except AssertionError as e:
                                logger.error(e)
                                logger.warning("Interpretting string for NumPy formatting")
                                raise e
                            assert type(frequency) is float, "_slurp | frequency field not cast as float"
                            
                            # if self.sorted is True:
                            #     self.profile[]
                            # else:
                            #     self.profile[j] = j
                            #     self.kmer_ids[j] = kmer_id
                            #     self.counts[kmer_id] = count
                            #     self.frequencies[kmer_id] = frequency
                        logger.debug("The {0}th line was kmer-id: {1} with an abundance of {2}".format(j, kmer_id, count))
                        i += 1

                        self.profile[j] = p
                        self.kmer_ids[j] = kmer_id
                        self.counts[kmer_id] = count
                        self.frequencies[kmer_id] = frequency
                        #self.all_kmer_metadata[i] = kmer_metadata)
                    except StopIteration as e:
                        if i == N:
                            logger.debug("Read {0} lines from the file...".format(i))
                            logger.warning("StopIteration was raised!!!")
                            raise e
                        else:
                            logger.error("Number of lines read: {0}".format(i))
                            logger.error("Number of k-mers (N) set globally.".format(N))
                            logger.error("Read only {0} lines from the file...".format(i))
                            logger.error("Profile must have been read before this point")
                            logger.error(e)
                            raise e
                    except Exception as e:
                        logger.error(e)
                        raise e
                    logger.debug("Validating shape prior to if unsorted is True")
                    logger.debug("Number of lines read: {0}".format(i))
                    logger.debug("Number of k-mer (N) set globally: {0}".format(N))
                    logger.debug("Profile shape: {0}".format(self.profile.size))
                    logger.debug("Kmer_ids shape: {0}".format(self.kmer_ids.size))
                    logger.debug("Counts shape: {0}".format(self.counts.size))
                    logger.debug("Frequencies shape: {0}".format(self.frequencies.size))
                    assert self.profile.size == self.kmer_ids.size, "Profile size has diverged from kmer_ids shape"
                    assert self.profile.size == self.counts.size, "Counts size has diverged from profile size/shape"
                    assert self.counts.size == self.frequencies.size, "Frequencies size has diverged from counts shape"
                elif sort is True:
                    i = 0
                    try:
                        for j in range(N):
                            #logger.debug("Reading {0}th line...".format(j))
                            try:
                                line = next(self)
                            except StopIteration as e:
                                logger.error("Finished loading initial profile through slurp-on-init")
                                raise e
                            if line is None:
                                logger.warning("Next was None... profile was sparse, breaking")
                                sys.exit(1)
                                break
                            # Don't forget to not parse the metadata column [:-1]

                            x, p, kmer_id, _count, _frequency, metadata = line.rstrip().split("\t")

                            logger.debug(line.rstrip())

                            logger.debug("{0}\t{1}\t{2}\t{3}".format(j, kmer_id, _count, _frequency))

                            _frequency = float(_frequency)
                            p = int(p)
                            kmer_id = int(kmer_id)

                            if isfloat(_count):
                                s = findall_float.match(_count)
                                if s is not None:
                                    count = float(_count)
                                else:
                                    raise TypeError("kmerdb.fileutil.KDBReader._slurp | couldn't properly detect float")
                                count = float(_count)
                                logger.info("Determining the type of data in the file1...")
                                logger.info("File counts at {0}, file frequency = {1}".format(count, _frequency))
                            else:
                                logger.info("Determining the type of data in the file2...")
                                logger.info("File counts at {0}, file frequency = {1}".format(count, _frequency))
                                logger.debug("kmer_id: {0}, count: {1}, frequency: {2}".format(kmer_id, count, _frequency))

                                s = is_integer.match(_count)
                                if s is not None:
                                    count = int(_count)
                                else:
                                    raise TypeError("kmerdb.fileutil.KDBReader._slurp | couldn't properly detect integer")

                            frequency = float(count)/N
                            # Trying to get the needle threaded here.
                            # The goal is for the profile to be a running index, similar to the row number
                            # But this could maintain lexical sort order
                            # 
                            self.profile[j] = p
                            self.kmer_ids[p] = kmer_id
                            self.counts[kmer_id] = count
                            self.frequencies[kmer_id] = _frequency
                            i+=1
                            #sys.stderr.write("::DEBUG::   |\-)(||||..... KMER_ID: {0} COUNT: {1}".format(kmer_id, count))
                            try:
                                if isfloat(_frequency):
                                    logger.info("Determining type of data in the file3...")
                                    logger.info("It looks like the frequency data is a double. interpretted frequency as {0}, frequency string: '{1}'".format(frequency, _frequency))
                                    assert float(frequency) == float(_frequency), "Frequency did not match expected value based on the count..."
                                else:
                                    logger.info("Determining type of data in the file4...")
                                    logger.info("It looks like the frequency data is a double. python frequency as {0}, frequency string: '{1}'".format(frequency, _frequency))
                                    #logger.info("It looks like the frequency data is a double. Caclulated frequency at {0}, frequency read: {1}".format(frequency, _frequency))
                                    assert float(frequency) == float(_frequency), "Frequency did not match expected value based on the count..."

                            except AssertionError as e:
                                logger.error(e)
                                logger.warning("Interpretting string for NumPy formatting")
                                #raise e
                            # if self.sorted is True:
                            #     self.profile[]
                            # else:
                            #     self.profile[j] = j
                            #     self.kmer_ids[j] = kmer_id
                            #     self.counts[kmer_id] = count
                            #     self.frequencies[kmer_id] = frequency
                        logger.debug("The {0}th line was kmer-id: {1} with an abundance of {2}".format(j, kmer_id, count))


                        logger.debug("Validating shape prior to if unsorted is True")
                        logger.debug("Number of lines read: {0}".format(i))
                        logger.debug("Number of k-mer (N) set globally: {0}".format(N))
                        logger.debug("Profile shape: {0}".format(self.profile.size))
                        logger.debug("Kmer_ids shape: {0}".format(self.kmer_ids.size))
                        logger.debug("Counts shape: {0}".format(self.counts.size))
                        logger.debug("Frequencies shape: {0}".format(self.frequencies.size))
                        assert self.profile.size == self.kmer_ids.size, "Profile size has diverged from kmer_ids shape"
                        assert self.profile.size == self.counts.size, "Counts size has diverged from profile size/shape"
                        assert self.counts.size == self.frequencies.size, "Frequencies size has diverged from counts shape"

                    except StopIteration as e:
                        if i == N:
                            logger.debug("Read {0} lines from the file...".format(i))
                            logger.warning("StopIteration was raised!!!")
                            raise e
                        else:
                            logger.error("Number of lines read: {0}".format(i))
                            logger.error("Number of k-mers (N) set globally.".format(N))
                            logger.error("Read only {0} lines from the file...".format(i))
                            logger.error("Profile must have been read before this point")
                            logger.error(e)
                            raise e
                    except Exception as e:
                        logger.error(e)
                        raise e
                logger.info("Read {0} lines from the file...".format(i))
                logger.debug("Rechecking shape... I can't help myself.")
                assert self.profile.size == self.kmer_ids.size, "Number of k-mer ids mismatched in count from row-index 'profile'"
                assert self.kmer_ids.size == self.counts.size, "Number of Counts is mismatched in count from row-index 'profile'"
                assert self.frequencies.size == self.counts.size, "Number of Frequencies is mismatched in count from row-index 'profile'"
                logger.debug("Check completed...")            
                self._handle.seek(0)
                self._load_block()

                logger.debug("Dammit, why can't i reset the Bio.bgzf filehandle...")
                if sort is True:
                    # If the file is sorted, do not sort
                    kmer_ids_sorted_by_count = np.lexsort((self.kmer_ids, self.counts))
                    reverse_kmer_ids_sorted_by_count = np.flipud(kmer_ids_sorted_by_count)
                    for i, idx in enumerate(kmer_ids_sorted_by_count): # This is right, not fixing this.
                        p = self.profile[i]
                        kmer_id = self.kmer_ids[p]
                        count = self.counts[idx]

                        logger.debug("{0}\t{1}\t{2}\t{3}".format(i, p,  kmer_id, count))
                        # I stand corrected
                        #self.profile[i] = profile[i]
                        #self.kmer_ids[i] = kmer_ids[profile[i]]
                        #self.counts[idx] = counts[idx]
                        ###self.frequencies[idx] = frequencies[idx]
                        logger.debug("Just in casey eggs and bakey...")

            else:
                logger.debug("Profile has been intiailized as zeros only, incorrect instantiation.")
                raise RuntimeError("Whoops, something went wrong there")
                sys.exit(1)
                                    
            if self.profile is not None:
                logger.info("Data types | profile: {0}, kmer_ids: {1}, counts: {2}, frequencies: {3}".format(type(self.profile[0]), type(self.kmer_ids[0]), type(self.counts[0]), type(self.frequencies[0])))

            else:
                logger.info("Data types | profile: {0}, kmer_ids: {1}, counts: {2}, frequencies: {3}".format(type(self.profile[0]), type(self.kmer_ids[0]), type(self.counts[0]), type(self.frequencies[0])))
                logger.warning("Profile, dickhead, it's already in memory. Dont rerun.")
                logger.warning("Profile is already loaded in memory! Did you remember to deallocate it when you were done?")
                #return self.profile
        else:
            raise OSError("The dimensionality at k={0} or 4^k = {1} exceeds the available amount of available memory (bytes) {2}".format(self.k, N, vmem.available))
        if self.profile.dtype != self.profile_dtype:
            raise TypeError("kmerdb.fileutil.KDBReader._slurp encountered an awful TypeError")
        elif self.kmer_ids.dtype != self.kmer_ids_dtype:
            raise TypeError("kmerdb.fileutil.KDBReader._slurp encountered an awful TypeError")
        elif self.counts.dtype != self.count_dtype:
            raise TypeError("kmerdb.fileutil.KDBReader._slurp encountered an awful TypeError")
        elif self.frequencies.dtype != self.frequencies_dtype:
            raise TypeError("kmerdb.fileutil.KDBReader._slurp encountered an awful TypeError")
        #self.profile = np.array(counts, dtype=suggested_dtype)
        #self.kmer_ids = np.array(kmer_ids, dtype=suggested_dtype)
        #self.kmer_ids = np.array(kmer_ids, dtype=column_dtypes)
        #self.profile = np.array(profile, dtype=column_dtypes)
        #self.counts = np.array(counts, dtype=count_dtypes)
        #self.frequencies = np.array(frequencies, dtype=frequencies_dtype)
        self._handle.seek(0)
        self._load_block()
        #print([x for x in np.ndindex(self.kmer_ids.flat) if x < ])
        logger.info("Data types read from file")
        logger.info(type(self.profile[0]))
        logger.info(type(self.kmer_ids[0]))
        logger.info(type(self.counts[0]))
        logger.info(type(self.frequencies[0]))


        logger.info("K-mer profile size (row-index): {0}".format(self.profile.size))
        logger.info("Kmer_ids size: {0}".format(self.kmer_ids.size))
        logger.info("Counts size: {0}".format(self.counts.size))
        logger.info("Frequencies size: {0}".format(self.frequencies.size))
        assert self.profile.size == self.kmer_ids.size, "Number of k-mer ids mismatched in count from row-index 'profile'"
        assert self.kmer_ids.size == self.counts.size, "Number of Counts is mismatched in count from row-index 'profile'"
        assert self.frequencies.size == self.counts.size, "Number of Frequencies is mismatched in count from row-index 'profile'"
        
        return self.counts

    def slurp(self, column_dtypes:str="uint64", count_dtypes:str="uint64", frequencies_dtype:str="float64", sort:bool=False):
        """
        A function to lazy-load an entire .kdb file into memory. 

        :param column_dtypes: a NumPy uint datatype
        :type column_dtypes: str
        :param count_dtypes: a NumPy uint datatype
        :type count_dtypes: str
        :param frequencies_dtype: a NumPy float datatype
        :type frequencies_dtype: str
        :param sort: Whether or not to sort the columns?
        :type sort: bool
        """

        if np.sum(self.counts) != 0:
            return self.counts
        else:
            return self._slurp(column_dtypes=column_dtypes, count_dtypes=count_dtypes, frequencies_dtype=frequencies_dtype, sort=sort)

    
    
class KDBWriter(bgzf.BgzfWriter):
    """
    A wrapper class around Bio.bgzf.BgzfWriter to write a kdb file to disk.

    :ivar metadata: OrderedDict
    :ivar filename: str
    :ivar mode: str
    :ivar fileobj: io.IOBase
    :ivar compresslevel: int
    """
    
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



class FileReader:

    def load_file(self, f):
        self.file = open(f, 'r', slurp=False)
        return self.file

    def slurp(self):
        return self.file.slurp()
        
