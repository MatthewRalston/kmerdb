import io
import sys
import os
import gzip
import tempfile
import yaml, json
from collections import deque

from builtins import open as _open

import jsonschema
from Bio import SeqIO, bgzf
import boto3
sys.path.append('..')

import config
from kdb import kmer, database

# Logging configuration
import logging
logger = logging.getLogger(__file__)
# S3 configuration
s3 = boto3.resource('s3')
s3client = boto3.client('s3')
s3prefix = "s3://"


def _s3_file_download(self, seqpath, temporary=True):
    """
    Note: the file will be downloaded into a temporary file that needs to be deleted afterwards
    It will create the temporary file with respect to the TMP bash variable 'export TMP=/some/temporary/location'
    :param seqpath: The s3 identifier of a object. 's3://bucket/example.fasta'
    :type seqpath: str
    :returns: The location of a downloaded gennomic Fasta file
    :rtype: str
    """
    if type(seqpath) is not str:
        raise TypeError("kdb.fileutil.SeqReader.__s3_file_download expects a str 'seqpath' as its first positional argument")
    elif seqpath[0:5] != s3prefix:
        raise TypeError("kdb.fileutil.SeqReader.__s3_file_download expects a s3 object reference its first positional argument. e.g. 's3://bucket/example.txt'")
    elif type(temporary) is not bool:
        raise TypeError("kdb.fileutil.SeqReader.__s3_file_download expects the keyword argument temporary to be a bool")
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


def open(filepath, mode="r"):
    if type(filepath) is not str:
        raise TypeError("kdb.fileutil.open expects a str as its first positional argument")
    elif type(mode) is not str:
        raise TypeError("kdb.fileutil.open expects the keyword argument 'mode' to be a str")
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
        return KDBReader(filepath, mode=mode)
    elif "w" in mode.lower() or "a" in mode.lower():
        return KDBWriter(filepath, mode=mode)
    else:
        raise ValueError("Bad mode %r" % mode)




class KDBReader(bgzf.BgzfReader):
    def __init__(self, filename:str=None, fileobj:io.IOBase=None, mode:str="r", max_cache:int=100):
        if fileobj is not None and not isinstance(fileobj, io.IOBase):
            raise TypeError("kdb.fileutil.KDBReader expects the keyword argument 'fileobj' to be a file object")
        elif filename is not None and type(filename) is not str:
            raise TypeError("kdb.fileutil.KDBReader expects the keyword argument 'filename' to be a str")
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
        self._load_block(handle.tell())
        header_data = yaml.safe_load(self._buffer)
        try:
            jsonschema.validate(instance=header_data, schema=config.header_schema)
            self.header = header_data
            self.k = self.header['k']
        except jsonschema.ValidationError as e:
            logger.debug(e)
            logger.error("kdb.fileutil.KDBReader couldn't validate the header YAML")
            raise e

        self._offsets = deque()
        for values in bgzf.BgzfBlocks(self._handle):
            #logger.debug("Raw start %i, raw length %i, data start %i, data length %i" % values)
            self._offsets.appendleft(values) # raw start, raw length, data start, data length
        if len(self._offsets) == 0:
            raise IOError("kdb.fileutil.KDBReader opened an empty file")
        # Skip the zeroth block
        self._load_block()
        # print(str(self._buffer)) # 1
        # print(self.readline())
        # self._load_block()
        # print(self._buffer) # 2
        # print(self.readline())

        #

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
    
    # def __next__(self):
    #     try:
    #         self._load_block()
    #     except StopIteration as e:
    #         logger.warning(e)
    #         self._handle.close()
    #         raise StopIteration
    #     return self._buffer
    



class KDBWriter(bgzf.BgzfWriter):
    def __init__(self, header:dict, filename=None, mode="w", fileobj=None, compresslevel=6):
        """Initilize the class."""
        if type(header) is not dict:
            raise TypeError("kdb.fileutil.KDBWriter expects a valid header object as its first positional argument")
        try:
            jsonschema.validate(instance=header, schema=config.header_schema)
            self.header = header
            self.k = self.header['k']
        except jsonschema.ValidationError as e:
            logger.debug(e)
            logger.error("kdb.fileutil.KDBReader couldn't validate the header YAML")
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
        self._buffer = b""
        self.compresslevel = compresslevel

        self._write_block(bgzf._as_bytes(yaml.dump(self.header)))
        self._buffer = b""
        self._handle.flush()
        
    def write_block(self, recs):
        if type(recs) is not list:
            raise TypeError("kdb.fileutil.KDBWriter expects a str as its first positional argument")
        self._block_size = bgzf._as_bytes(self._buffer).__sizeof__()
        while len(recs):
            while self._block_size < 65530:
                newline = "\t".join(recs.pop()) + "\n"

                if self._block_size + bgzf.as_bytes(newline).__sizeof__() < 65530:
                    self._buffer += newline
                    self._block_size += bgzf._as_bytes(newline).__sizeof__()
                else:
                    self._write_block(bgzf._as_bytes(self._buffer))
                    self._buffer = b""
                    self._handle.flush()
                    self._buffer = newline
                    self._block_size = bgzf._as_bytes(newline).__sizeof__()
        self._write_block(bgzf._as_bytes(self._buffer))
        self._buffer = b""
        self._handle.flush()
        
