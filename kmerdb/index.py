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


import sys
import os
import copy
import gzip
import yaml
import numpy as np
#import jsonschema
sys.path.append('..')

from kmerdb import fileutil


import logging
logger = logging.getLogger(__file__)



def open(filepath:str=None, indexfilepath:str=None, idx:list=None, mode:str="u", force:bool=False, DEBUG:bool=False):
    if filepath is not None and type(filepath) is not str:
        raise TypeError("kmerdb.index.open expects a str for keyword argument 'filepath'")
    elif indexfilepath is not None and type(indexfilepath) is not str:
        raise TypeError("kmerdb.index.open expects a str for argument 'indexfilepath'")
    elif mode is None or type(mode) is not str:
        raise TypeError("kmerdb.index.open expects a str for argument 'mode'")
    elif force is not None and type(force) is not bool:
        raise TypeError("kmerdb.index.open expects a bool for argument 'force'")
    elif DEBUG is not None and type(DEBUG) is not bool:
        raise TypeError("kmerdb.index.open expects a bool for argument 'DEBUG'")

    if idx is None and ("x" not in mode and "w" not in mode):
        raise ValueError("Index data cannot be None if in create/write/overwrite mode")
    elif idx is not None and type(idx) is not list:
        raise TypeError("kmerdb.index.open expects a list for argument 'idx'")
    elif idx is not None and not all(type(i) is tuple for i in idx):
        raise TypeError("kmerdb.index.open expects a list of tuples for argument 'idx'")
        
    if "w" in mode and idx is not None and type(idx) is not np.ndarray:
        raise TypeError("kmerdb.index.open in create mode expects an int as its first positional argument")
    elif ("w" in mode or "x" in mode) and (kdbfile is None or not isinstance(kdbfile, fileutil.KDBReader)):
        raise ValueError("kmerdb.index.open could not create KDBReader")
            
    modes = set(list(mode.lower()))
    if modes - set("xrwbtu") or len(mode) > len(modes):
        raise ValueError("invalid mode: {}".format(mode))

    unknown = "u" in modes
    creating = "x" in modes
    reading  = "r" in modes
    writing  = "w" in modes
    binary   = "b" in modes
    text     = "t" in modes




    if "u" in mode.lower():
        if filepath is None or type(filepath) is not str:
            raise TypeError("kmerdb.index.open expects a str as its first positional argument")
        elif os.access(filepath, os.R_OK):
            logger.debug("Existing .kdb file found at '{0}'".format(filepath))
            
            if os.access(indexfilepath, os.R_OK):
                logger.debug("Additionally, index file (.kdbi) found at '{0}'".format(indexfilepath))
                if force is True:
                    logger.warning("--force flag provided. Forcing overwrite of index file!")
                    mode = "xt"
                else:
                    mode = "r"
                    raise ValueError("kmerdb.index.open: refusing to write over existing index data")

            else:
                logger.debug("No index data found. Assuming create mode")
                mode = "xt"
        elif not os.access(filepath, os.R_OK):
            raise ValueError("No .kdb data found for indexing.")
            
        else:
            raise ValueError("Could not determine mode from file data and invocation. Internal Error.")
    elif "w" in mode.lower():
        logger.error("Write mode is deprecated. Running in 'create' mode")
        logger.info("Externally provided index. Using write mode (deprecated)")
        raise RuntimeError("Internal error. Write mode is not supported")

        
    if text and binary:
        raise ValueError("can't have text and binary mode simultaneously")
    elif (creating or writing) and reading:
        raise ValueError("must have exactly one or read/write")

        
    if "r" in mode.lower():
        return IndexReader(filepath) # Will attempt to read this as an index file
    elif "x" in mode.lower() or set(list("xw")).intersection(modes):
        logger.info("Building index from file and constructing (mode='x') new .kdbi file")
        write_index(idx, indexfilepath, DEBUG=DEBUG) 



class IndexReader:
    def __init__(self, indexfile:str):
        if type(indexfile) is not str:
            raise TypeError("kmerdb.index.IndexReader expects a str as its first positional argument")
        elif not os.path.exists(indexfile):
            raise TypeError("kmerdb.index.IndexReader expects an existing .kdbi index file as its first positional argument")
        elif not is_gz_file(indexfile):
            raise IOError("kmerdb.index.IndexReader expects a gzip compressed index file as its first positional argument")
        
        idx = []

        #i = 0
        with gzip.open(indexfile, 'rt') as ifile:

            i = 0
            for line in ifile:
                j, kmer_id, line_offset = [int(i) for i in line.rstrip().split("\t")]
                idx.append(line_offset)
                i += 1
        self.index = np.array(idx, dtype="uint64")

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        del self.index
        return


def write_index(index:list, indexfile:str, DEBUG=False):
    if index is None or not(type(t) is tuple for t in index) or not(all(type(i) is int for i in t) for t in index):
        raise TypeError("kmerdb.index.write_index expects .kdb index data as its first positional argument")
    elif type(indexfile) is not str:
        raise TypeError("kmerdb.index.write_index expects a str as its second positional argument")
    elif os.path.exists(indexfile):
        logger.warning("kmerdb.index.write_index is overwriting an existing index '{0}'...".format(indexfile))
        mode = "wt"
    else:
        mode = "xt"

    # print(index)
    # raise RuntimeError("Printed index data")
        
    with gzip.open(indexfile, mode) as ofile:
        for i, offset_tuple in enumerate(index):
            kmer_id, block_offset = offset_tuple
            # File format is index, kmer_id, block_offset
            ofile.write("{0}\t{1}\t{2}\n".format(i, kmer_id, block_offset))

